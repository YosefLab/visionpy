import logging
import random
from re import compile, match
from typing import Dict, List, Literal, Optional, Sequence, Tuple, Union
import time
import numpy as np
import pandas as pd
from anndata import AnnData
from scanpy.metrics._gearys_c import _gearys_c
from .utils import _get_mean_var
from scipy import sparse
from scipy.sparse import csr_matrix, issparse
from scipy.stats import chi2_contingency
from sklearn.cluster import KMeans
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)

DOWN_SIG_KEY = "DN"
UP_SIG_KEY = "UP"


# ---------------------------------------------------------------------------
# Sparse-preserving, GPU-accelerated signature scoring (batchSigEvalNorm)
# ---------------------------------------------------------------------------

def _resolve_torch_device(device: str) -> Optional[str]:
    """Return a torch device string if a GPU path is available, else None.

    ``"auto"`` selects CUDA if available, then MPS (Apple Silicon), then
    falls back to ``None`` so the scipy path is used for pure-CPU runs.
    """
    try:
        import torch
    except ImportError:
        return None

    if device == "auto":
        if torch.cuda.is_available():
            return "cuda"
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            return "mps"
        return None  # no GPU → caller will use scipy path

    return device  # honour explicit override ("cpu", "cuda", "cuda:1", "mps", …)


def _sig_to_dense_f32(sig_matrix) -> np.ndarray:
    """Convert any signature-matrix format to a dense float32 ndarray."""
    if isinstance(sig_matrix, pd.DataFrame):
        return sig_matrix.values.astype(np.float32)
    if sparse.issparse(sig_matrix):
        return sig_matrix.toarray().astype(np.float32)
    return np.asarray(sig_matrix, dtype=np.float32)


def _sig_to_scipy_csc(sig_matrix) -> "sparse.csc_matrix":
    """Convert any signature-matrix format to scipy CSC float64."""
    if isinstance(sig_matrix, pd.DataFrame):
        return sparse.csc_matrix(sig_matrix.values)
    if sparse.issparse(sig_matrix):
        return sig_matrix.tocsc()
    return sparse.csc_matrix(np.asarray(sig_matrix, dtype=float))


def _batch_sig_eval_norm_scipy(norm_data, sig_matrix, batch_size: int) -> np.ndarray:
    """Scipy sparse CPU implementation of ``batch_sig_eval_norm``."""
    X = norm_data.data
    X = X.tocsr() if sparse.issparse(X) else sparse.csr_matrix(X)

    gof = norm_data.gene_offsets          # (n_genes,)
    gsf = norm_data.gene_scale_factors    # (n_genes,)
    cof = norm_data.cell_offsets          # (n_cells,)
    csf = norm_data.cell_scale_factors    # (n_cells,)

    S_all = _sig_to_scipy_csc(sig_matrix)  # (n_genes, n_sigs)
    n_sigs = S_all.shape[1]
    chunks = []

    for start in range(0, n_sigs, batch_size):
        S = S_all[:, start: start + batch_size]           # (n_genes, bs)

        S_gsf = S.multiply(gsf[:, None])                  # (n_genes, bs), sparse
        t1 = (X @ S_gsf).toarray()                         # (n_cells, bs), dense
        t2 = np.asarray(S_gsf.multiply(gof[:, None]).sum(axis=0)).ravel()  # (bs,)
        sig_sums = np.asarray(S.sum(axis=0)).ravel()      # (bs,)
        t3 = np.outer(cof, sig_sums)                      # (n_cells, bs)

        result = (t1 + t2[None, :] + t3) * csf[:, None]

        denom = np.asarray(np.abs(S).sum(axis=0)).ravel()
        denom = np.where(denom > 0, denom, 1.0)
        chunks.append(result / denom[None, :])

    return np.concatenate(chunks, axis=1)                  # (n_cells, n_sigs)


def _batch_sig_eval_norm_torch(
    norm_data, sig_matrix, device: str, batch_size: int
) -> np.ndarray:
    """PyTorch implementation of ``batch_sig_eval_norm`` (CUDA / MPS / CPU).

    CUDA:  expression matrix held as sparse CSR on the GPU;
           ``torch.sparse.mm`` drives the heavy sparse × dense multiply.
    MPS:   sparse tensors are not supported by Metal — X is densified and
           moved to the MPS device; use only when n_cells × n_genes fits
           in GPU VRAM.
    CPU:   falls through to the scipy path via the public dispatcher.
    """
    import torch

    dev = torch.device(device)
    use_cuda = dev.type == "cuda"
    use_mps = dev.type == "mps"

    # ── expression matrix ───────────────────────────────────────────────────
    X = norm_data.data
    if use_cuda:
        # Sparse CSR is well-supported on CUDA
        Xc = (X if sparse.issparse(X) else sparse.csr_matrix(X)).tocsr().astype(np.float32)
        X_pt = torch.sparse_csr_tensor(
            torch.from_numpy(Xc.indptr.copy().astype(np.int64)).to(dev),
            torch.from_numpy(Xc.indices.copy().astype(np.int64)).to(dev),
            torch.from_numpy(Xc.data.copy()).to(dev),
            size=tuple(Xc.shape),
        )
    else:
        # MPS / CPU: use a plain dense tensor
        X_dense = (X.toarray() if sparse.issparse(X) else np.asarray(X)).astype(np.float32)
        X_pt = torch.as_tensor(X_dense, device=dev)

    # ── normalisation parameters ────────────────────────────────────────────
    gof = torch.as_tensor(norm_data.gene_offsets.astype(np.float32), device=dev)
    gsf = torch.as_tensor(norm_data.gene_scale_factors.astype(np.float32), device=dev)
    cof = torch.as_tensor(norm_data.cell_offsets.astype(np.float32), device=dev)
    csf = torch.as_tensor(norm_data.cell_scale_factors.astype(np.float32), device=dev)

    # ── signature matrix as dense float32 on device ─────────────────────────
    # Signature matrices are typically small (n_genes × n_sigs) even when n_sigs
    # is large; dense is fine here; batching caps peak memory per iteration.
    S_full = torch.as_tensor(_sig_to_dense_f32(sig_matrix), device=dev)  # (n_genes, n_sigs)
    n_sigs = S_full.shape[1]

    chunks = []
    for start in range(0, n_sigs, batch_size):
        S = S_full[:, start: start + batch_size]           # (n_genes, bs)
        S_gsf = S * gsf.unsqueeze(1)                       # (n_genes, bs)

        # T1: (n_cells, n_genes) × (n_genes, bs) → (n_cells, bs)
        if use_cuda:
            t1 = torch.sparse.mm(X_pt, S_gsf)
        else:
            t1 = X_pt @ S_gsf

        # T2: per-sig constant broadcast over cells
        t2 = (S_gsf * gof.unsqueeze(1)).sum(dim=0)        # (bs,)
        # T3: per-(cell, sig) outer product
        t3 = torch.outer(cof, S.sum(dim=0))               # (n_cells, bs)

        result = (t1 + t2.unsqueeze(0) + t3) * csf.unsqueeze(1)
        denom = S.abs().sum(dim=0).clamp(min=1e-10)
        chunks.append((result / denom.unsqueeze(0)).cpu())

    return torch.cat(chunks, dim=1).numpy()                # (n_cells, n_sigs)


def batch_sig_eval_norm(
    norm_data,
    sig_matrix: Union[np.ndarray, "sparse.spmatrix", pd.DataFrame],
    device: str = "auto",
    batch_size: int = 1200,
) -> np.ndarray:
    """Evaluate signature scores via NormData without materialising the dense expression matrix.

    Mirrors R VISION's ``batchSigEvalNorm`` / ``innerEvalSignatureBatchNorm``.
    The lazy NormData normalisation is expanded analytically into three
    sparse-computable terms, so the full n_cells × n_genes matrix is never
    formed in memory::

        Z[c, g] = ((X[c,g] + gof[g]) * gsf[g] + cof[c]) * csf[c]
        score(s,c) = Σ_g sig[g,s] · Z[c,g]  /  Σ_g |sig[g,s]|

    Expansion (no dense n_cells × n_genes intermediate)::

        T1 = X @ (S * gsf)               sparse × dense
        T2 = Σ_g sig[g,s]*gsf[g]*gof[g]  per-sig scalar, broadcast
        T3 = cof ⊗ Σ_g sig[g,s]          outer product
        result = (T1 + T2 + T3) * csf[:,None] / denom

    Parameters
    ----------
    norm_data : NormData
        Lazy normalisation container produced by
        :func:`~visionpy.normalization.get_normalized_copy_sparse`.
    sig_matrix : array-like of shape (n_genes, n_sigs)
        Signature weight matrix (+1 up, −1 down, 0 absent).
    device : str, optional
        Compute device.  ``"auto"`` (default) selects CUDA when available,
        then MPS (Apple Silicon), then falls back to scipy on CPU.  Pass
        ``"cuda"``, ``"cuda:1"``, ``"mps"``, or ``"cpu"`` to override.
    batch_size : int, optional
        Signatures per batch (bounds peak GPU/RAM usage), by default 1200.

    Returns
    -------
    ndarray of shape (n_cells, n_sigs)
        Per-cell signature scores.

    Notes
    -----
    On MPS the expression matrix is densified before transfer to GPU.  For
    very large datasets (> 100 k cells × 30 k genes) prefer ``"cuda"`` or
    the default scipy path.
    """
    torch_device = _resolve_torch_device(device)

    if torch_device is not None:
        try:
            return _batch_sig_eval_norm_torch(norm_data, sig_matrix, torch_device, batch_size)
        except Exception as exc:
            logger.warning(
                "PyTorch path failed (%s); falling back to scipy sparse.", exc
            )

    return _batch_sig_eval_norm_scipy(norm_data, sig_matrix, batch_size)


# def compute_obs_df_scores(adata):
#     """Computes Geary's C for numerical data."""
#     numerical_df = adata.obs._get_numeric_data()
#     num_cols = numerical_df.columns.tolist()
#     cols = adata.obs.columns.tolist()
#     cat_cols = list(set(cols) - set(num_cols))

#     # first handle numerical data with geary's
#     # c = gearys_c(adata=adata, vals=numerical_df.to_numpy().transpose())
#     weights = adata.obsp["normalized_connectivities"]
#     gearys_c = _gearys_c(weights, numerical_df.to_numpy().transpose())
#     pvals = np.zeros_like(gearys_c)
#     fdr = np.zeros_like(gearys_c)

#     # categorical data
#     cramers_v = []
#     for c in cat_cols:
#         one_hot_cat = csr_matrix(pd.get_dummies(adata.obs[c]).to_numpy())
#         # cells by categories
#         c_hat_im = weights @ one_hot_cat
#         # categories by categories
#         x_lm = (c_hat_im.T @ (one_hot_cat)).toarray()

#         # Cramer's V
#         # https://www.statology.org/cramers-v-in-python/
#         chi_sq, pval, _, _ = chi2_contingency(x_lm)
#         n = np.sum(x_lm)
#         min_dim = min(x_lm.shape) - 1
#         cramers_v.append(np.sqrt((chi_sq / n) / min_dim))

#     # aggregate results
#     res = pd.DataFrame(index=adata.obs.columns)
#     res["c_prime"] = 0
#     res["pvals"] = 0
#     res["fdr"] = 0
#     res.loc[num_cols, "c_prime"] = 1 - gearys_c
#     res.loc[num_cols, "pvals"] = pvals
#     res.loc[num_cols, "fdr"] = fdr

#     res.loc[cat_cols, "c_prime"] = cramers_v

#     return res


def compute_obs_df_scores(adata):
    """Computes Geary's C for numerical data."""
    start = time.time()
    print("Computing observation scores...")

    numerical_df = adata.obs._get_numeric_data()
    num_cols = numerical_df.columns.tolist()
    cols = adata.obs.columns.tolist()
    cat_cols = list(set(cols) - set(num_cols))

    # aggregate results — float columns so pandas doesn't warn on partial assignment
    res = pd.DataFrame(0.0, index=adata.obs.columns, columns=["c_prime", "pvals", "fdr"])

    weights = adata.obsp["weights"].tocsr()

    # first handle numerical data with geary's
    if len(num_cols) > 0:
        gearys_c = _gearys_c(weights, numerical_df.to_numpy().transpose())
        c_prime = 1-gearys_c
        pvals_num = compute_num_var_pvals(c_prime, weights, numerical_df)
        fdr_num = multipletests(pvals_num, method="fdr_bh")[1]

        res.loc[num_cols, "c_prime"] = c_prime
        res.loc[num_cols, "pvals"] = pvals_num
        res.loc[num_cols, "fdr"] = fdr_num

    # categorical data
    if len(cat_cols) > 0:
        cramers_v = []
        pvals_cat = []
        for c in cat_cols:
            one_hot_cat = csr_matrix(pd.get_dummies(adata.obs[c]).to_numpy())
            # cells by categories
            c_hat_im = weights @ one_hot_cat
            # categories by categories
            x_lm = (c_hat_im.T @ (one_hot_cat)).toarray()

            # Cramer's V
            # https://www.statology.org/cramers-v-in-python/
            chi_sq, pval, _, _ = chi2_contingency(x_lm)
            n = np.sum(x_lm)
            min_dim = min(x_lm.shape) - 1
            cramers_v.append(np.sqrt((chi_sq / n) / min_dim))
            pvals_cat.append(pval)
        fdr_cat = multipletests(pvals_cat, method="fdr_bh")[1]

        res.loc[cat_cols, "c_prime"] = cramers_v
        res.loc[cat_cols, "pvals"] = pvals_cat
        res.loc[cat_cols, "fdr"] = fdr_cat

    print("Finished computing observation scores in %.3f seconds" %(time.time()-start))
    
    return res


def compute_signature_scores(adata, norm_data_key, signature_varm_key):
    """Computes Geary's C for all the signatures."""
    start = time.time()
    print("Computing Geary's C for all the signatures...")
    
    df = adata.obsm["vision_signatures"]
    # first handle numerical data with geary's
    weights = adata.obsp["weights"]
    weights = weights.tocsr()
    gearys_c = _gearys_c(weights, df.to_numpy().transpose())
    c_prime = 1-gearys_c

    print("Generating the null distribution...")

    random_sig_df, clusters, random_clusters = generate_permutations_null(adata, norm_data_key, signature_varm_key)
    adata.varm["random_signatures"] = random_sig_df

    random_df = compute_signatures_anndata(
        adata,
        norm_data_key=norm_data_key,
        signature_varm_key='random_signatures',
        signature_names_uns_key=None,
        )

    del adata.varm["random_signatures"]
    adata.obsm["vision_signatures"] = df  # restore real scores overwritten by random run

    random_gearys_c = _gearys_c(weights, random_df.to_numpy().transpose())
    random_c_prime = 1-random_gearys_c

    pvals = list(map(lambda i:compute_sig_pvals(c_prime, i, random_c_prime, clusters.tolist(), random_clusters), range(0, len(c_prime))))
    fdr = multipletests(pvals, method="fdr_bh")[1]

    # aggregate results
    res = pd.DataFrame(index=df.columns)
    res["c_prime"] = c_prime
    res["pvals"] = pvals
    res["fdr"] = fdr

    print("Finished computing Geary's C for all the signatures in %.3f seconds" %(time.time()-start))

    return res


# def compute_signature_scores(adata):
#     """Computes Geary's C for numerical data."""
#     df = adata.obsm["vision_signatures"]
#     # first handle numerical data with geary's
#     # c = gearys_c(adata=adata, vals=numerical_df.to_numpy().transpose())
#     weights = adata.obsp["normalized_connectivities"]
#     gearys_c = _gearys_c(weights, df.to_numpy().transpose())
#     np.zeros_like(gearys_c)
#     np.zeros_like(gearys_c)

#     # aggregate results
#     res = pd.DataFrame(index=df.columns)
#     res["c_prime"] = 1 - gearys_c
#     res["pvals"] = 0
#     res["fdr"] = 0

#     return res


def compute_num_var_pvals(c_prime, weights, numerical_df):

    num = 3000
    rand_c_prime_null = []

    for i in range(0, num):
        rand_numerical_df = numerical_df.sample(frac=1).copy()
        rand_gearys_c = _gearys_c(weights, rand_numerical_df.to_numpy().transpose())
        rand_c_prime = 1-rand_gearys_c
        rand_c_prime_null.append(rand_c_prime)

    rand_c_prime_null = np.stack(rand_c_prime_null, axis=0)
    x = np.sum(c_prime <= rand_c_prime_null, axis=0)
    p_values = (x+1) / (num+1)

    return p_values


def compute_sig_pvals(c_prime, i, random_c_prime, clusters, random_clusters):
    cluster = clusters[i]
    random_c_prime_clust = random_c_prime[random_clusters == cluster]
    x = (random_c_prime_clust >= c_prime[i]).sum()
    n = len(random_c_prime_clust)
    pval = (x+1)/(n+1)

    return pval


# def signatures_from_gmt(
#     gmt_files: Sequence[str],
#     adata: AnnData,
#     use_raw: bool = False,
# ) -> pd.DataFrame:
#     """
#     Compute signature scores from .gmt files.

#     Parameters
#     ----------
#     gmt_files
#         List of .gmt files to use for signature computation.
#     adata
#         AnnData object to compute signatures for.
#     use_raw
#         Whether to use adata.raw.X for signature computation.
#     layer
#         Layer to use for signature computation. Only used if
#         use_raw is False.
#     scale
#         Whether to scale the signature scores (using VISION-style
#         correction). Not used if use_decoupler is True.
#     use_decoupler
#         Use decouplerpy wmean for computing signature scores.
#         Otherwise, closed-form VISION computation is performed.

#     Returns
#     -------
#     Genes by signatures dataframe and cells by signatures dataframe
#     with scores. Index is aligned to genes from adata.
#     """
#     sig_dict = {}
#     for gmt_file in gmt_files:
#         sig_dict.update(read_gmt(gmt_file))
#     index = adata.raw.var.index if use_raw else adata.var_names
#     columns = list(sig_dict.keys())
#     data = np.zeros((len(index), len(columns)))
#     # Genes by signatures
#     sig_df = pd.DataFrame(index=index, columns=columns, data=data)
#     sig_df.index = sig_df.index.str.lower()
#     for sig_name, genes_up_down in sig_dict.items():
#         for key in [UP_SIG_KEY, DOWN_SIG_KEY]:
#             genes = genes_up_down.get(key, None)
#             if genes is not None:
#                 print(genes)
#                 genes = pd.Index(genes).str.lower()
#                 genes = genes.intersection(sig_df.index)
#                 sig_df.loc[genes, sig_name] = 1.0 if key == UP_SIG_KEY else -1.0
#     # Put back the old index
#     sig_df.index = index
#     sig_df.columns = sig_df.columns.str.lower()

#     return sig_df


# def read_gmt(
#     gmt_file: str,
#     up_suffix: Tuple[str] = "(_up|_plus)",
#     down_suffix: str = "(_down|_dn|_minus|_mius)",
#     verbose: bool = False,
# ) -> Dict[str, Dict[str, List[str]]]:
#     """Read gmt file to extract signed genesets.

#     Code adapted from https://bedapub.github.io/besca/sig/besca.tl.sig.read_GMT_sign.html
#     """
#     with open(gmt_file) as f:
#         text_gmt = f.read().split("\n")
#     signed_sign = {}
#     # Here \S is used as signature might have '-' in their name
#     #  (\w is not suficient if number in signature for EX.)
#     pattern_down = compile(r"(\S+)" + down_suffix + "$")
#     pattern_up = compile(r"(\S+)" + up_suffix + "$")
#     # TODO: remove this for loop.
#     for i in range(len(text_gmt)):
#         temp_split = text_gmt[i].split("\t")
#         signature_full_name = temp_split[0]
#         if len(temp_split) < 3:
#             if verbose:
#                 print("Skipping empty entry; less than 3 fields for " + signature_full_name)
#             continue
#         # Skipping empty lines in gmt files
#         if len(signature_full_name):
#             z = match(pattern_down, signature_full_name.lower())
#             if z:
#                 signature_name = z.groups()[0]
#                 direction = DOWN_SIG_KEY
#             else:
#                 z = match(pattern_up, signature_full_name.lower())
#                 if z:
#                     signature_name = z.groups()[0]
#                     direction = UP_SIG_KEY
#                 else:
#                     signature_name = signature_full_name
#                     direction = UP_SIG_KEY
#             # Get the gene names removing empty entry
#             initial_value = temp_split[2 : len(temp_split)]
#             gene_list = [x for x in initial_value if len(x)]
#             if signature_name in signed_sign.keys():
#                 signed_sign[signature_name][direction] = gene_list
#             else:
#                 signed_sign[signature_name] = {direction: gene_list}
#             if verbose:
#                 print(i, ": ", signature_full_name)

#     return signed_sign


# def compute_signatures_anndata(
#     adata: AnnData,
#     norm_data_key: Union[Literal["use_raw"], str],
#     signature_varm_key: str,
#     signature_names_uns_key: str,
# ) -> None:
#     """
#     Compute signatures for each cell.

#     Parameters
#     ----------
#     adata
#         AnnData object to compute signatures for.
#     norm_data_key
#         Key for adata.layers to use for signature computation. If "use_raw", use adata.raw.X.
#     signature_varm_key
#         Key in `adata.varm` for signatures. If `None` (default), no signatures. Matrix
#         should encode positive genes with 1, negative genes with -1, and all other genes with 0
#     signature_names_uns_key
#         Key in `adata.uns` for signature names. If `None`, attempts to read columns if `signature_varm_key`
#         is a pandas DataFrame. Otherwise, uses `Signature_1`, `Signature_2`, etc.
#     """
#     use_raw_for_signatures = norm_data_key == "use_raw"
#     if norm_data_key is None:
#         gene_expr = adata.X
#     elif norm_data_key == "use_raw":
#         gene_expr = adata.raw.X
#     else:
#         gene_expr = adata.layers[norm_data_key]
#     sig_matrix = adata.varm[signature_varm_key] if not use_raw_for_signatures else adata.raw.varm[signature_varm_key]
#     if signature_names_uns_key is not None:
#         cols = adata.uns[signature_names_uns_key]
#     elif isinstance(sig_matrix, pd.DataFrame):
#         cols = sig_matrix.columns
#     else:
#         cols = [f"signature_{i}" for i in range(sig_matrix.shape[1])]
#     if not issparse(sig_matrix):
#         if isinstance(sig_matrix, pd.DataFrame):
#             sig_matrix = sig_matrix.to_numpy()
#         sig_matrix = csr_matrix(sig_matrix)
#     cell_signature_matrix = (gene_expr @ sig_matrix).toarray()
#     sig_df = pd.DataFrame(data=cell_signature_matrix, columns=cols, index=adata.obs_names)
#     # normalize
#     mean, var = _get_mean_var(gene_expr, axis=1)
#     n = np.asarray((sig_matrix > 0).sum(0))
#     m = np.asarray((sig_matrix < 0).sum(0))
#     sig_df = sig_df / (n + m)
#     # cells by signatures
#     sig_mean = np.outer(mean, ((n - m) / (n + m)))
#     # cells by signatures
#     sig_std = np.sqrt(np.outer(var, 1 / (n + m)))
#     sig_df = (sig_df - sig_mean) / sig_std

#     adata.obsm["vision_signatures"] = sig_df
    
    
def signatures_from_file(
    adata: AnnData,
    use_raw: bool = False,
    species: Optional[Union[Literal["mouse"], Literal["human"]]] = None,
    gmt_files: Optional[Sequence[str]] = None,
    dicts: Optional[Sequence[dict]] = None,
    sig_names: Optional[Sequence[str]] = None,
):
    """Compute signature scores from .gmt files.

    Parameters
    ----------
    adata
        AnnData object to compute signatures for.
    use_raw
        Whether to use adata.raw.X for signature computation.
    species
        Species identity to select one (or more) of the signatures downloaded from MSigDB.
    gmt_files
        List of .gmt files to use for signature computation.
    dicts
        List of dictionaries to use for signature computation.
    sig_names
        List of signature file names downloaded from MSigDB to use for signature computation.

    Returns
    -------
    Genes by signatures dataframe and cells by signatures dataframe
    with scores. Index is aligned to genes from adata.

    """
    signature_paths = {
        'human': {
            'BioCarta_Human': 'https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.2.Hs/c2.cp.biocarta.v2023.2.Hs.symbols.gmt',
            'Reactome_Human': 'https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.2.Hs/c2.cp.reactome.v2023.2.Hs.symbols.gmt',
            'KEGG_Human': 'https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.2.Hs/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt',
            'GO:BP_Human': 'https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.2.Hs/c5.go.bp.v2023.2.Hs.symbols.gmt',
        },
        'mouse': {
            'BioCarta_Mouse': 'https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.2.Mm/m2.cp.biocarta.v2023.2.Mm.symbols.gmt',
            'Reactome_Mouse': 'https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.2.Mm/m2.cp.reactome.v2023.2.Mm.symbols.gmt',
            'GO:BP_Mouse': 'https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.2.Mm/m5.go.bp.v2023.2.Mm.symbols.gmt',
        }
    }

    if dicts is None and species not in ['human', 'mouse'] and gmt_files is None:
        raise ValueError(f'species type: {species} is not supported currently. You should choose either "human" or "mouse".')

    if dicts is None and gmt_files is None and (sig_names is None or sig_names not in list(signature_paths[species].keys())):
        raise ValueError(f'Please provide either your own signature list or select one (or more) of: {list(signature_paths[species].keys())}')

    if dicts is not None and gmt_files is not None:
        files = dicts + gmt_files
    elif dicts is None and gmt_files is not None:
        files = gmt_files
    elif dicts is not None and gmt_files is None:
        files = dicts
    
    if type(sig_names) is list:
        if type(gmt_files) is list:
            for sig_name in sig_names:
                if sig_name not in list(signature_paths[species].keys()):
                    raise ValueError(f'{sig_name} is not available. Please select one (or more) of: {list(signature_paths[species].keys())}')
                else:
                    sig_name_path = signature_paths[species][sig_name]
                    files.append(sig_name_path)
        else:
            files = []
            for sig_name in sig_names:
                sig_name_path = signature_paths[species][sig_name]
                files.append(sig_name_path)

    sig_dict = {}
    for file in files:
        if isinstance(file, dict):
            sig_dict.update(read_dict(file))
        else:
            sig_dict.update(read_gmt(file))
    index = adata.raw.var.index if use_raw else adata.var_names
    columns = list(sig_dict.keys())
    data = np.zeros((len(index), len(columns)))
    sig_df = pd.DataFrame(index=index, columns=columns, data=data)
    sig_df.index = sig_df.index.str.lower()
    for sig_name, genes_up_down in sig_dict.items():
        for key in [UP_SIG_KEY, DOWN_SIG_KEY]:
            genes = genes_up_down.get(key, None)
            if genes is not None:
                genes = pd.Index(genes).str.lower()
                genes = genes.intersection(sig_df.index)
                sig_df.loc[genes, sig_name] = 1.0 if key == UP_SIG_KEY else -1.0
    sig_df.index = index
    sig_df = sig_df.loc[:, (sig_df!=0).any(axis=0)]

    adata.varm["signatures"] = sig_df

    return


def read_gmt(
    gmt_file: str,
    up_suffix: Tuple[str] = "(_up|_plus)",
    down_suffix: str = "(_down|_dn|_minus|_mius)",
    verbose: bool = False,
) -> Dict[str, Dict[str, List[str]]]:
    """Read gmt file to extract signed genesets.

    Code adapted from https://bedapub.github.io/besca/sig/besca.tl.sig.read_GMT_sign.html
    """
    with open(gmt_file) as f:
        text_gmt = f.read().split("\n")
    signed_sign = {}
    # Here \S is used as signature might have '-' in their name
    #  (\w is not suficient if number in signature for EX.)
    pattern_down = compile(r"(\S+)" + down_suffix + "$")
    pattern_up = compile(r"(\S+)" + up_suffix + "$")
    # TODO: remove this for loop.
    for i in range(len(text_gmt)):
        temp_split = text_gmt[i].split("\t")
        signature_full_name = temp_split[0]
        if len(temp_split) < 3:
            if verbose:
                print("Skipping empty entry; less than 3 fields for " + signature_full_name)
            continue
        # Skipping empty lines in gmt files
        if len(signature_full_name):
            z = match(pattern_down, signature_full_name.lower())
            if z:
                signature_name = z.groups()[0]
                direction = DOWN_SIG_KEY
            else:
                z = match(pattern_up, signature_full_name.lower())
                if z:
                    signature_name = z.groups()[0]
                    direction = UP_SIG_KEY
                else:
                    signature_name = signature_full_name
                    direction = UP_SIG_KEY
            # Get the gene names removing empty entry
            initial_value = temp_split[2 : len(temp_split)]
            gene_list = [x for x in initial_value if len(x)]
            if signature_name in signed_sign.keys():
                signed_sign[signature_name][direction] = gene_list
            else:
                signed_sign[signature_name] = {direction: gene_list}
            if verbose:
                print(i, ": ", signature_full_name)

    return signed_sign


def read_dict(
    dict: dict,
    up_suffix: Tuple[str] = "(_up|_plus)",
    down_suffix: str = "(_down|_dn|_minus|_mius)",
) -> Dict[str, Dict[str, List[str]]]:
    """Read dictionary to extract signed genesets.
    """
    signed_sign = {}

    pattern_down = compile(r"(\S+)" + down_suffix + "$")
    pattern_up = compile(r"(\S+)" + up_suffix + "$")
    
    for signature_full_name, gene_list in dict.items():
        
        z = match(pattern_down, signature_full_name.lower())
        if z:
            signature_name = z.groups()[0]
            direction = DOWN_SIG_KEY
        else:
            z = match(pattern_up, signature_full_name.lower())
            if z:
                signature_name = z.groups()[0]
                direction = UP_SIG_KEY
            else:
                signature_name = signature_full_name
                direction = UP_SIG_KEY

        if signature_name in signed_sign.keys():
            signed_sign[signature_name][direction] = gene_list
        else:
            signed_sign[signature_name] = {direction: gene_list}

    return signed_sign


def compute_signatures_anndata(
    adata: AnnData,
    norm_data_key: Union[Literal["use_raw"], str],
    signature_varm_key: str,
    signature_names_uns_key: Optional[str] = None,
    sig_norm_method: str = "znorm_columns",
    device: str = "auto",
    batch_size: int = 1200,
) -> pd.DataFrame:
    """Compute per-cell signature scores.

    For sparse expression matrices the computation is performed entirely
    without densifying the n_cells × n_genes matrix.  A
    :class:`~visionpy.normalization.NormData` object is built from the
    expression data and :func:`batch_sig_eval_norm` is called, selecting
    CUDA → MPS → scipy automatically.

    For dense inputs the existing VISION z-score formula is retained.

    Parameters
    ----------
    adata
        AnnData object.
    norm_data_key
        Key in ``adata.layers`` for the expression matrix, or ``"use_raw"``
        to read from ``adata.raw.X``.  Pass ``None`` to use ``adata.X``.
    signature_varm_key
        Key in ``adata.varm`` for the (n_genes × n_sigs) signature weight
        matrix (+1 up-regulated, −1 down-regulated, 0 absent).
    signature_names_uns_key
        Key in ``adata.uns`` for column names.  When ``None``, column names
        are read from the DataFrame index (if applicable) or auto-generated.
    sig_norm_method : str, optional
        Normalisation method passed to
        :func:`~visionpy.normalization.get_normalized_copy_sparse` for the
        sparse path.  Matches R VISION's ``sig_norm_method`` parameter;
        default ``"znorm_columns"`` (z-normalise each cell across genes).
    device : str, optional
        Compute device for :func:`batch_sig_eval_norm`; ``"auto"`` by default.
    batch_size : int, optional
        Signatures processed per batch, by default 1200.
    """
    from .normalization import get_normalized_copy_sparse

    t0 = time.time()
    logger.info("Computing signatures for each cell...")

    adata.uns["norm_data_key"] = norm_data_key
    adata.uns["signature_varm_key"] = signature_varm_key

    use_raw = norm_data_key == "use_raw"
    if norm_data_key is None:
        gene_expr = adata.X
    elif use_raw:
        gene_expr = adata.raw.X
    else:
        gene_expr = adata.layers[norm_data_key]

    sig_mat_raw = adata.varm[signature_varm_key] if not use_raw else adata.raw.varm[signature_varm_key]
    if signature_names_uns_key is not None:
        cols = adata.uns[signature_names_uns_key]
    elif isinstance(sig_mat_raw, pd.DataFrame):
        cols = sig_mat_raw.columns.tolist()
    else:
        cols = [f"signature_{i}" for i in range(sig_mat_raw.shape[1])]

    if issparse(gene_expr):
        # ── sparse path: NormData + GPU-accelerated batch scoring ──────────
        norm_data = get_normalized_copy_sparse(gene_expr, method=sig_norm_method)
        scores = batch_sig_eval_norm(
            norm_data, sig_mat_raw, device=device, batch_size=batch_size
        )
        sig_df = pd.DataFrame(data=scores, columns=cols, index=adata.obs_names)
    else:
        # ── dense path: existing VISION z-score formula ─────────────────────
        gene_expr = np.asarray(gene_expr, dtype=float)
        sig_matrix = csr_matrix(
            sig_mat_raw.to_numpy() if isinstance(sig_mat_raw, pd.DataFrame)
            else (sig_mat_raw.toarray() if issparse(sig_mat_raw) else sig_mat_raw)
        )
        # Avoid densifying: sparse @ sparse → small dense (n_cells × n_sigs)
        cell_signature_matrix = (gene_expr @ sig_matrix).toarray() \
            if issparse(gene_expr @ sig_matrix) \
            else gene_expr @ sig_matrix
        sig_df = pd.DataFrame(data=cell_signature_matrix, columns=cols, index=adata.obs_names)
        mean, var = _get_mean_var(gene_expr, axis=1)
        n = np.asarray((sig_matrix > 0).sum(0))
        m = np.asarray((sig_matrix < 0).sum(0))
        denom = (n + m)
        sig_df = sig_df / np.where(denom > 0, denom, 1.0)
        sig_mean = np.outer(mean, np.where(denom > 0, (n - m) / denom, 0.0))
        sig_std = np.sqrt(np.outer(var, np.where(denom > 0, 1.0 / denom, 1.0)))
        sig_std[sig_std == 0] = 1.0
        sig_df = (sig_df - sig_mean) / sig_std

    adata.obsm["vision_signatures"] = sig_df
    logger.info("Signatures computed in %.3f s.", time.time() - t0)
    return sig_df


def generate_permutations_null(adata, norm_data_key, signature_varm_key):

    use_raw = norm_data_key == "use_raw"
    exp_genes = adata.raw.var.index if use_raw else adata.var_names

    sig_matrix = adata.varm[signature_varm_key] if not use_raw else adata.raw.varm[signature_varm_key]
    if not isinstance(sig_matrix, pd.DataFrame):
        sig_matrix = pd.DataFrame(sig_matrix)

    sig_size = sig_matrix.apply(lambda x: (x!=0).sum(), axis=0)

    sig_size = np.log10(sig_size)
    sig_balance = sig_matrix.apply(lambda x: (x>0).sum()/(x!=0).sum(), axis=0)

    # sig_balance[sig_balance < 0.5] = 1 - sig_balance[sig_balance < 0.5]

    sig_size.name = 'sig_size'
    sig_balance.name = 'sig_balance'
    sig_vars = pd.concat([sig_size, sig_balance],axis=1)

    n_components = 5

    if sig_vars.shape[0] <= n_components:
        n_components = sig_vars.shape[0]
        centers = sig_vars
        clusters = pd.Series(list(range(0, sig_vars.shape[0])))
        clusters = clusters.astype('category')
        clusters.index = sig_vars.index
    else:
        if sig_vars.drop_duplicates().shape[0] <= n_components:
            n_components = sig_vars.drop_duplicates().shape[0]

        kmeans = KMeans(init="random", n_clusters=n_components, random_state=42)
        kmeans.fit(sig_vars)

        centers = kmeans.cluster_centers_
        centers = pd.DataFrame(centers, columns=sig_vars.columns)

        levels = list(range(0, n_components))
        levels = [str(x) for x in levels]
        clusters = kmeans.labels_

    # # Re-order the centers
    # row_i = centers.sort_values("sig_size").index

    # centers = centers.reindex(row_i)
    # # levels(clusters) = as.character(order(row_i))
    # centers.index = levels

    # undo the log scaling
    centers["sig_size"] = round(10 ** centers["sig_size"])

    num = 3000
    gene_list = exp_genes.tolist()
    n_genes = len(gene_list)
    gene_to_idx = {g: i for i, g in enumerate(gene_list)}

    random_clusters = []
    new_sig_names = []
    # Pre-allocate the full matrix to avoid per-column fragmentation
    n_total = num * len(centers.index)
    mat = np.zeros((n_genes, n_total), dtype=np.float32)

    col_idx = 0
    for cluster_i in centers.index:
        size = centers.loc[cluster_i, "sig_size"]
        balance = centers.loc[cluster_i, "sig_balance"]

        for j in range(0, num):
            new_sig_genes = random.sample(gene_list, int(min(size, n_genes)))

            up_genes = round(balance * size)
            remainder = (balance * size) % 1
            if np.random.uniform(low=0, high=1, size=1) < remainder:
                up_genes = up_genes + 1
            new_sig_signs = [1] * up_genes + [-1] * int(size - up_genes)

            gene_idxs = [gene_to_idx[g] for g in new_sig_genes]
            mat[gene_idxs, col_idx] = new_sig_signs
            new_sig_names.append("RANDOM_BG_" + str(cluster_i) + "_" + str(j))
            random_clusters.append(cluster_i)
            col_idx += 1

    random_sig_df = pd.DataFrame(mat, index=exp_genes, columns=new_sig_names)

    random_clusters = pd.Series(random_clusters, index=new_sig_names, dtype="category")

    return random_sig_df, clusters, random_clusters
