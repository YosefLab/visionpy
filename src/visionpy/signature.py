from re import compile, match
from typing import Dict, List, Literal, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scanpy.metrics._gearys_c import _gearys_c
from scanpy.preprocessing._utils import _get_mean_var
from scipy.sparse import csr_matrix, issparse
from scipy.stats import chi2_contingency

DOWN_SIG_KEY = "DN"
UP_SIG_KEY = "UP"


def compute_obs_df_scores(adata):
    """Computes Geary's C for numerical data."""
    numerical_df = adata.obs._get_numeric_data()
    num_cols = numerical_df.columns.tolist()
    cols = adata.obs.columns.tolist()
    cat_cols = list(set(cols) - set(num_cols))

    # first handle numerical data with geary's
    # c = gearys_c(adata=adata, vals=numerical_df.to_numpy().transpose())
    weights = adata.obsp["normalized_connectivities"]
    gearys_c = _gearys_c(weights, numerical_df.to_numpy().transpose())
    pvals = np.zeros_like(gearys_c)
    fdr = np.zeros_like(gearys_c)

    # categorical data
    cramers_v = []
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

    # aggregate results
    res = pd.DataFrame(index=adata.obs.columns)
    res["c_prime"] = 0
    res["pvals"] = 0
    res["fdr"] = 0
    res.loc[num_cols, "c_prime"] = 1 - gearys_c
    res.loc[num_cols, "pvals"] = pvals
    res.loc[num_cols, "fdr"] = fdr

    res.loc[cat_cols, "c_prime"] = cramers_v

    return res


def compute_signature_scores(adata):
    """Computes Geary's C for numerical data."""
    df = adata.obsm["vision_signatures"]
    # first handle numerical data with geary's
    # c = gearys_c(adata=adata, vals=numerical_df.to_numpy().transpose())
    weights = adata.obsp["normalized_connectivities"]
    gearys_c = _gearys_c(weights, df.to_numpy().transpose())
    np.zeros_like(gearys_c)
    np.zeros_like(gearys_c)

    # aggregate results
    res = pd.DataFrame(index=df.columns)
    res["c_prime"] = 1 - gearys_c
    res["pvals"] = 0
    res["fdr"] = 0

    return res


def signatures_from_gmt(
    gmt_files: Sequence[str],
    adata: AnnData,
    use_raw: bool = False,
) -> pd.DataFrame:
    """
    Compute signature scores from .gmt files.

    Parameters
    ----------
    gmt_files
        List of .gmt files to use for signature computation.
    adata
        AnnData object to compute signatures for.
    use_raw
        Whether to use adata.raw.X for signature computation.
    layer
        Layer to use for signature computation. Only used if
        use_raw is False.
    scale
        Whether to scale the signature scores (using VISION-style
        correction). Not used if use_decoupler is True.
    use_decoupler
        Use decouplerpy wmean for computing signature scores.
        Otherwise, closed-form VISION computation is performed.

    Returns
    -------
    Genes by signatures dataframe and cells by signatures dataframe
    with scores. Index is aligned to genes from adata.
    """
    sig_dict = {}
    for gmt_file in gmt_files:
        sig_dict.update(read_gmt(gmt_file))
    index = adata.raw.var.index if use_raw else adata.var_names
    columns = list(sig_dict.keys())
    data = np.zeros((len(index), len(columns)))
    # Genes by signatures
    sig_df = pd.DataFrame(index=index, columns=columns, data=data)
    sig_df.index = sig_df.index.str.lower()
    for sig_name, genes_up_down in sig_dict.items():
        for key in [UP_SIG_KEY, DOWN_SIG_KEY]:
            genes = genes_up_down.get(key, None)
            if genes is not None:
                print(genes)
                genes = pd.Index(genes).str.lower()
                genes = genes.intersection(sig_df.index)
                sig_df.loc[genes, sig_name] = 1.0 if key == UP_SIG_KEY else -1.0
    # Put back the old index
    sig_df.index = index
    sig_df.columns = sig_df.columns.str.lower()

    return sig_df


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


def compute_signatures_anndata(
    adata: AnnData,
    norm_data_key: Union[Literal["use_raw"], str],
    signature_varm_key: str,
    signature_names_uns_key: str,
) -> None:
    """
    Compute signatures for each cell.

    Parameters
    ----------
    adata
        AnnData object to compute signatures for.
    norm_data_key
        Key for adata.layers to use for signature computation. If "use_raw", use adata.raw.X.
    signature_varm_key
        Key in `adata.varm` for signatures. If `None` (default), no signatures. Matrix
        should encode positive genes with 1, negative genes with -1, and all other genes with 0
    signature_names_uns_key
        Key in `adata.uns` for signature names. If `None`, attempts to read columns if `signature_varm_key`
        is a pandas DataFrame. Otherwise, uses `Signature_1`, `Signature_2`, etc.
    """
    use_raw_for_signatures = norm_data_key == "use_raw"
    if norm_data_key is None:
        gene_expr = adata.X
    elif norm_data_key == "use_raw":
        gene_expr = adata.raw.X
    else:
        gene_expr = adata.layers[norm_data_key]
    sig_matrix = adata.varm[signature_varm_key] if not use_raw_for_signatures else adata.raw.varm[signature_varm_key]
    if signature_names_uns_key is not None:
        cols = adata.uns[signature_names_uns_key]
    elif isinstance(sig_matrix, pd.DataFrame):
        cols = sig_matrix.columns
    else:
        cols = [f"signature_{i}" for i in range(sig_matrix.shape[1])]
    if not issparse(sig_matrix):
        if isinstance(sig_matrix, pd.DataFrame):
            sig_matrix = sig_matrix.to_numpy()
        sig_matrix = csr_matrix(sig_matrix)
    cell_signature_matrix = (gene_expr @ sig_matrix).toarray()
    sig_df = pd.DataFrame(data=cell_signature_matrix, columns=cols, index=adata.obs_names)
    # normalize
    mean, var = _get_mean_var(gene_expr, axis=1)
    n = np.asarray((sig_matrix > 0).sum(0))
    m = np.asarray((sig_matrix < 0).sum(0))
    sig_df = sig_df / (n + m)
    # cells by signatures
    sig_mean = np.outer(mean, ((n - m) / (n + m)))
    # cells by signatures
    sig_std = np.sqrt(np.outer(var, 1 / (n + m)))
    sig_df = (sig_df - sig_mean) / sig_std

    adata.obsm["vision_signatures"] = sig_df
