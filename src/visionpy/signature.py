import random
from re import compile, match
from typing import Dict, List, Literal, Optional, Sequence, Tuple, Union
import time
import numpy as np
import pandas as pd
from anndata import AnnData
from scanpy.metrics._gearys_c import _gearys_c
from scanpy.preprocessing._utils import _get_mean_var
from scipy.sparse import csr_matrix, issparse
from scipy.stats import chi2_contingency
from sklearn.cluster import KMeans
from statsmodels.stats.multitest import multipletests

DOWN_SIG_KEY = "DN"
UP_SIG_KEY = "UP"


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

    # aggregate results
    res = pd.DataFrame(index=adata.obs.columns)
    res["c_prime"] = 0
    res["pvals"] = 0
    res["fdr"] = 0

    # first handle numerical data with geary's
    if len(num_cols) > 0:
        weights = adata.obsp["weights"]
        weights = weights.tocsr()
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
) -> None:
    """Compute signatures for each cell.

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
    start = time.time()
    print("Computing signatures for each cell...")
    
    adata.uns['norm_data_key'] = norm_data_key
    adata.uns['signature_varm_key'] = signature_varm_key

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
    is_sparse = issparse(gene_expr)
    # gene_expr = gene_expr.A if is_sparse else gene_expr
    gene_expr = gene_expr.toarray() if is_sparse else gene_expr
    cell_signature_matrix = gene_expr @ sig_matrix
    # cell_signature_matrix = (gene_expr @ sig_matrix).toarray()
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

    print("Finished computing signatures for each cell in %.3f seconds" %(time.time()-start))

    return sig_df


def generate_permutations_null(adata, norm_data_key, signature_varm_key):

    use_raw = norm_data_key == "use_raw"
    exp_genes = adata.raw.var.index if use_raw else adata.var_names

    sig_matrix = adata.varm[signature_varm_key] if not use_raw else adata.raw.varm[signature_varm_key]

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

    random_sig_df = pd.DataFrame(index=exp_genes)

    num = 3000
    random_clusters = []
    new_sig_names = []

    for cluster_i in centers.index:

        size = centers.loc[cluster_i, "sig_size"]
        balance = centers.loc[cluster_i, "sig_balance"]

        for j in range(0, num):

            new_sig_genes = random.sample(exp_genes.tolist(), int(min(size, len(exp_genes))))

            up_genes = round(balance * size)
            remainder = (balance * size) % 1
            if np.random.uniform(low=0, high=1, size=1) < remainder:
                up_genes = up_genes + 1
            new_sig_signs = [1]*up_genes + [-1]*int(size - up_genes)
            new_sig_name = "RANDOM_BG_" + str(cluster_i) + "_" + str(j)
            new_sig_names.append(new_sig_name)

            random_sig_df[new_sig_name] = 0
            random_sig_df.loc[new_sig_genes, new_sig_name] = new_sig_signs

            random_clusters.append(cluster_i)

    random_clusters = pd.Series(random_clusters)
    random_clusters = random_clusters.astype('category')
    random_clusters.index = new_sig_names

    return random_sig_df, clusters, random_clusters
