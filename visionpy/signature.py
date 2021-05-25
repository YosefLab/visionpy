import numpy as np
import pandas as pd
from scanpy.metrics._gearys_c import _gearys_c
from scipy.sparse import csr_matrix
from scipy.stats import chi2_contingency


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
    pvals = np.zeros_like(gearys_c)
    fdr = np.zeros_like(gearys_c)

    # aggregate results
    res = pd.DataFrame(index=df.columns)
    res["c_prime"] = 1 - gearys_c
    res["pvals"] = 0
    res["fdr"] = 0

    return res
