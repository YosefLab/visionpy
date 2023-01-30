import math
from typing import Literal, Optional, Union

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from scipy.sparse import issparse
from scipy.stats import chisquare, pearsonr

from .diffexp import rank_genes_groups
from .signature import compute_obs_df_scores, compute_signature_scores


class AnnDataAccessor:
    def __init__(
        self,
        adata: Optional[anndata.AnnData] = None,
        norm_data_key: Optional[Union[Literal["use_raw"], str]] = None,
        protein_obsm_key: Optional[str] = None,
        signature_varm_key: Optional[str] = None,
        signature_names_uns_key: Optional[str] = None,
    ) -> None:
        """Accessor to AnnData object.

        Parameters
        ----------
        adata
            AnnData object
        norm_data_key
            Key for layer with log library size normalized data. If
            `None` (default), uses `adata.X`
        protein_obsm_key
            Location for protein data
        signature_varm_key
            Location for genes by signature matrix
        """
        self._adata = adata
        self._protein_obsm_key = protein_obsm_key
        self._norm_data_key = norm_data_key
        self._signature_varm_key = signature_varm_key
        self._signature_names_uns_key = signature_names_uns_key
        self._cells_selections = {}

    @property
    def adata(self):
        return self._adata

    @property
    def var_names(self):
        if self._norm_data_key == "use_raw":
            return self.adata.raw.var_names
        else:
            return self.adata.var_names

    @property
    def cells_selections(self):
        return self._cells_selections.keys()

    def add_cells_selection(self, key, val):
        self._cells_selections[key] = val

    def get_cells_selection(self, key):
        return self._cells_selections[key]

    @adata.setter
    def adata(self, adata: anndata.AnnData):
        self._adata = adata
        num_cols = adata.obs._get_numeric_data().columns.tolist()
        cols = adata.obs.columns.tolist()
        cat_vars = list(set(cols) - set(num_cols))
        self.cat_obs_cols = cat_vars
        self.numeric_obs_cols = num_cols

    @property
    def protein_obsm_key(self):
        return self._protein_obsm_key

    @protein_obsm_key.setter
    def protein_obsm_key(self, key: str):
        self._protein_obsm_key = key

    @property
    def signature_varm_key(self):
        return self._signature_varm_key

    @signature_varm_key.setter
    def signature_varm_key(self, key: str):
        self._signature_varm_key = key

    @property
    def signature_names_uns_key(self):
        return self._signature_names_uns_key

    @signature_names_uns_key.setter
    def signature_names_uns_key(self, key: str):
        self._signature_names_uns_key = key

    @property
    def norm_data_key(self):
        return self._norm_data_key

    @norm_data_key.setter
    def norm_data_key(self, key: str):
        self._norm_data_key = key

    def get_gene_expression(self, gene: str, return_list=True) -> list:
        if self.adata is None:
            raise ValueError("Accessor not populated with anndata.")
        if self.norm_data_key == "use_raw":
            data = self.adata.raw[:, gene].X
        elif self.norm_data_key is None:
            data = self.adata[:, gene].X
        else:
            data = self.adata[:, gene].layers[self.norm_data_key]

        if scipy.sparse.issparse(data):
            data = data.toarray()

        if return_list:
            return data.ravel().tolist()
        else:
            return data

    def get_genes_by_signature(self, sig_name: str) -> pd.DataFrame:
        """Df of genes in index, sign as values."""
        index = np.where(np.asarray(self.adata.uns[self.signature_names_uns_key]) == sig_name)[0][0]
        if self._norm_data_key == "use_raw":
            matrix = self.adata.raw.varm[self.signature_varm_key]
        else:
            matrix = self.adata.varm[self.signature_varm_key]

        if isinstance(matrix, pd.DataFrame):
            matrix = matrix.to_numpy()

        matrix = matrix[:, index]
        if issparse(matrix):
            matrix = matrix.toarray().ravel()

        mask = matrix != 0
        sig_df = pd.DataFrame(index=self.var_names[mask], data=matrix[mask])

        return sig_df

    def compute_obs_df_scores(self):
        self.adata.uns["vision_obs_df_scores"] = compute_obs_df_scores(self.adata)

    def compute_signature_scores(self):
        self.adata.uns["vision_signature_scores"] = compute_signature_scores(self.adata)

    def compute_one_vs_all_signatures(self):
        sig_adata = anndata.AnnData(self.adata.obsm["vision_signatures"])
        sig_adata.obs = self.adata.obs.loc[:, self.cat_obs_cols].copy()
        for c in self.cat_obs_cols:
            rank_genes_groups(
                sig_adata,
                groupby=c,
                key_added=f"rank_genes_groups_{c}",
                method="wilcoxon",
            )
        self.sig_adata = sig_adata

    def compute_one_vs_all_obs_cols(self):
        # log for scanpy de
        obs_adata = anndata.AnnData(np.log1p(self.adata.obs._get_numeric_data().copy()))
        obs_adata.obs = self.adata.obs.loc[:, self.cat_obs_cols].copy()
        for c in self.cat_obs_cols:
            try:
                rank_genes_groups(
                    obs_adata,
                    groupby=c,
                    key_added=f"rank_genes_groups_{c}",
                    method="wilcoxon",
                )
            # one category only has one obs
            except ValueError:
                # TODO: Log it
                self.cat_obs_cols = [c_ for c_ in self.cat_obs_cols if c_ != c]
                continue

            for g in categories(obs_adata.obs[c]):
                mask = (obs_adata.obs[c] == g).to_numpy()
                obs_pos_masked = obs_adata.obs.iloc[mask]
                obs_neg_masked = obs_adata.obs.iloc[~mask]
                for j in obs_pos_masked.columns:
                    pos_freq = obs_pos_masked[j].value_counts(normalize=False)
                    neg_freq = obs_neg_masked[j].value_counts(normalize=False)
                    freqs = pd.concat([pos_freq, neg_freq], axis=1).fillna(0)
                    # TODO: cramer's v might be incorrect
                    grand_total = np.sum(freqs.to_numpy())
                    r = len(freqs) - 1
                    try:
                        stat, pval = chisquare(
                            freqs.iloc[:, 0].to_numpy().ravel(),
                            freqs.iloc[:, 1].to_numpy().ravel(),
                        )
                    except ValueError:
                        stat = grand_total * r  # so that v is 1
                        pval = 0
                    if math.isinf(pval) or math.isnan(pval):
                        pval = 1
                    if math.isinf(stat) or math.isnan(stat):
                        v = 1
                    else:
                        v = np.sqrt(stat / (grand_total * r))
                    obs_adata.uns[f"chi_sq_{j}_{g}"] = {
                        "stat": v,
                        "pval": pval,
                    }

        self.obs_adata = obs_adata

    # TODO: refactor this function
    def compute_gene_score_per_signature(self):
        gene_score_sig = {}
        sig_names = self.adata.uns[self.signature_names_uns_key]
        for s in sig_names:
            gene_score_sig[s] = {"genes": [], "values": []}
            df = self.get_genes_by_signature(s)
            gene_names = df.index
            # cells by genes
            expr = np.array(self.get_gene_expression(gene_names, return_list=False))
            # cells
            sign = df.to_numpy().ravel()
            gene_score_sig[s]["signs"] = sign.tolist()

            # TODO: Make faster
            for i, (g, sign_) in enumerate(zip(gene_names, sign)):
                gene_score_sig[s]["values"].append(
                    sign_ * pearsonr(expr[:, i], self.adata.obsm["vision_signatures"][s])[0]
                )
                gene_score_sig[s]["genes"].append(g)

        for s in sig_names:
            info = gene_score_sig[s]
            gene_score_sig[s]["geneImportance"] = {g: v for g, v in zip(info["genes"], info["values"])}
            gene_score_sig[s]["sigDict"] = {g: v for g, v in zip(info["genes"], info["signs"])}

        self.gene_score_sig = gene_score_sig

    def compute_one_vs_one_de(self, key: str, group1: str, group2: str):
        rank_genes_groups(
            self.adata,
            groupby=key,
            groups=[group1],
            reference=group2,
            key_added=f"rank_genes_groups_{key}",
            method="wilcoxon",
            use_raw=self.norm_data_key == "use_raw",
            layer=self.norm_data_key if self.norm_data_key != "use_raw" else None,
        )
        return sc.get.rank_genes_groups_df(self.adata, group1, key=f"rank_genes_groups_{key}")


def categories(col):
    return col.astype("category").cat.categories
