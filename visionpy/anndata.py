from typing import Optional, Union

import anndata
import scipy
from scipy.sparse import issparse
from scipy.stats import pearsonr
import scanpy as sc
import pandas as pd
import numpy as np
from ._compat import Literal
from .signature import compute_obs_df_scores, compute_signature_scores


class AnnDataAccessor(object):
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

    @property
    def adata(self):
        return self._adata

    @property
    def var_names(self):
        if self._norm_data_key == "use_raw":
            return self.adata.raw.var_names
        else:
            return self.adata.var_names

    @adata.setter
    def adata(self, adata: anndata.AnnData):
        self._adata = adata
        num_cols = adata.obs._get_numeric_data().columns.tolist()
        cols = adata.obs.columns.tolist()
        cat_vars = list(set(cols) - set(num_cols))
        self.cat_obs_cols = cat_vars

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
        index = np.where(
            np.asarray(self.adata.uns[self.signature_names_uns_key]) == sig_name
        )[0][0]
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
            sc.tl.rank_genes_groups(
                sig_adata,
                groupby=c,
                key_added="rank_genes_groups_{}".format(c),
                method="wilcoxon",
            )
        self.sig_adata = sig_adata

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
                    sign_
                    * pearsonr(expr[:, i], self.adata.obsm["vision_signatures"][s])[0]
                )
                gene_score_sig[s]["genes"].append(g)

        for s in sig_names:
            info = gene_score_sig[s]
            gene_score_sig[s]["geneImportance"] = {
                g: v for g, v in zip(info["genes"], info["values"])
            }
            gene_score_sig[s]["sigDict"] = {
                g: v for g, v in zip(info["genes"], info["signs"])
            }

        self.gene_score_sig = gene_score_sig
