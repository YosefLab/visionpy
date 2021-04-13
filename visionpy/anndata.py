from typing import Optional, Union

import anndata
import scipy
import scanpy as sc

from ._compat import Literal
from .signature import compute_obs_df_scores, compute_signature_scores


class AnnDataAccessor(object):
    def __init__(
        self,
        adata: Optional[anndata.AnnData] = None,
        norm_data_key: Optional[Union[Literal["use_raw"], str]] = None,
        protein_obsm_key: Optional[str] = None,
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
        """
        self._adata = adata
        self._protein_obsm_key = protein_obsm_key
        self._norm_data_key = norm_data_key

    @property
    def adata(self):
        return self._adata

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
    def norm_data_key(self):
        return self._norm_data_key

    @norm_data_key.setter
    def norm_data_key(self, key: str):
        self._norm_data_key = key

    def get_gene_expression(self, gene: str) -> list:
        if self.adata is None:
            raise ValueError("Accessor not populated with anndata.")
        if self.norm_data_key == "use_raw":
            data = self.adata.raw[:, gene].X
        elif self.norm_data_key is None:
            data = self.adata[:, gene].X
        else:
            data = self.adata[:, gene].layers[self.norm_data_key]

        if scipy.sparse.issparse(data):
            data = data.toarray().ravel()

        return data.tolist()

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
