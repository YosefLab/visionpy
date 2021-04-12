from typing import Optional, Union

import anndata
import scipy

from ._compat import Literal


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
        if self._adata is None:
            raise ValueError("Accessor not populated with anndata.")
        if self.norm_data_key == "use_raw":
            data = self._adata.raw[:, gene].X
        elif self.norm_data_key is None:
            data = self._adata[:, gene].X
        else:
            data = self._adata[:, gene].layers[self.norm_data_key]

        if scipy.sparse.issparse(data):
            data = data.toarray().ravel()

        return data.tolist()
