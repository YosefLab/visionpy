from typing import Optional

import anndata


class AnnDataAccessor(object):
    def __init__(
        self,
        adata: Optional[anndata.AnnData] = None,
        protein_obsm_key: Optional[str] = None,
    ) -> None:
        self._adata = adata
        self._protein_obsm_key = protein_obsm_key

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
