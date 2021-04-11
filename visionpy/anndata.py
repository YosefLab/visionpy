from typing import Optional

import anndata


class AnnDataAccessor(object):
    def __init__(self, adata: Optional[anndata.AnnData] = None) -> None:
        self._adata = adata

    @property
    def adata(self):
        return self._adata

    @adata.setter
    def adata(self, adata: anndata.AnnData):
        self._adata = adata
