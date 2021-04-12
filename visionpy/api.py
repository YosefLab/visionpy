from typing import Optional, Union

import anndata
import scanpy as sc
from visionpy import create_app, data_accessor

from ._compat import Literal


def start_vision(
    adata: Union[str, anndata.AnnData],
    norm_data_key: Optional[Union[Literal["use_raw"], str]] = None,
    compute_neighbors_on_key: Optional[str] = None,
    debug: bool = False,
):
    """Wrapper function to start VISION server.

    Parameters
    ----------
    adata
        AnnData object
    compute_neighbors_on_key
        Key in `adata.obsm` to use for computing neighbors. If `None`, use
        neighbors stored in `adata`. If no neighbors have been previously
        computed an error will be raised.
    norm_data_key
        Key for layer with log library size normalized data. If
        `None` (default), uses `adata.X`
    """

    if isinstance(adata, str):
        adata = anndata.read(str)

    if compute_neighbors_on_key is not None:
        sc.pp.neighbors(adata, use_rep=compute_neighbors_on_key)
    else:
        # TODO: check that neighbors is in anndata
        pass

    data_accessor.adata = adata
    data_accessor.norm_data_key = norm_data_key

    app = create_app()
    app.run(threaded=False, processes=1, debug=debug)
