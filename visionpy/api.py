from typing import Optional, Sequence, Union

import anndata
import scanpy as sc
import numpy as np
import pandas as pd
import click
from sklearn.preprocessing import normalize
from visionpy import create_app, data_accessor
from scanpy.preprocessing._utils import _get_mean_var
from ._compat import Literal


def start_vision(
    adata: Union[str, anndata.AnnData],
    name: str,
    norm_data_key: Optional[Union[Literal["use_raw"], str]] = None,
    compute_neighbors_on_key: Optional[str] = None,
    signature_varm_key: Optional[str] = None,
    signature_names_uns_key: Optional[str] = None,
    signatures_files: Sequence[str] = None,
    use_raw_for_signatures: bool = False,
    debug: bool = False,
):
    """Wrapper function to start VISION server.

    Parameters
    ----------
    adata
        AnnData object
    name
        Name for the VISION session
    norm_data_key
        Key for layer with log library size normalized data. If
        `None` (default), uses `adata.X`
    compute_neighbors_on_key
        Key in `adata.obsm` to use for computing neighbors. If `None`, use
        neighbors stored in `adata`. If no neighbors have been previously
        computed an error will be raised.
    """

    if isinstance(adata, str):
        adata = anndata.read(str)

    if compute_neighbors_on_key is not None:
        sc.pp.neighbors(adata, use_rep=compute_neighbors_on_key)
    else:
        # TODO: check that neighbors is in anndata
        pass

    # row normalize connectivity
    adata.obsp["normalized_connectivities"] = normalize(
        adata.obsp["connectivities"], norm="l1", axis=1
    )
    adata.uns["vision_session_name"] = name

    # compute signatures
    if signature_varm_key is not None:
        compute_signature(
            adata,
            norm_data_key,
            use_raw_for_signatures,
            signature_varm_key,
            signature_names_uns_key,
        )

    data_accessor.adata = adata
    data_accessor.norm_data_key = norm_data_key

    data_accessor.compute_obs_df_scores()
    data_accessor.compute_signature_scores()
    data_accessor.compute_one_vs_all_signatures()

    app = create_app()
    app.run(threaded=False, processes=1, debug=debug)


def compute_signature(
    adata: anndata.AnnData,
    norm_data_key: str,
    use_raw_for_signatures: bool,
    signature_varm_key: str,
    signature_names_uns_key: str,
):
    if norm_data_key is None:
        gene_expr = adata.X
    elif norm_data_key == "use_raw":
        gene_expr = adata.raw.X
    else:
        gene_expr = adata.layers[norm_data_key]
    sig_matrix = (
        adata.varm[signature_varm_key]
        if not use_raw_for_signatures
        else adata.raw.varm[signature_varm_key]
    )
    cell_signature_matrix = np.asarray(gene_expr @ sig_matrix)
    if signature_names_uns_key is not None:
        cols = adata.uns[signature_names_uns_key]
    else:
        cols = None
    sig_df = pd.DataFrame(
        data=cell_signature_matrix, columns=cols, index=adata.obs_names
    )
    # TODO: might not be most efficient
    total = np.asarray((sig_matrix != 0).sum(0)).ravel()
    sig_df = sig_df / total
    # normalize
    mean, var = _get_mean_var(gene_expr, axis=1)
    mean = mean.reshape(-1, 1) * np.asarray(sig_matrix.sum(0)).ravel()
    mean = mean / total
    var = var.reshape(-1, 1) / total
    sig_df = (sig_df - mean) / np.sqrt(var)

    adata.obsm["vision_signatures"] = sig_df


@click.command()
@click.option("--adata", prompt="Path to anndata.", help="Test")
@click.option("--name", prompt="Name of VISION session.")
@click.option("--norm_data_key", default=None, prompt="Name of VISION session.")
@click.option(
    "--compute_neighbors_on_key", default=None, prompt="Name of VISION session."
)
@click.option("--debug", default=False, prompt="Name of VISION session.")
def _start_vision_cli(adata, name, norm_data_key, compute_neighbors_on_key, debug):
    start_vision(adata, name, norm_data_key, compute_neighbors_on_key, debug)
