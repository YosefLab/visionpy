from typing import Literal, List, Optional, Sequence, Union

import anndata
import click
from sklearn.preprocessing import normalize

from visionpy import create_app, data_accessor

from .knn import compute_knn_weights_anndata
from .microclusters import apply_micro_clustering
from .projections import _ALL_PROJECTION_METHODS, compute_latent_space, generate_projections
from .signature import compute_signatures_anndata


def _prepare_vision(
    adata: Union[str, anndata.AnnData],
    name: str,
    norm_data_key: Optional[Union[Literal["use_raw"], str]] = None,
    compute_neighbors_on_key: Optional[str] = None,
    signature_varm_key: Optional[str] = None,
    signature_names_uns_key: Optional[str] = None,
    # Preprocessing (mirrors R VISION's analyze() pipeline)
    run_pca: bool = False,
    pca_max_components: int = 30,
    use_permutation_wpca: bool = False,
    run_projections: bool = False,
    projection_methods: Optional[Sequence[str]] = None,
    run_micro_clustering: bool = False,
    cells_per_partition: int = 10,
    # GPU / batch parameters for signature scoring
    device: str = "auto",
    batch_size: int = 1200,
):
    if isinstance(adata, str):
        adata = anndata.read_h5ad(adata)

    # ------------------------------------------------------------------
    # 1. Optional PCA (must run before KNN weights if latent space absent)
    # ------------------------------------------------------------------
    if run_pca:
        compute_latent_space(
            adata,
            norm_data_key=norm_data_key if norm_data_key != "use_raw" else None,
            max_components=pca_max_components,
            use_permutation_wpca=use_permutation_wpca,
        )

    # ------------------------------------------------------------------
    # 2. VISION exponential-kernel KNN weights → obsp["weights"]
    #    Geary's C in signature.py reads this key.
    # ------------------------------------------------------------------
    if compute_neighbors_on_key is not None:
        compute_knn_weights_anndata(adata, obsm_key=compute_neighbors_on_key)
    elif "X_pca" in adata.obsm:
        compute_knn_weights_anndata(adata, obsm_key="X_pca")
    elif "connectivities" in adata.obsp:
        # Fallback: L1-normalise scanpy connectivities into the weights slot
        adata.obsp["weights"] = normalize(adata.obsp["connectivities"], norm="l1", axis=1)
    else:
        raise ValueError(
            "No latent space or neighbor graph found in adata. "
            "Pass run_pca=True, set compute_neighbors_on_key to an obsm key, "
            "or pre-compute neighbors with sc.pp.neighbors()."
        )

    # ------------------------------------------------------------------
    # 3. Optional micro-clustering
    # ------------------------------------------------------------------
    if run_micro_clustering:
        apply_micro_clustering(adata, cells_per_partition=cells_per_partition)

    # ------------------------------------------------------------------
    # 4. Optional 2D projections
    # ------------------------------------------------------------------
    if run_projections:
        methods = projection_methods if projection_methods is not None else _ALL_PROJECTION_METHODS
        generate_projections(adata, projection_methods=methods)

    # ------------------------------------------------------------------
    # 5. Wire up data accessor and compute analysis scores
    # ------------------------------------------------------------------
    adata.uns["vision_session_name"] = name

    data_accessor.adata = adata
    data_accessor.norm_data_key = norm_data_key
    data_accessor.signature_varm_key = signature_varm_key
    data_accessor.signature_names_uns_key = signature_names_uns_key
    data_accessor.compute_obs_df_scores()
    data_accessor.compute_one_vs_all_obs_cols()

    # ------------------------------------------------------------------
    # 6. Optional signature scoring pipeline
    # ------------------------------------------------------------------
    if signature_varm_key is not None:
        compute_signatures_anndata(
            adata,
            norm_data_key,
            signature_varm_key,
            signature_names_uns_key,
            device=device,
            batch_size=batch_size,
        )
        data_accessor.compute_signature_scores()
        data_accessor.compute_one_vs_all_signatures()
        data_accessor.compute_gene_score_per_signature()


def start_vision(
    adata: Union[str, anndata.AnnData],
    name: str,
    norm_data_key: Optional[Union[Literal["use_raw"], str]] = None,
    compute_neighbors_on_key: Optional[str] = None,
    signature_varm_key: Optional[str] = None,
    signature_names_uns_key: Optional[str] = None,
    port: Optional[int] = None,
    debug: bool = False,
    run_pca: bool = False,
    pca_max_components: int = 30,
    use_permutation_wpca: bool = False,
    run_projections: bool = False,
    projection_methods: Optional[Sequence[str]] = None,
    run_micro_clustering: bool = False,
    cells_per_partition: int = 10,
    device: str = "auto",
    batch_size: int = 1200,
):
    """Start the VISION web server after running the analysis pipeline.

    Parameters
    ----------
    adata
        AnnData object or path to a saved ``.h5ad`` file.
    name
        Display name for the VISION session.
    norm_data_key
        Key in ``adata.layers`` containing log-library-size-normalised counts.
        Use ``"use_raw"`` to pull from ``adata.raw``. If ``None`` (default),
        ``adata.X`` is used directly.
    compute_neighbors_on_key
        Key in ``adata.obsm`` to use as the latent space when building the
        VISION KNN weight graph.  If ``None`` and ``adata.obsm["X_pca"]``
        exists, that is used automatically.  Falls back to L1-normalising
        ``adata.obsp["connectivities"]`` if available.
    signature_varm_key
        Key in ``adata.varm`` for the genes × signatures matrix (1 = positive
        gene, -1 = negative gene, 0 = absent).  If ``None`` (default), no
        signature analysis is run.
    signature_names_uns_key
        Key in ``adata.uns`` holding signature names.  If ``None``, column
        names from the ``varm`` DataFrame are used, or auto-generated labels.
    port
        Port for the Flask web server.  Defaults to 5000.
    debug
        Enable Flask debug mode.
    run_pca
        If ``True``, run :func:`~visionpy.projections.compute_latent_space`
        before building the KNN graph.  Ignored when the latent space is
        already present in ``adata.obsm``.
    pca_max_components
        Maximum number of PCs when *run_pca* is ``True``.  Default 30.
    use_permutation_wpca
        Use permutation WPCA (automatic PC selection) instead of fixed
        *pca_max_components*.  Only effective when *run_pca* is ``True``.
    run_projections
        If ``True``, generate 2D visualisation projections (tSNE, UMAP, …)
        via :func:`~visionpy.projections.generate_projections`.
    projection_methods
        Subset of projection methods to run.  Defaults to all five methods
        when *run_projections* is ``True``.
    run_micro_clustering
        If ``True``, partition cells into micro-clusters via
        :func:`~visionpy.microclusters.apply_micro_clustering`.
    cells_per_partition
        Target cells per micro-cluster when *run_micro_clustering* is ``True``.
        Default 10.
    device
        Device for GPU-accelerated signature scoring: ``"auto"`` (default,
        tries CUDA then MPS then CPU), ``"cuda"``, ``"mps"``, or ``"cpu"``.
    batch_size
        Number of signatures processed per GPU batch.  Default 1200.
    """
    _prepare_vision(
        adata=adata,
        name=name,
        norm_data_key=norm_data_key,
        compute_neighbors_on_key=compute_neighbors_on_key,
        signature_varm_key=signature_varm_key,
        signature_names_uns_key=signature_names_uns_key,
        run_pca=run_pca,
        pca_max_components=pca_max_components,
        use_permutation_wpca=use_permutation_wpca,
        run_projections=run_projections,
        projection_methods=projection_methods,
        run_micro_clustering=run_micro_clustering,
        cells_per_partition=cells_per_partition,
        device=device,
        batch_size=batch_size,
    )
    app = create_app()
    app.run(threaded=False, processes=1, debug=debug, port=port)


@click.command()
@click.option("--adata", required=True, help="Path to .h5ad file.")
@click.option("--name", required=True, help="Display name for the VISION session.")
@click.option("--norm_data_key", default=None, show_default=True, help="Layer key for normalised counts (or 'use_raw').")
@click.option("--compute_neighbors_on_key", default=None, show_default=True, help="obsm key for latent space used to build KNN weights.")
@click.option("--signature_varm_key", default=None, show_default=True, help="varm key for the genes × signatures matrix.")
@click.option("--signature_names_uns_key", default=None, show_default=True, help="uns key for signature names.")
@click.option("--run_pca", is_flag=True, default=False, help="Compute PCA latent space before analysis.")
@click.option("--pca_max_components", default=30, show_default=True, help="Max PCs when --run_pca is set.")
@click.option("--run_projections", is_flag=True, default=False, help="Generate 2D projections (tSNE, UMAP, …).")
@click.option("--run_micro_clustering", is_flag=True, default=False, help="Partition cells into micro-clusters.")
@click.option("--cells_per_partition", default=10, show_default=True, help="Target cells per micro-cluster.")
@click.option("--device", default="auto", show_default=True, help="Device for signature scoring: auto | cuda | mps | cpu.")
@click.option("--batch_size", default=1200, show_default=True, help="Signatures per GPU batch.")
@click.option("--port", default=None, type=int, help="Web server port (default 5000).")
@click.option("--debug", is_flag=True, default=False, help="Enable Flask debug mode.")
def _start_vision_cli(
    adata, name, norm_data_key, compute_neighbors_on_key,
    signature_varm_key, signature_names_uns_key,
    run_pca, pca_max_components, run_projections,
    run_micro_clustering, cells_per_partition,
    device, batch_size, port, debug,
):
    start_vision(
        adata=adata,
        name=name,
        norm_data_key=norm_data_key,
        compute_neighbors_on_key=compute_neighbors_on_key,
        signature_varm_key=signature_varm_key,
        signature_names_uns_key=signature_names_uns_key,
        run_pca=run_pca,
        pca_max_components=pca_max_components,
        run_projections=run_projections,
        run_micro_clustering=run_micro_clustering,
        cells_per_partition=cells_per_partition,
        device=device,
        batch_size=batch_size,
        port=port,
        debug=debug,
    )
