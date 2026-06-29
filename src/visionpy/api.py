import os
from typing import Literal, List, Optional, Sequence, Union

import anndata
import numpy as np
import pandas as pd
import click
from sklearn.preprocessing import normalize
from statsmodels.stats.multitest import multipletests

from visionpy import create_app, data_accessor

from .knn import (
    compute_knn_weights_anndata,
    compute_knn_weights_from_tree_anndata,
    compute_knn_weights_from_tree_lca_anndata,
)
from .microclusters import apply_micro_clustering
from .projections import _ALL_PROJECTION_METHODS, compute_latent_space, generate_projections
from .signature import compute_signatures_anndata, split_signed_signatures


def _recompute_joint_fdr(adata: anndata.AnnData) -> None:
    """BH-correct p-values for signatures and metadata jointly.

    R VISION applies ``p.adjust(method="BH")`` over all tests (signatures +
    numerical meta + categorical meta) in a single call.  Reproducing that here
    prevents the FDR denominator from differing between the two tools.
    """
    obs_key = "vision_obs_df_scores"
    sig_key = "vision_signature_scores"

    obs_df  = adata.uns.get(obs_key)   # may be absent if only meta was requested
    sig_df  = adata.uns.get(sig_key)

    if obs_df is None or sig_df is None:
        return   # can't pool if one is missing

    pvals_obs = obs_df["pvals"].values
    pvals_sig = sig_df["pvals"].values

    all_pvals = np.concatenate([pvals_obs, pvals_sig])
    all_fdr   = multipletests(all_pvals, method="fdr_bh")[1]

    n_obs = len(pvals_obs)
    obs_df = obs_df.copy()
    sig_df = sig_df.copy()
    obs_df["fdr"] = all_fdr[:n_obs]
    sig_df["fdr"] = all_fdr[n_obs:]

    adata.uns[obs_key] = obs_df
    adata.uns[sig_key] = sig_df


def _pearsonr_matrix(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Pearson r between every column of A (n×p) and every column of B (n×q) → (p, q)."""
    A = A - A.mean(axis=0)
    B = B - B.mean(axis=0)
    norm_A = np.sqrt((A ** 2).sum(axis=0))
    norm_B = np.sqrt((B ** 2).sum(axis=0))
    norm_A[norm_A == 0] = 1.0
    norm_B[norm_B == 0] = 1.0
    return (A.T @ B) / np.outer(norm_A, norm_B)


def _linkage_to_newick(Z, labels: list) -> str:
    """Convert scipy linkage matrix to newick string."""
    n = len(labels)
    nodes = {i: labels[i] for i in range(n)}
    for i, (left, right, dist, _) in enumerate(Z):
        left, right = int(left), int(right)
        nodes[n + i] = f"({nodes[left]}:{dist / 2:.4f},{nodes[right]}:{dist / 2:.4f})"
    return nodes[2 * n - 2] + ";"


def _annotate_latent_components(adata):
    """Pearson correlations between PCA components and signature/numeric-obs scores."""
    if "X_pca" not in adata.obsm or "vision_signatures" not in adata.obsm:
        return
    from scipy.stats import t as t_dist

    latent = np.asarray(adata.obsm["X_pca"], dtype=float)   # (n_cells, n_pcs)
    n = latent.shape[0]
    pc_labels = [f"PC{i + 1}" for i in range(latent.shape[1])]

    sig_df = adata.obsm["vision_signatures"]
    corr = _pearsonr_matrix(latent, sig_df.to_numpy().astype(float)).T  # (n_sigs, n_pcs)
    t_stat = corr * np.sqrt(n - 2) / np.sqrt(np.maximum(1 - corr ** 2, 1e-10))
    pvals = 2 * t_dist.sf(np.abs(t_stat), df=n - 2)
    adata.uns["vision_lca"] = {
        "sig_labels": sig_df.columns.tolist(),
        "proj_labels": pc_labels,
        "zscores": corr.tolist(),
        "pvals": pvals.tolist(),
    }

    numeric_data = adata.obs._get_numeric_data()
    if numeric_data.shape[1] > 0:
        meta_corr = _pearsonr_matrix(latent, numeric_data.to_numpy().astype(float)).T  # (n_meta, n_pcs)
        t_meta = meta_corr * np.sqrt(n - 2) / np.sqrt(np.maximum(1 - meta_corr ** 2, 1e-10))
        adata.uns["vision_lca_meta"] = {
            "sig_labels": numeric_data.columns.tolist(),
            "proj_labels": pc_labels,
            "zscores": meta_corr.tolist(),
            "pvals": (2 * t_dist.sf(np.abs(t_meta), df=n - 2)).tolist(),
        }


def _compute_signature_dendrogram(adata):
    """Hierarchically cluster signatures by score vectors, store newick in uns."""
    if "vision_signatures" not in adata.obsm:
        return
    from scipy.cluster.hierarchy import linkage

    sig_df = adata.obsm["vision_signatures"]
    sig_names = sig_df.columns.tolist()
    if len(sig_names) < 3:
        return
    Z = linkage(sig_df.T.to_numpy(), method="ward", metric="euclidean")
    adata.uns["vision_dendrogram"] = _linkage_to_newick(Z, sig_names)


def _compute_signature_clusters(adata):
    """Cut signature dendrogram into groups for the filter-panel sidebar."""
    if "vision_signatures" not in adata.obsm:
        return
    from scipy.cluster.hierarchy import linkage, fcluster

    sig_df = adata.obsm["vision_signatures"]
    sig_names = sig_df.columns.tolist()
    n_sigs = len(sig_names)
    if n_sigs < 3:
        adata.uns["vision_sig_clusters"] = {s: 1 for s in sig_names}
        return
    Z = linkage(sig_df.T.to_numpy(), method="ward", metric="euclidean")
    n_groups = min(max(2, int(np.sqrt(n_sigs))), 10)
    labels = fcluster(Z, t=n_groups, criterion="maxclust")
    adata.uns["vision_sig_clusters"] = {s: int(labels[i]) for i, s in enumerate(sig_names)}


def _cluster_cells(adata, obsm_key=None, num_neighbors=None, exact_knn=False):
    """Louvain clustering matching R VISION's clusterCells().

    R VISION recomputes a fresh KNN with K=min(numNeighbors, 30) — separate from
    the K=sqrt(n) weight matrix used for Geary's C.  It builds a directed
    adjacency graph, converts to undirected with mode='each' (multi-edges for
    reciprocal pairs), and runs Louvain without explicit edge weights.

    When obsm_key is None (scanpy connectivities fallback path), falls back to
    the existing weight matrix topology.
    """
    import igraph as ig
    from .knn import find_knn

    if obsm_key is not None and obsm_key in adata.obsm:
        latent = np.asarray(adata.obsm[obsm_key], dtype=float)
        n = latent.shape[0]
        if num_neighbors is None:
            num_neighbors = max(1, round(np.sqrt(n)))
        K = min(num_neighbors, 30)          # R VISION caps clustering K at 30
        idx, _ = find_knn(latent, K, exact=exact_knn)
        edges = [(i, int(idx[i, j])) for i in range(n) for j in range(K)]
        g = ig.Graph(n=n, edges=edges, directed=True)
    else:
        # Connectivities fallback: use weight matrix topology
        W = adata.obsp["weights"].tocoo()
        edges = list(zip(W.row.tolist(), W.col.tolist()))
        g = ig.Graph(n=adata.n_obs, edges=edges, directed=True)

    # mode="each": each directed edge → one undirected edge; reciprocal pairs
    # become multi-edges, matching R's igraph::as.undirected(mode="each").
    g = g.as_undirected(mode="each")
    membership = g.community_multilevel().membership   # unweighted, matching R
    adata.obs["VISION_Clusters"] = pd.Categorical([str(m) for m in membership])


def _prepare_vision(
    adata: Union[str, anndata.AnnData],
    name: str,
    norm_data_key: Optional[Union[Literal["use_raw"], str]] = None,
    compute_neighbors_on_key: Optional[str] = None,
    tree: Optional[str] = None,
    lca_knn: bool = False,
    lca_min_size: int = 20,
    signature_varm_key: Optional[str] = None,
    signature_names_uns_key: Optional[str] = None,
    protein_obsm_key: Optional[str] = None,
    unnormalized_layer_key: Optional[str] = None,
    # Preprocessing (mirrors R VISION's analyze() pipeline)
    run_pca: bool = False,
    pca_max_components: int = 30,
    use_permutation_wpca: bool = False,
    run_projections: bool = False,
    projection_methods: Optional[Sequence[str]] = None,
    num_neighbors: Optional[int] = None,
    exact_knn: bool = False,
    run_micro_clustering: Union[bool, Literal["auto"]] = "auto",
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
    _neighbors_key = None
    if tree is not None:
        # PhyloVision mode: KNN weights from tree structure
        newick_str = tree
        if os.path.isfile(tree):
            with open(tree) as _f:
                newick_str = _f.read().strip()
        if lca_knn:
            compute_knn_weights_from_tree_lca_anndata(
                adata, newick_str, min_size=lca_min_size
            )
        else:
            compute_knn_weights_from_tree_anndata(adata, newick_str, K=num_neighbors)
        adata.uns["vision_tree"] = newick_str
    elif compute_neighbors_on_key is not None:
        compute_knn_weights_anndata(adata, obsm_key=compute_neighbors_on_key, K=num_neighbors, exact=exact_knn)
        _neighbors_key = compute_neighbors_on_key
    elif "X_pca" in adata.obsm:
        compute_knn_weights_anndata(adata, obsm_key="X_pca", K=num_neighbors, exact=exact_knn)
        _neighbors_key = "X_pca"
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
    # 2b. Louvain clustering on a fresh K=min(K,30) KNN → VISION_Clusters
    #     Skipped in PhyloVision mode unless a coordinate obsm key is available.
    # ------------------------------------------------------------------
    if _neighbors_key is not None:
        _cluster_cells(adata, obsm_key=_neighbors_key, num_neighbors=num_neighbors, exact_knn=exact_knn)

    # ------------------------------------------------------------------
    # 3. Optional micro-clustering ("auto" pools when n_obs > 100 000)
    # ------------------------------------------------------------------
    _should_pool = run_micro_clustering is True or (
        run_micro_clustering == "auto" and adata.n_obs > 100_000
    )
    if _should_pool:
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
    data_accessor.protein_obsm_key = protein_obsm_key
    # Store unnormalized key for blueprint display (raw counts for gene expression view)
    if unnormalized_layer_key is not None:
        adata.uns["vision_unnormalized_key"] = unnormalized_layer_key
    data_accessor.compute_obs_df_scores()
    data_accessor.compute_one_vs_all_obs_cols()
    data_accessor.persist_meta_differential()
    # Protein analysis (CITE-seq ADT) — no p-values, matching R VISION's fbConsistencyScores
    if protein_obsm_key is not None:
        data_accessor.compute_protein_autocorrelation()
        data_accessor.compute_protein_differential()

    # ------------------------------------------------------------------
    # 6. Optional signature scoring pipeline
    # ------------------------------------------------------------------
    if signature_varm_key is not None:
        # Auto-expand bidirectional signatures into _UP / _DOWN sub-columns
        split_signed_signatures(adata, varm_key=signature_varm_key)
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
        data_accessor.persist_signature_differential()
        data_accessor.compute_gene_score_per_signature()
        data_accessor.persist_gene_importance()
        _annotate_latent_components(adata)
        _compute_signature_dendrogram(adata)
        _compute_signature_clusters(adata)
        # Re-apply BH correction jointly over all tests (sigs + meta) to match R VISION
        _recompute_joint_fdr(adata)


def start_vision(
    adata: Union[str, anndata.AnnData],
    name: str,
    norm_data_key: Optional[Union[Literal["use_raw"], str]] = None,
    compute_neighbors_on_key: Optional[str] = None,
    tree: Optional[str] = None,
    lca_knn: bool = False,
    lca_min_size: int = 20,
    signature_varm_key: Optional[str] = None,
    signature_names_uns_key: Optional[str] = None,
    protein_obsm_key: Optional[str] = None,
    unnormalized_layer_key: Optional[str] = None,
    port: Optional[int] = None,
    debug: bool = False,
    run_pca: bool = False,
    pca_max_components: int = 30,
    use_permutation_wpca: bool = False,
    run_projections: bool = False,
    projection_methods: Optional[Sequence[str]] = None,
    num_neighbors: Optional[int] = None,
    run_micro_clustering: Union[bool, Literal["auto"]] = "auto",
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
        before building the KNN graph.
    pca_max_components
        Maximum number of PCs when *run_pca* is ``True``.  Default 30.
    use_permutation_wpca
        Use permutation WPCA (automatic PC selection) instead of fixed
        *pca_max_components*.  Only effective when *run_pca* is ``True``.
    run_projections
        If ``True``, generate 2D visualisation projections (tSNE, UMAP, …).
    projection_methods
        Subset of projection methods to run.  Defaults to all methods.
    num_neighbors
        Number of KNN neighbors for the VISION weight graph.  Defaults to
        ``round(sqrt(n_cells))``, matching R VISION.
    run_micro_clustering
        ``True`` always pools; ``False`` never pools; ``"auto"`` (default)
        pools automatically when ``n_obs > 100 000``, matching R VISION's
        ``pool="auto"`` behaviour.
    cells_per_partition
        Target cells per micro-cluster.  Default 10.
    device
        Device for GPU-accelerated signature scoring: ``"auto"``, ``"cuda"``,
        ``"mps"``, or ``"cpu"``.
    batch_size
        Number of signatures processed per GPU batch.  Default 1200.
    """
    if isinstance(adata, str):
        adata = anndata.read_h5ad(adata)

    if "vision_session_name" in adata.uns:
        # Already prepared — just wire the data accessor and launch.
        data_accessor.adata = adata
        data_accessor.norm_data_key = adata.uns.get("norm_data_key", norm_data_key)
        data_accessor.signature_varm_key = adata.uns.get("signature_varm_key", signature_varm_key)
        data_accessor.signature_names_uns_key = signature_names_uns_key
        data_accessor.protein_obsm_key = protein_obsm_key
    else:
        _prepare_vision(
            adata=adata,
            name=name,
            norm_data_key=norm_data_key,
            compute_neighbors_on_key=compute_neighbors_on_key,
            tree=tree,
            lca_knn=lca_knn,
            lca_min_size=lca_min_size,
            signature_varm_key=signature_varm_key,
            signature_names_uns_key=signature_names_uns_key,
            run_pca=run_pca,
            pca_max_components=pca_max_components,
            use_permutation_wpca=use_permutation_wpca,
            run_projections=run_projections,
            projection_methods=projection_methods,
            num_neighbors=num_neighbors,
            run_micro_clustering=run_micro_clustering,
            cells_per_partition=cells_per_partition,
            device=device,
            batch_size=batch_size,
            protein_obsm_key=protein_obsm_key,
            unnormalized_layer_key=unnormalized_layer_key,
        )
    app = create_app()
    app.run(host="0.0.0.0", threaded=True, processes=1, debug=debug, port=port)


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
@click.option("--num_neighbors", default=None, type=int, show_default=True, help="KNN neighbors for weight graph (default: sqrt(n_cells)).")
@click.option("--pool", "pool", default="auto", show_default=True,
              type=click.Choice(["true", "false", "auto"]),
              help="Micro-cluster cells: true | false | auto (pool when >100k cells).")
@click.option("--cells_per_partition", default=10, show_default=True, help="Target cells per micro-cluster.")
@click.option("--device", default="auto", show_default=True, help="Device for signature scoring: auto | cuda | mps | cpu.")
@click.option("--batch_size", default=1200, show_default=True, help="Signatures per GPU batch.")
@click.option("--port", default=None, type=int, help="Web server port (default 5000).")
@click.option("--debug", is_flag=True, default=False, help="Enable Flask debug mode.")
def _start_vision_cli(
    adata, name, norm_data_key, compute_neighbors_on_key,
    signature_varm_key, signature_names_uns_key,
    run_pca, pca_max_components, run_projections,
    num_neighbors, pool, cells_per_partition,
    device, batch_size, port, debug,
):
    pool_map = {"true": True, "false": False, "auto": "auto"}
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
        num_neighbors=num_neighbors,
        run_micro_clustering=pool_map[pool],
        cells_per_partition=cells_per_partition,
        device=device,
        batch_size=batch_size,
        port=port,
        debug=debug,
    )


# ---------------------------------------------------------------------------
# Post-construction helpers
# ---------------------------------------------------------------------------

def add_trajectory(
    adata: anndata.AnnData,
    milestone_network: pd.DataFrame,
    progressions: pd.DataFrame,
    K: Optional[int] = None,
) -> None:
    """Add a dynverse-format trajectory and compute trajectory autocorrelation.

    After calling this function the VISION UI will show a **Tree** tab with
    milestone-graph projections and per-signature trajectory autocorrelation.

    Parameters
    ----------
    adata
        AnnData already prepared by :func:`start_vision`.
    milestone_network
        DataFrame with columns ``from``, ``to``, ``length`` describing the
        milestone graph edges.
    progressions
        DataFrame indexed by cell_id (or with a ``cell_id`` column) with
        columns ``from``, ``to``, ``position`` (fraction ∈ [0, 1] along the
        from→to edge).
    K
        Number of trajectory nearest neighbours for the Geary's C weight
        graph.  Defaults to ``round(sqrt(n_cells))``.
    """
    from .trajectory import (
        _build_adj_matrix,
        calculate_trajectory_distances,
        compute_trajectory_knn_weights,
        create_trajectory_metadata,
        generate_trajectory_projections,
    )
    from .signature import _gearysc_for_dataframe

    # Normalise progressions index
    if "cell_id" in progressions.columns:
        progressions = progressions.set_index("cell_id")

    # Align progressions to adata cell order (inner join)
    common = progressions.index.intersection(adata.obs_names)
    if len(common) < len(progressions):
        import warnings
        warnings.warn(
            f"{len(progressions) - len(common)} cells in progressions not found in "
            "adata — they will be dropped.",
            stacklevel=2,
        )
    progressions = progressions.loc[common]

    adj, milestone_ids = _build_adj_matrix(milestone_network)

    # Geodesic distances → KNN weights
    dist_mat = calculate_trajectory_distances(progressions, adj, milestone_ids)
    traj_weights = compute_trajectory_knn_weights(dist_mat, K=K)
    adata.obsp["trajectory_weights"] = traj_weights

    # 2-D projections
    adata.uns["vision_trajectory_projections"] = generate_trajectory_projections(
        progressions, adj, milestone_ids
    )

    # Trajectory metadata → adata.obs
    traj_meta = create_trajectory_metadata(progressions, adj, milestone_ids)
    for col in traj_meta.columns:
        adata.obs[col] = np.nan
        adata.obs.loc[traj_meta.index, col] = traj_meta[col]

    # Autocorrelation (same Geary's C but with trajectory weights)
    traj_autocorr: dict = {}

    if "vision_signatures" in adata.obsm:
        sig_df = adata.obsm["vision_signatures"].loc[common]
        traj_autocorr["Signatures"] = _gearysc_for_dataframe(
            traj_weights, sig_df, compute_pvals=True
        )

    # All obs metadata (numeric + categorical handled inside _gearysc_for_dataframe)
    numeric_obs = adata.obs.loc[common]._get_numeric_data()
    if numeric_obs.shape[1] > 0:
        traj_autocorr["Meta"] = _gearysc_for_dataframe(
            traj_weights, numeric_obs, compute_pvals=True
        )

    # Protein (CITE-seq), no permutation p-values — matching R VISION
    prot_key = adata.uns.get("vision_protein_differential_key")
    if prot_key and prot_key in adata.obsm:
        mat = adata.obsm[prot_key]
        prot_names = adata.uns.get("vision_protein_autocorr", pd.DataFrame()).index.tolist()
        if isinstance(mat, pd.DataFrame):
            protein_df = mat.loc[common]
        else:
            protein_df = pd.DataFrame(
                np.asarray(mat)[adata.obs_names.get_indexer(common)],
                index=common,
                columns=prot_names or [f"Protein_{i}" for i in range(mat.shape[1])],
            )
        traj_autocorr["Proteins"] = _gearysc_for_dataframe(
            traj_weights, protein_df, compute_pvals=False
        )

    adata.uns["vision_trajectory_autocorr"] = traj_autocorr


def add_projection(
    adata: anndata.AnnData,
    name: str,
    coordinates: np.ndarray,
) -> None:
    """Add a precomputed 2-D projection to an already-prepared AnnData.

    Parameters
    ----------
    adata
        AnnData that has already been processed by :func:`start_vision` /
        :func:`_prepare_vision`.
    name
        Label shown in the VISION UI (stored as ``X_<name>`` in ``obsm``).
    coordinates
        Array of shape ``(n_cells, 2)`` with the projection coordinates.
    """
    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.shape != (adata.n_obs, 2):
        raise ValueError(
            f"coordinates must have shape ({adata.n_obs}, 2), got {coordinates.shape}"
        )
    key = f"X_{name}"
    adata.obsm[key] = coordinates


def add_signatures(
    adata: anndata.AnnData,
    signatures: Union[str, List[str]],
    norm_data_key: Optional[Union[Literal["use_raw"], str]] = None,
    min_signature_genes: int = 5,
    sig_gene_threshold: float = 0.001,
    device: str = "auto",
    batch_size: int = 1200,
) -> None:
    """Score additional gene signatures and append them to an existing VISION session.

    New signatures are merged with any already stored in
    ``adata.obsm["vision_signatures"]`` and the autocorrelation / differential /
    dendrogram / LCA are recomputed to include them.

    Parameters
    ----------
    adata
        AnnData already prepared by :func:`start_vision`.
    signatures
        Path(s) to ``*.gmt``, ``*.csv``, or ``*.txt`` signature files, or
        dictionaries mapping signature names to ``{UP: [...], DN: [...]}`` dicts.
    norm_data_key
        Layer key for normalised data (or ``'use_raw'``).  Defaults to ``None``
        (uses ``adata.X``).
    min_signature_genes, sig_gene_threshold
        Passed to :func:`~visionpy.signature.signatures_from_file`.
    device, batch_size
        Compute device and batch size for signature scoring.
    """
    from .anndata import AnnDataAccessor
    from .signature import signatures_from_file

    if isinstance(signatures, str):
        signatures = [signatures]

    # Preserve existing varm signatures (if any) so we can merge them later
    _TMP_VARM_KEY = "_vision_new_sigs"
    existing_varm = adata.varm.get("signatures")

    gmt_files = [s for s in signatures if isinstance(s, str)]
    dicts = [s for s in signatures if isinstance(s, dict)]

    signatures_from_file(
        adata,
        use_raw=(norm_data_key == "use_raw"),
        gmt_files=gmt_files or None,
        dicts=dicts or None,
        min_signature_genes=min_signature_genes,
        sig_gene_threshold=sig_gene_threshold,
    )
    new_varm = adata.varm["signatures"]  # (n_genes, n_new_sigs) DataFrame

    # Merge with pre-existing varm signatures; new sigs override on name collision
    if existing_varm is not None:
        if isinstance(existing_varm, pd.DataFrame):
            old_df = existing_varm
        else:
            old_df = pd.DataFrame(existing_varm, index=adata.var_names)
        merged_varm = old_df.copy()
        for col in new_varm.columns:
            merged_varm[col] = new_varm[col]
        adata.varm["signatures"] = merged_varm
    # else: new_varm is already in adata.varm["signatures"]

    # Score all signatures (new + existing) in one pass
    compute_signatures_anndata(
        adata,
        norm_data_key=norm_data_key,
        signature_varm_key="signatures",
        device=device,
        batch_size=batch_size,
    )

    # Recompute autocorrelation, differential, dendrogram, LCA
    accessor = AnnDataAccessor(adata)
    accessor.sig_adata = None  # reset cached adata
    accessor.compute_signature_scores()
    accessor.compute_one_vs_all_signatures()
    accessor.persist_signature_differential()
    accessor.compute_gene_score_per_signature()
    accessor.persist_gene_importance()

    _recompute_joint_fdr(adata)
    _annotate_latent_components(adata)
    _compute_signature_dendrogram(adata)
    _compute_signature_clusters(adata)


def add_latent_space(
    adata: anndata.AnnData,
    name: str,
    latent_space: np.ndarray,
    generate_projections_methods: Optional[Sequence[str]] = None,
) -> None:
    """Store a new latent space and optionally regenerate 2-D projections from it.

    Parameters
    ----------
    adata
        AnnData already prepared by :func:`start_vision`.
    name
        Label for the latent space (stored as ``X_<name>`` in ``obsm``).
    latent_space
        Array of shape ``(n_cells, n_dims)``.
    generate_projections_methods
        If provided, generate 2-D projections using these methods
        (subset of :data:`~visionpy.projections._ALL_PROJECTION_METHODS`).
        Pass an empty list or ``None`` to skip projection generation.
    """
    latent_space = np.asarray(latent_space, dtype=float)
    if latent_space.shape[0] != adata.n_obs:
        raise ValueError(
            f"latent_space must have {adata.n_obs} rows, got {latent_space.shape[0]}"
        )
    key = f"X_{name}"
    adata.obsm[key] = latent_space

    if generate_projections_methods:
        new_coords = generate_projections(latent_space, methods=generate_projections_methods)
        for method, coords in new_coords.items():
            proj_key = f"X_{method}"
            adata.obsm[proj_key] = coords
