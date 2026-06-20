"""Smoke tests for functions added in the scalability + integration pass.

All tests use small synthetic AnnData objects (no network downloads) so they
run in CI without hitting external servers.
"""

import numpy as np
import pandas as pd
import pytest
import anndata
from scipy.sparse import csr_matrix


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_adata(n_cells: int = 60, n_genes: int = 80, seed: int = 0) -> anndata.AnnData:
    rng = np.random.default_rng(seed)
    X = csr_matrix(rng.poisson(1.0, (n_cells, n_genes)).astype(np.float32))
    adata = anndata.AnnData(X)
    adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
    adata.obs_names = [f"Cell_{i}" for i in range(n_cells)]
    adata.obs["cluster"] = pd.Categorical(["A"] * (n_cells // 2) + ["B"] * (n_cells - n_cells // 2))
    return adata


# ---------------------------------------------------------------------------
# microclusters
# ---------------------------------------------------------------------------

def test_apply_micro_clustering_returns_pools():
    from visionpy.microclusters import apply_micro_clustering

    adata = _make_adata()
    pools = apply_micro_clustering(adata, cells_per_partition=10)

    assert isinstance(pools, dict)
    assert len(pools) >= 1
    # Every cell must appear in exactly one pool
    all_cells = [c for cells in pools.values() for c in cells]
    assert len(all_cells) == adata.n_obs
    assert set(all_cells) == set(adata.obs_names)
    # Result stored in uns
    assert "micro_clusters" in adata.uns


def test_pool_matrix_shape():
    from visionpy.microclusters import apply_micro_clustering, pool_matrix

    adata = _make_adata()
    pools = apply_micro_clustering(adata, cells_per_partition=10)
    pooled = pool_matrix(adata.X, pools, list(adata.obs_names))

    assert pooled.shape == (len(pools), adata.n_vars)
    # Means should be non-negative (Poisson counts)
    assert (pooled >= 0).all()


# ---------------------------------------------------------------------------
# projections
# ---------------------------------------------------------------------------

def test_compute_latent_space():
    from visionpy.projections import compute_latent_space

    adata = _make_adata()
    compute_latent_space(adata, max_components=5)

    assert "X_pca" in adata.obsm
    assert adata.obsm["X_pca"].shape[0] == adata.n_obs


def test_generate_projections_tsne():
    from visionpy.projections import compute_latent_space, generate_projections

    adata = _make_adata()
    compute_latent_space(adata, max_components=5)
    results = generate_projections(adata, projection_methods=["tSNE10"])

    assert "tSNE10" in results
    assert results["tSNE10"].shape == (adata.n_obs, 2)
    assert "X_tSNE10" in adata.obsm


# ---------------------------------------------------------------------------
# signature scoring
# ---------------------------------------------------------------------------

def _make_sig_adata(n_cells: int = 60, n_genes: int = 80, n_sigs: int = 5) -> anndata.AnnData:
    adata = _make_adata(n_cells=n_cells, n_genes=n_genes)
    sig_matrix = np.zeros((n_genes, n_sigs), dtype=np.float32)
    sig_matrix[:10, :] = 1.0   # first 10 genes positive
    sig_matrix[-5:, :] = -1.0  # last 5 genes negative
    adata.varm["sigs"] = sig_matrix
    adata.uns["sig_names"] = [f"Sig_{i}" for i in range(n_sigs)]
    return adata


def test_compute_signatures_anndata():
    from visionpy.signature import compute_signatures_anndata

    adata = _make_sig_adata()
    scores = compute_signatures_anndata(
        adata,
        norm_data_key=None,
        signature_varm_key="sigs",
        signature_names_uns_key="sig_names",
        device="cpu",
    )

    assert scores.shape == (adata.n_obs, 5)
    assert "vision_signatures" in adata.obsm


def test_batch_sig_eval_norm_cpu():
    from visionpy.normalization import get_normalized_copy_sparse
    from visionpy.signature import batch_sig_eval_norm

    adata = _make_sig_adata()
    X = adata.X  # sparse
    sig_matrix = adata.varm["sigs"]

    norm_data = get_normalized_copy_sparse(X, method="znorm_columns")
    scores = batch_sig_eval_norm(norm_data, sig_matrix, device="cpu", batch_size=3)

    assert scores.shape == (adata.n_obs, 5)
    assert np.isfinite(scores).all()


# ---------------------------------------------------------------------------
# pipeline integration
# ---------------------------------------------------------------------------

def test_prepare_vision_runs():
    from visionpy.api import _prepare_vision

    adata = _make_sig_adata()
    # Need a latent space for KNN weights
    from visionpy.projections import compute_latent_space
    compute_latent_space(adata, max_components=5)

    _prepare_vision(
        adata,
        name="test",
        norm_data_key=None,
        signature_varm_key="sigs",
        signature_names_uns_key="sig_names",
    )

    assert "vision_signatures" in adata.obsm
    assert "vision_signature_scores" in adata.uns
    assert "weights" in adata.obsp


def test_prepare_vision_with_preprocessing():
    from visionpy.api import _prepare_vision

    adata = _make_sig_adata()

    _prepare_vision(
        adata,
        name="test",
        norm_data_key=None,
        signature_varm_key="sigs",
        signature_names_uns_key="sig_names",
        run_pca=True,
        pca_max_components=5,
        run_micro_clustering=True,
        cells_per_partition=10,
    )

    assert "X_pca" in adata.obsm
    assert "micro_clusters" in adata.uns
    assert "weights" in adata.obsp
    assert "vision_signatures" in adata.obsm
