"""KNN graph construction with exponential distance weighting.

Mirrors R VISION's ``computeKNNWeights`` (matrix method) in ``Projections.R``
and ``find_knn_parallel`` in ``Utilities.R``.

The weight kernel is W_ij = exp(−d²_ij / σ²_i), where σ_i is the distance
from cell i to its K-th (farthest) neighbor.  Rows are then L1-normalised so
each cell's weights sum to 1.  The result is stored as a ``(n_cells, n_cells)``
CSR sparse matrix in ``adata.obsp``, replacing the L1-normalised scanpy
connectivities that visionpy previously used as a proxy.
"""

from __future__ import annotations

import logging
from typing import Optional, Tuple

import numpy as np
from anndata import AnnData
from scipy.sparse import csr_matrix
from sklearn.neighbors import NearestNeighbors

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Low-level KNN search
# ---------------------------------------------------------------------------

def find_knn(
    data: np.ndarray,
    K: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Find K nearest neighbors for every point in *data*.

    Self-distances are excluded from the result (mirrors R VISION's
    ``find_knn_parallel``, which queries K+1 neighbors and drops the first
    column).

    Parameters
    ----------
    data : ndarray of shape (n_cells, n_dims)
        Coordinates in the latent space.
    K : int
        Number of neighbors to return per cell.

    Returns
    -------
    indices : ndarray of shape (n_cells, K)
        0-based column indices of the K nearest neighbors for each cell.
    distances : ndarray of shape (n_cells, K)
        Euclidean distances to each of the K nearest neighbors.
    """
    # Query K+1 so we can drop the self-hit at position 0
    nbrs = NearestNeighbors(n_neighbors=K + 1, algorithm="auto").fit(data)
    distances, indices = nbrs.kneighbors(data)

    # Column 0 is the point itself (distance 0); drop it
    return indices[:, 1:], distances[:, 1:]


# ---------------------------------------------------------------------------
# Exponential weight kernel + row normalisation
# ---------------------------------------------------------------------------

def compute_knn_weights(
    latent_space: np.ndarray,
    K: Optional[int] = None,
) -> csr_matrix:
    """Build an exponential-kernel KNN weight matrix from a latent space.

    Mirrors R VISION's ``computeKNNWeights`` (matrix S4 method):

    1. Find the K nearest neighbors in Euclidean space.
    2. Per-cell bandwidth: σ_i = max distance among the K neighbors.
    3. Weight kernel: W_ij = exp(−d²_ij / σ²_i).
    4. Row-normalise so each row sums to 1.

    Parameters
    ----------
    latent_space : ndarray of shape (n_cells, n_dims)
        Cell coordinates (e.g. PCA scores stored in ``adata.obsm["X_pca"]``).
    K : int, optional
        Number of neighbors.  Defaults to ``round(sqrt(n_cells))``,
        matching R VISION's default.

    Returns
    -------
    scipy.sparse.csr_matrix of shape (n_cells, n_cells)
        Row-normalised sparse weight matrix.  Entry (i, j) is the weight
        cell i assigns to cell j; only KNN connections are non-zero.

    Notes
    -----
    Unlike R VISION, which returns raw ``(indices, weights)`` arrays and
    assembles the full sparse matrix lazily inside callers such as
    ``sigsVsProjection_pcf``, this function builds the complete CSR matrix
    immediately.  This matches the format already expected by visionpy's
    existing ``adata.obsp["weights"]`` slot.
    """
    latent_space = np.asarray(latent_space, dtype=float)
    n_cells = latent_space.shape[0]

    if K is None:
        K = max(1, round(np.sqrt(n_cells)))

    logger.info("Building KNN graph: n_cells=%d, K=%d.", n_cells, K)

    idx, d = find_knn(latent_space, K)  # both (n_cells, K)

    # Per-cell bandwidth σ_i = distance to the K-th (farthest) neighbor
    sigma = d.max(axis=1, keepdims=True)   # (n_cells, 1)
    sigma[sigma == 0] = 1.0                # guard: all K neighbors co-located

    # Exponential kernel: W_ij = exp(−d²_ij / σ²_i)
    W = np.exp(-(d ** 2) / (sigma ** 2))  # (n_cells, K)

    # Row-normalise
    row_sums = W.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    W = W / row_sums

    # Assemble CSR sparse matrix
    rows = np.repeat(np.arange(n_cells), K)  # (n_cells * K,)
    cols = idx.ravel()                        # (n_cells * K,)
    vals = W.ravel()

    weights = csr_matrix((vals, (rows, cols)), shape=(n_cells, n_cells))
    logger.info("KNN weight matrix built — nnz=%d.", weights.nnz)
    return weights


# ---------------------------------------------------------------------------
# AnnData wrapper
# ---------------------------------------------------------------------------

def compute_knn_weights_anndata(
    adata: AnnData,
    obsm_key: str = "X_pca",
    K: Optional[int] = None,
    obsp_key: str = "weights",
) -> AnnData:
    """Compute VISION exponential KNN weights and store in ``adata.obsp``.

    Reads the latent space from ``adata.obsm[obsm_key]`` and writes the
    resulting ``(n_cells, n_cells)`` sparse weight matrix to
    ``adata.obsp[obsp_key]``.  This replaces the L1-normalised scanpy
    connectivities previously used by visionpy as a proxy for VISION weights.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.  Must have ``obsm_key`` in ``adata.obsm``
        (e.g. produced by :func:`~visionpy.projections.compute_latent_space`).
    obsm_key : str, optional
        Key in ``adata.obsm`` for the latent-space coordinates,
        by default ``"X_pca"``.
    K : int, optional
        Number of neighbors.  Defaults to ``round(sqrt(n_cells))``.
    obsp_key : str, optional
        Key under which to store the weight matrix in ``adata.obsp``,
        by default ``"weights"``.

    Returns
    -------
    AnnData
        The input *adata* modified in-place with
        ``adata.obsp[obsp_key]`` populated.

    Raises
    ------
    KeyError
        If *obsm_key* is not present in ``adata.obsm``.
    """
    if obsm_key not in adata.obsm:
        raise KeyError(
            f"'{obsm_key}' not found in adata.obsm. "
            "Run compute_latent_space() first or pass a valid obsm_key."
        )

    latent_space = np.asarray(adata.obsm[obsm_key], dtype=float)
    weights = compute_knn_weights(latent_space, K=K)
    adata.obsp[obsp_key] = weights

    logger.info(
        "KNN weights stored in adata.obsp['%s'] — shape %s.",
        obsp_key,
        weights.shape,
    )
    return adata
