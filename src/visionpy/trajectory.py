"""Trajectory analysis for visionpy — port of R VISION's trajectory module.

Accepts dynverse/dynwrap-format trajectories (``milestone_network`` +
``progressions``) and provides:
  - Geodesic distance computation between cells
  - Trajectory KNN weights (Gaussian kernel, same formula as the latent-space KNN)
  - 2-D layout projections (FR / Tree / MDS / DH)
  - Per-edge trajectory metadata for adata.obs
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import igraph as ig
import numpy as np
import pandas as pd
import scipy.sparse
from scipy.sparse.csgraph import shortest_path


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _build_adj_matrix(
    milestone_network: pd.DataFrame,
) -> Tuple[np.ndarray, List[str]]:
    """Build a symmetric weighted adjacency matrix from *milestone_network*.

    Returns
    -------
    adj : (n_milestones, n_milestones) float64 array
    milestone_ids : sorted list of milestone name strings
    """
    milestone_ids = sorted(
        set(milestone_network["from"]) | set(milestone_network["to"])
    )
    n = len(milestone_ids)
    idx = {m: i for i, m in enumerate(milestone_ids)}
    adj = np.zeros((n, n), dtype=float)
    for _, row in milestone_network.iterrows():
        i, j = idx[row["from"]], idx[row["to"]]
        adj[i, j] = adj[j, i] = float(row["length"])
    return adj, milestone_ids


def _get_edges(adj: np.ndarray) -> List[Tuple[int, int]]:
    """Upper-triangle non-zero entries of *adj* as (i, j) index pairs."""
    rows, cols = np.where(np.triu(adj) > 0)
    return list(zip(rows.tolist(), cols.tolist()))


def _calc_intra_edge_dist_mat(edge_len: float, edge_pos: np.ndarray) -> np.ndarray:
    """(n+2) × (n+2) absolute-difference matrix including edge endpoints.

    Row/column 0 = 'from' endpoint (position 0).
    Row/column -1 = 'to' endpoint (position 1).
    Rows/columns 1..-2 = cells.
    """
    pos = np.concatenate([[0.0], edge_pos, [1.0]]) * edge_len
    return np.abs(pos[:, None] - pos[None, :])


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def calculate_trajectory_distances(
    progressions: pd.DataFrame,
    adj: np.ndarray,
    milestone_ids: List[str],
    chunk_size: int = 2000,
) -> np.ndarray:
    """Pairwise geodesic distances along the trajectory graph.

    Ported from R VISION's ``calculateTrajectoryDistances``.  Uses a chunked
    computation so the full n × n matrix is never materialised in memory.

    Parameters
    ----------
    progressions
        DataFrame indexed by cell_id with columns ``from``, ``to``,
        ``position`` (fraction ∈ [0, 1] along the from→to edge).
    adj
        (n_milestones, n_milestones) symmetric weighted adjacency matrix.
    milestone_ids
        Ordered milestone name list (matching ``adj`` row/col order).
    chunk_size
        Number of rows computed per chunk.  Tune to fit in RAM.

    Returns
    -------
    (n_cells, n_cells) float32 distance matrix.
    """
    id_to_idx = {m: i for i, m in enumerate(milestone_ids)}
    n = len(progressions)

    # All-pairs shortest paths in the milestone graph
    node_dist = shortest_path(adj, directed=False)  # (M, M)

    from_idx = progressions["from"].map(id_to_idx).values.astype(int)  # (n,)
    to_idx   = progressions["to"].map(id_to_idx).values.astype(int)    # (n,)
    pos      = progressions["position"].values.astype(float)            # (n,)
    edge_len = adj[from_idx, to_idx]                                    # (n,)

    dist_to_from = pos * edge_len          # distance to 'from' endpoint
    dist_to_to   = (1.0 - pos) * edge_len  # distance to 'to' endpoint

    # Pre-fetch rows of node_dist for each cell's endpoints
    nd_from = node_dist[from_idx]  # (n, M)
    nd_to   = node_dist[to_idx]    # (n, M)

    dist_mat = np.empty((n, n), dtype=np.float32)

    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)

        # Distance from each cell in [start:end] to every other cell (j)
        # via 4 endpoint-pair paths:
        #   (from_i → from_j), (from_i → to_j),
        #   (to_i  → from_j), (to_i  → to_j)
        nd_ff = nd_from[start:end][:, from_idx]  # (chunk, n): dist(from_i, from_j)
        nd_ft = nd_from[start:end][:, to_idx]    # (chunk, n): dist(from_i, to_j)
        nd_tf = nd_to  [start:end][:, from_idx]  # (chunk, n): dist(to_i, from_j)
        nd_tt = nd_to  [start:end][:, to_idx]    # (chunk, n): dist(to_i, to_j)

        dtf_i = dist_to_from[start:end, None]  # (chunk, 1)
        dtt_i = dist_to_to  [start:end, None]  # (chunk, 1)
        dtf_j = dist_to_from[None, :]          # (1, n)
        dtt_j = dist_to_to  [None, :]          # (1, n)

        chunk_dist = np.minimum.reduce([
            dtf_i + nd_ff + dtf_j,
            dtf_i + nd_ft + dtt_j,
            dtt_i + nd_tf + dtf_j,
            dtt_i + nd_tt + dtt_j,
        ]).astype(np.float32)  # (chunk, n)

        # Override same-edge pairs with exact intra-edge distance
        fi_c = from_idx[start:end]
        ti_c = to_idx[start:end]
        pi_c = pos[start:end]
        li_c = edge_len[start:end]

        same_fwd = (fi_c[:, None] == from_idx[None, :]) & (ti_c[:, None] == to_idx[None, :])
        same_rev = (fi_c[:, None] == to_idx[None, :])   & (ti_c[:, None] == from_idx[None, :])

        chunk_dist[same_fwd] = (
            np.abs(pi_c[:, None] - pos[None, :]) * li_c[:, None]
        ).astype(np.float32)[same_fwd]
        chunk_dist[same_rev] = (
            np.abs(pi_c[:, None] - (1.0 - pos[None, :])) * li_c[:, None]
        ).astype(np.float32)[same_rev]

        # Zero out self-distances
        local_self = np.arange(end - start)
        chunk_dist[local_self, local_self + start] = 0.0

        dist_mat[start:end] = chunk_dist

    # Symmetrise (pmax in R)
    return np.fmax(dist_mat, dist_mat.T)


def compute_trajectory_knn_weights(
    dist_mat: np.ndarray,
    K: Optional[int] = None,
) -> scipy.sparse.csr_matrix:
    """Gaussian-kernel KNN weights from trajectory geodesic distances.

    Mirrors R VISION's ``computeKNNWeights(object = Trajectory, K)``.

    The kernel is ``exp(-d² / σ²)`` where σ = max distance to the K-th
    nearest neighbour (per row), row-normalised to sum to 1.
    """
    n = dist_mat.shape[0]
    if K is None:
        K = max(1, round(np.sqrt(n)))

    d = dist_mat.astype(float).copy()
    np.fill_diagonal(d, np.inf)

    # Indices of K nearest neighbours per row
    knn_idx  = np.argpartition(d, K, axis=1)[:, :K]   # (n, K)
    rows     = np.arange(n)[:, None]
    knn_dist = d[rows, knn_idx]                        # (n, K)

    sigma = knn_dist.max(axis=1, keepdims=True)
    sigma[sigma == 0] = 1.0

    w = np.exp(-knn_dist ** 2 / sigma ** 2)
    w /= w.sum(axis=1, keepdims=True).clip(min=1e-12)

    I = np.repeat(np.arange(n), K)
    J = knn_idx.ravel()
    return scipy.sparse.csr_matrix((w.ravel(), (I, J)), shape=(n, n))


def _translate_cell_positions(
    progressions: pd.DataFrame,
    v_data: np.ndarray,
    milestone_ids: List[str],
    adj: np.ndarray,
) -> np.ndarray:
    """Map fractional edge positions into 2-D layout coordinates.

    Mirrors R VISION's ``translateCellPositions``.
    """
    id_to_idx = {m: i for i, m in enumerate(milestone_ids)}
    n = len(progressions)
    p_data = np.zeros((n, 2), dtype=float)

    fv = progressions["from"].values
    tv = progressions["to"].values
    pv = progressions["position"].values

    for mi, mj in _get_edges(adj):
        m_from, m_to = milestone_ids[mi], milestone_ids[mj]
        fc, tc = v_data[mi], v_data[mj]
        edge_dist = float(np.linalg.norm(tc - fc)) or 1.0

        mask_fwd = (fv == m_from) & (tv == m_to)
        mask_rev = (fv == m_to)   & (tv == m_from)
        all_idx  = np.where(mask_fwd | mask_rev)[0]
        if len(all_idx) == 0:
            continue

        # Position along edge (0 = m_from side)
        all_pos = np.where(mask_fwd[all_idx], pv[all_idx], 1.0 - pv[all_idx])

        x = all_pos * edge_dist
        y = np.random.normal(0.0, 0.1, len(all_idx))

        # Rotation matrix aligning (1, 0) to edge direction
        dx, dy = tc[0] - fc[0], tc[1] - fc[1]
        c, s = dx / edge_dist, dy / edge_dist
        # (local @ R^T) + offset   where R = [[c, -s], [s, c]]
        R_T = np.array([[c, s], [-s, c]])
        local = np.column_stack([x, y])
        p_data[all_idx] = local @ R_T + fc

    return p_data


def generate_trajectory_projections(
    progressions: pd.DataFrame,
    adj: np.ndarray,
    milestone_ids: List[str],
) -> Dict[str, Dict]:
    """Compute 2-D trajectory projections for the VISION UI.

    Tries FR, Tree, MDS, and DH (in that order).  Each projection is a dict
    with keys ``p_data`` (cell coords), ``v_data`` (milestone coords),
    ``adj_binary``, ``cell_ids``, ``milestone_ids``.
    """
    adj_binary = (adj > 0).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        adj_inv = np.where(adj > 0, 1.0 / adj, 0.0)

    n_m     = len(milestone_ids)
    edges   = _get_edges(adj)
    weights = [float(adj[i, j]) for i, j in edges]
    inv_w   = [float(adj_inv[i, j]) for i, j in edges]

    def _graph(w):
        g = ig.Graph(n=n_m, edges=edges)
        g.es["weight"] = w
        g.vs["name"]   = milestone_ids
        return g

    g      = _graph(weights)
    g_inv  = _graph(inv_w)

    layout_specs = [
        ("FR",   g_inv, lambda gr: gr.layout_fruchterman_reingold(weights=gr.es["weight"])),
        ("Tree", g,     lambda gr: gr.layout_reingold_tilford()),
        ("MDS",  g,     lambda gr: gr.layout_mds()),
        ("DH",   g,     lambda gr: gr.layout_davidson_harel()),
    ]

    projections: Dict[str, Dict] = {}
    cell_ids = progressions.index.tolist()

    for name, graph, layout_fn in layout_specs:
        try:
            layout  = layout_fn(graph)
            v_data  = np.array(layout.coords)  # (n_m, 2)
            p_data  = _translate_cell_positions(progressions, v_data, milestone_ids, adj)
            projections[name] = {
                "p_data":       p_data.tolist(),
                "v_data":       v_data.tolist(),
                "adj_binary":   adj_binary.tolist(),
                "cell_ids":     cell_ids,
                "milestone_ids": milestone_ids,
            }
        except Exception:
            pass  # skip layouts unavailable in this igraph build

    return projections


def create_trajectory_metadata(
    progressions: pd.DataFrame,
    adj: np.ndarray,
    milestone_ids: List[str],
) -> pd.DataFrame:
    """Build per-cell trajectory metadata (edge label + fractional position).

    Mirrors R VISION's ``createTrajectoryMetaData``.
    """
    meta = pd.DataFrame(index=progressions.index)
    meta["TrajectoryEdge"] = ""

    fv = progressions["from"].values
    tv = progressions["to"].values
    pv = progressions["position"].values

    for mi, mj in _get_edges(adj):
        m_from, m_to = milestone_ids[mi], milestone_ids[mj]
        edge_name    = f"{m_from}->{m_to}"
        pos_col      = f"Position: {edge_name}"

        mask_fwd = (fv == m_from) & (tv == m_to)
        mask_rev = (fv == m_to)   & (tv == m_from)

        meta.loc[mask_fwd, "TrajectoryEdge"] = edge_name
        meta.loc[mask_rev, "TrajectoryEdge"] = edge_name

        meta[pos_col] = np.nan
        meta.loc[mask_fwd, pos_col] = pv[mask_fwd]
        meta.loc[mask_rev, pos_col] = 1.0 - pv[mask_rev]

    meta["TrajectoryEdge"] = meta["TrajectoryEdge"].astype("category")
    return meta
