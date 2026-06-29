"""PhyloVision utilities: Fitch-Hartigan parsimony and plasticity scores.

Mirrors R VISION's ``computePlasticityScores`` /
``computeFitchHartiganParsimonyPerNode`` from ``FitchParsimony.R`` and
``AnalysisFunctions.R``.
"""

from __future__ import annotations

import logging
from collections import Counter
from typing import Dict, FrozenSet, List, Optional

import numpy as np
import pandas as pd
from anndata import AnnData

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Fitch-Hartigan algorithm helpers
# ---------------------------------------------------------------------------

def _build_parent_map(root) -> Dict[int, object]:
    """Return {id(clade): parent_clade} for every node in the tree."""
    parent: Dict[int, object] = {}

    def _dfs(clade, par) -> None:
        parent[id(clade)] = par
        for child in clade.clades:
            _dfs(child, clade)

    _dfs(root, None)
    return parent


def _fitch_bottom_up(
    root,
    metadata: Dict[str, str],
) -> Dict[int, FrozenSet[str]]:
    """Bottom-up Fitch-Hartigan pass for the subtree rooted at *root*.

    For each leaf: possible = {metadata[leaf]}.
    For each internal node: possible = mode of children's possible sets
    (union of most-frequent labels across children).

    Matches R VISION's ``bottomUpFitchHartigan``.
    """
    possible: Dict[int, FrozenSet[str]] = {}

    for node in root.find_clades(order="postorder"):
        if node.is_terminal():
            label = metadata.get(node.name)
            possible[id(node)] = frozenset([label] if label is not None else [])
        else:
            counter: Counter = Counter()
            for child in node.clades:
                for label in possible.get(id(child), frozenset()):
                    counter[label] += 1
            if counter:
                max_freq = max(counter.values())
                possible[id(node)] = frozenset(
                    lbl for lbl, freq in counter.items() if freq == max_freq
                )
            else:
                possible[id(node)] = frozenset()

    return possible


def _fitch_top_down(
    root,
    possible: Dict[int, FrozenSet[str]],
) -> Dict[int, Optional[str]]:
    """Top-down Fitch-Hartigan pass; returns a maximum-parsimony assignment.

    Tie-breaking is deterministic (lexicographic minimum), unlike R VISION
    which uses ``sample()`` (random).  The parsimony score is identical for
    any valid max-parsimony assignment.

    Matches R VISION's ``topDownFitchHartigan``.
    """
    assignments: Dict[int, Optional[str]] = {}

    def _assign(node, parent_label: Optional[str]) -> None:
        poss = possible.get(id(node), frozenset())
        if not poss:
            assignments[id(node)] = parent_label
        elif parent_label is not None and parent_label in poss:
            assignments[id(node)] = parent_label
        else:
            assignments[id(node)] = min(poss)          # deterministic tie-break
        for child in node.clades:
            _assign(child, assignments[id(node)])

    _assign(root, None)
    return assignments


def _score_parsimony(root, assignments: Dict[int, Optional[str]]) -> int:
    """Count label transitions across all edges in the subtree.

    Matches R VISION's ``scoreParsimony``.
    """
    parsimony = 0
    for node in root.find_clades():
        par_label = assignments.get(id(node))
        for child in node.clades:
            if par_label != assignments.get(id(child)):
                parsimony += 1
    return parsimony


def _normalized_parsimony(root, metadata: Dict[str, str]) -> float:
    """Normalized Fitch-Hartigan parsimony for the subtree rooted at *root*.

    = parsimony / (n_nodes - 1), where n_nodes is the number of nodes
    in the subtree (leaves + internal).

    Matches R VISION's ``computeNormalizedFitchHartiganParsimony(tree,
    metaData, source=node)``.
    """
    possible = _fitch_bottom_up(root, metadata)
    assignments = _fitch_top_down(root, possible)
    parsimony = _score_parsimony(root, assignments)
    n_nodes = sum(1 for _ in root.find_clades())
    return parsimony / (n_nodes - 1) if n_nodes > 1 else 0.0


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_plasticity_scores(
    adata: AnnData,
    newick: str,
) -> None:
    """Compute per-cell Fitch-Hartigan parsimony plasticity scores.

    For every categorical variable in ``adata.obs``, a ``{var}_plasticity``
    column is added to ``adata.obs``.  The score for cell *c* is the average
    normalised Fitch-Hartigan parsimony across all internal nodes on the path
    from the tree root to *c*.

    - **High plasticity**: the metadata label switches often along the
      lineage path → the trait is labile / plastic.
    - **Low plasticity**: the label is conserved along the lineage path.

    Mirrors R VISION's ``computePlasticityScores()``.

    Parameters
    ----------
    adata : AnnData
        Must have ``obs_names`` matching tree leaf labels.
    newick : str
        Newick-format tree string (or the ``adata.uns["vision_tree"]`` value).
    """
    from Bio import Phylo
    import io

    tree = Phylo.read(io.StringIO(newick), "newick")
    terminals = tree.get_terminals()

    obs_names_set = set(adata.obs_names)
    name_to_term = {t.name: t for t in terminals}
    parent = _build_parent_map(tree.root)

    # Identify categorical / object columns
    cat_cols: List[str] = [
        col
        for col in adata.obs.columns
        if pd.api.types.is_categorical_dtype(adata.obs[col])
        or adata.obs[col].dtype == object
    ]

    if not cat_cols:
        logger.info("No categorical obs columns found; skipping plasticity scores.")
        return

    # Collect all internal nodes (needed for per-node parsimony computation)
    internal_nodes = [c for c in tree.find_clades() if not c.is_terminal()]
    logger.info(
        "Computing plasticity scores for %d categorical variable(s) "
        "across %d internal nodes.",
        len(cat_cols),
        len(internal_nodes),
    )

    for col in cat_cols:
        # Build metadata dict: leaf name → label string
        metadata: Dict[str, str] = {
            name: str(adata.obs.at[name, col])
            for name in adata.obs_names
            if name in name_to_term
        }
        if not metadata:
            continue

        # Compute normalised parsimony score for each internal node's subtree
        node_score: Dict[int, float] = {
            id(node): _normalized_parsimony(node, metadata)
            for node in internal_nodes
        }

        # Per-leaf plasticity = average node score along root-to-leaf path
        plasticity: Dict[str, float] = {}
        for leaf_name in adata.obs_names:
            if leaf_name not in name_to_term:
                plasticity[leaf_name] = np.nan
                continue

            leaf = name_to_term[leaf_name]
            ancestor_scores: List[float] = []

            # Walk from leaf's parent up to (and including) root
            node = leaf
            while parent[id(node)] is not None:
                node = parent[id(node)]
                s = node_score.get(id(node))
                if s is not None:
                    ancestor_scores.append(s)

            plasticity[leaf_name] = (
                float(np.mean(ancestor_scores)) if ancestor_scores else 0.0
            )

        out_col = f"{col}_plasticity"
        adata.obs[out_col] = [plasticity.get(name, np.nan) for name in adata.obs_names]
        logger.info("Plasticity scores stored in adata.obs['%s'].", out_col)
