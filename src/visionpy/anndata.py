import math
from typing import Literal, Optional, Union

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from scipy.sparse import issparse
from scipy.stats import chi2_contingency

from .diffexp import rank_genes_groups
from .signature import (
    compute_obs_df_scores,
    compute_signature_scores,
    _gearysc_for_dataframe,
)


def _pearson_cols(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Pearson correlation between each column of A (n×p) and vector b (n,)."""
    A = A - A.mean(axis=0)        # (n, p) centered
    b = b - b.mean()              # (n,)  centered
    denom_A = np.sqrt((A ** 2).sum(axis=0))   # (p,)
    denom_b = float(np.sqrt((b ** 2).sum()))
    denom_A[denom_A == 0] = 1.0
    if denom_b == 0:
        return np.zeros(A.shape[1])
    return (b @ A) / (denom_b * denom_A)      # (p,)


class AnnDataAccessor:
    def __init__(
        self,
        adata: Optional[anndata.AnnData] = None,
        norm_data_key: Optional[Union[Literal["use_raw"], str]] = None,
        protein_obsm_key: Optional[str] = None,
        signature_varm_key: Optional[str] = None,
        signature_names_uns_key: Optional[str] = None,
    ) -> None:
        """Accessor to AnnData object.

        Parameters
        ----------
        adata
            AnnData object
        norm_data_key
            Key for layer with log library size normalized data. If
            `None` (default), uses `adata.X`
        protein_obsm_key
            Location for protein data
        signature_varm_key
            Location for genes by signature matrix
        """
        self._adata = adata
        self._protein_obsm_key = protein_obsm_key
        self._norm_data_key = norm_data_key
        self._signature_varm_key = signature_varm_key
        self._signature_names_uns_key = signature_names_uns_key
        self._cells_selections = {}
        self.protein_adata = None

    @property
    def adata(self):
        return self._adata

    @property
    def var_names(self):
        if self._norm_data_key == "use_raw":
            return self.adata.raw.var_names
        else:
            return self.adata.var_names

    @property
    def cells_selections(self):
        return self._cells_selections.keys()

    def add_cells_selection(self, key, val):
        self._cells_selections[key] = val

    def get_cells_selection(self, key):
        return self._cells_selections[key]

    @adata.setter
    def adata(self, adata: anndata.AnnData):
        self._adata = adata
        num_cols = adata.obs._get_numeric_data().columns.tolist()
        cols = adata.obs.columns.tolist()
        cat_vars = list(set(cols) - set(num_cols))
        self.cat_obs_cols = cat_vars
        self.numeric_obs_cols = num_cols

    @property
    def protein_obsm_key(self):
        return self._protein_obsm_key

    @protein_obsm_key.setter
    def protein_obsm_key(self, key: str):
        self._protein_obsm_key = key

    @property
    def signature_varm_key(self):
        return self._signature_varm_key

    @signature_varm_key.setter
    def signature_varm_key(self, key: str):
        self._signature_varm_key = key

    @property
    def signature_names_uns_key(self):
        return self._signature_names_uns_key

    @signature_names_uns_key.setter
    def signature_names_uns_key(self, key: str):
        self._signature_names_uns_key = key

    @property
    def norm_data_key(self):
        return self._norm_data_key

    @norm_data_key.setter
    def norm_data_key(self, key: str):
        self._norm_data_key = key

    def get_gene_expression(self, gene: str, return_list=True) -> list:
        if self.adata is None:
            raise ValueError("Accessor not populated with anndata.")
        if self.norm_data_key == "use_raw":
            data = self.adata.raw[:, gene].X
        elif self.norm_data_key is None:
            data = self.adata[:, gene].X
        else:
            data = self.adata[:, gene].layers[self.norm_data_key]

        if scipy.sparse.issparse(data):
            data = data.toarray()

        if return_list:
            return data.ravel().tolist()
        else:
            return data

    def get_gene_expression_raw(self, gene: str) -> list:
        """Return unnormalized (raw count) expression for *gene*.

        Uses the layer stored in ``adata.uns["vision_unnormalized_key"]`` when
        available, otherwise falls back to ``adata.raw.X`` or ``adata.X``.
        """
        if self.adata is None:
            raise ValueError("Accessor not populated with anndata.")
        unnorm_key = self.adata.uns.get("vision_unnormalized_key")
        if unnorm_key is not None and unnorm_key in self.adata.layers:
            data = self.adata[:, gene].layers[unnorm_key]
        elif self.adata.raw is not None:
            data = self.adata.raw[:, gene].X
        else:
            data = self.adata[:, gene].X
        if scipy.sparse.issparse(data):
            data = data.toarray()
        return data.ravel().tolist()

    def get_genes_by_signature(self, sig_name: str) -> pd.DataFrame:
        """Df of genes in index, sign as values."""
        
        if self.signature_names_uns_key is not None:
            index = np.where(np.asarray(self.adata.uns[self.signature_names_uns_key]) == sig_name)[0][0]
        else:
            index = np.where(np.asarray(self.adata.obsm["vision_signatures"].columns) == sig_name)[0][0]
        
        if self._norm_data_key == "use_raw":
            matrix = self.adata.raw.varm[self.signature_varm_key]
        else:
            matrix = self.adata.varm[self.signature_varm_key]

        if isinstance(matrix, pd.DataFrame):
            matrix = matrix.to_numpy()

        matrix = matrix[:, index]
        if issparse(matrix):
            matrix = matrix.toarray().ravel()

        mask = matrix != 0
        sig_df = pd.DataFrame(index=self.var_names[mask], data=matrix[mask])

        return sig_df

    def compute_obs_df_scores(self):
        self.adata.uns["vision_obs_df_scores"] = compute_obs_df_scores(self.adata)

    def compute_protein_autocorrelation(self):
        """Geary's C for CITE-seq protein data.

        Mirrors R VISION's ``fbConsistencyScores``: rank-transforms protein
        columns and computes Geary's C using the cell KNN weight graph.
        No permutation p-values (matching R VISION's ``computePval=FALSE``).
        Results stored in ``adata.uns["vision_protein_autocorr"]``.
        """
        if self._protein_obsm_key is None:
            return
        mat = self.adata.obsm[self._protein_obsm_key]
        if isinstance(mat, pd.DataFrame):
            protein_df = mat
        else:
            n_prot = mat.shape[1]
            cols = (
                [f"Protein_{i}" for i in range(n_prot)]
                if not hasattr(mat, "columns")
                else mat.columns.tolist()
            )
            protein_df = pd.DataFrame(
                np.asarray(mat, dtype=float),
                index=self.adata.obs_names,
                columns=cols,
            )
        weights = self.adata.obsp["weights"].tocsr()
        result = _gearysc_for_dataframe(weights, protein_df, compute_pvals=False)
        self.adata.uns["vision_protein_autocorr"] = result

    def compute_protein_differential(self):
        """One-vs-all Wilcoxon differential for proteins across cluster levels.

        Mirrors the signature differential pipeline applied to protein expression.
        Results stored in ``adata.uns["vision_protein_differential"]``.
        """
        if self._protein_obsm_key is None:
            return
        mat = self.adata.obsm[self._protein_obsm_key]
        if isinstance(mat, pd.DataFrame):
            prot_arr = mat.to_numpy()
            prot_names = mat.columns.tolist()
        else:
            prot_arr = np.asarray(mat, dtype=float)
            prot_names = [f"Protein_{i}" for i in range(prot_arr.shape[1])]

        prot_adata = anndata.AnnData(
            prot_arr,
            obs=self.adata.obs.loc[:, self.cat_obs_cols].copy(),
        )
        prot_adata.var_names = prot_names

        for c in self.cat_obs_cols:
            try:
                rank_genes_groups(
                    prot_adata,
                    groupby=c,
                    key_added=f"rank_genes_groups_{c}",
                    method="wilcoxon",
                )
            except ValueError:
                continue
        self.protein_adata = prot_adata
        self.adata.uns["vision_protein_differential_key"] = self._protein_obsm_key

    def compute_signature_scores(self):
        self.adata.uns["vision_signature_scores"] = compute_signature_scores(
            self.adata, self.norm_data_key, self.signature_varm_key
        )

    def compute_one_vs_all_signatures(self):
        sig_adata = anndata.AnnData(self.adata.obsm["vision_signatures"])
        sig_adata.obs = self.adata.obs.loc[:, self.cat_obs_cols].copy()
        for c in self.cat_obs_cols:
            rank_genes_groups(
                sig_adata,
                groupby=c,
                key_added=f"rank_genes_groups_{c}",
                method="wilcoxon",
            )
        self.sig_adata = sig_adata

    def compute_one_vs_all_obs_cols(self):
        numeric_data = self.adata.obs._get_numeric_data().copy()
        n_numeric = numeric_data.shape[1]
        # AnnData requires ≥1 variable; use a zeros placeholder when there are
        # no numeric metadata columns so categorical chi-squared still runs.
        obs_X = np.log1p(numeric_data.to_numpy()) if n_numeric > 0 else np.zeros((self.adata.n_obs, 1))
        obs_var_names = list(numeric_data.columns) if n_numeric > 0 else ["_placeholder"]
        obs_adata = anndata.AnnData(
            X=obs_X,
            var=pd.DataFrame(index=obs_var_names),
            obs=self.adata.obs.loc[:, self.cat_obs_cols].copy(),
        )
        for c in self.cat_obs_cols:
            try:
                if n_numeric > 0:
                    rank_genes_groups(
                        obs_adata,
                        groupby=c,
                        key_added=f"rank_genes_groups_{c}",
                        method="wilcoxon",
                    )
            # one category only has one obs
            except ValueError:
                # TODO: Log it
                self.cat_obs_cols = [c_ for c_ in self.cat_obs_cols if c_ != c]
                continue

            for g in categories(obs_adata.obs[c]):
                mask = (obs_adata.obs[c] == g).to_numpy()
                obs_pos_masked = obs_adata.obs.iloc[mask]
                obs_neg_masked = obs_adata.obs.iloc[~mask]
                for j in obs_pos_masked.columns:
                    pos_freq = obs_pos_masked[j].value_counts(normalize=False)
                    neg_freq = obs_neg_masked[j].value_counts(normalize=False)
                    freqs = pd.concat([pos_freq, neg_freq], axis=1).fillna(0)
                    # TODO: cramer's v might be incorrect
                    grand_total = np.sum(freqs.to_numpy())
                    r = len(freqs) - 1
                    try:
                        # freqs shape: (n_categories, 2) — rows=categories, cols=pos/neg group
                        chi2, pval, _, _ = chi2_contingency(freqs.to_numpy())
                        stat = chi2
                    except ValueError:
                        stat = grand_total  # so that v is 1
                        pval = 0
                    if math.isinf(pval) or math.isnan(pval):
                        pval = 1
                    if math.isinf(stat) or math.isnan(stat):
                        v = 1
                    else:
                        # Cramér's V for a 2-column table: V = sqrt(chi2 / n)
                        # min(k-1, r-1) = min(2-1, r-1) = 1 for r ≥ 2
                        v = np.sqrt(stat / grand_total)
                    obs_adata.uns[f"chi_sq_{j}_{g}"] = {
                        "stat": v,
                        "pval": pval,
                    }

        self.obs_adata = obs_adata

    def compute_gene_score_per_signature(self):
        gene_score_sig = {}
                
        if self.signature_names_uns_key is not None:
            sig_names = self.adata.uns[self.signature_names_uns_key]
        else:
            sig_names = self.adata.obsm["vision_signatures"].columns
        
        for s in sig_names:
            gene_score_sig[s] = {"genes": [], "values": []}
            df = self.get_genes_by_signature(s)
            gene_names = df.index
            # cells by genes
            expr = np.array(self.get_gene_expression(gene_names, return_list=False))
            # cells
            sign = df.to_numpy().ravel()
            gene_score_sig[s]["signs"] = sign.tolist()

            # Vectorized: correlate all signature genes at once
            sig_scores = self.adata.obsm["vision_signatures"][s].to_numpy().ravel()
            corrs = _pearson_cols(expr, sig_scores)   # (n_genes,)
            gene_score_sig[s]["values"] = (sign * corrs).tolist()
            gene_score_sig[s]["genes"] = gene_names.tolist()

        for s in sig_names:
            info = gene_score_sig[s]
            gene_score_sig[s]["geneImportance"] = {g: v for g, v in zip(info["genes"], info["values"])}
            gene_score_sig[s]["sigDict"] = {g: v for g, v in zip(info["genes"], info["signs"])}

        self.gene_score_sig = gene_score_sig

    def persist_signature_differential(self):
        """Persist one-vs-all Wilcoxon signature results to adata.uns.

        Writes ``adata.uns["vision_signature_differential"]`` as a dict keyed
        by the groupby column name.  Each value is a tidy DataFrame with columns
        ``["group", "names", "scores", "pvals", "pvals_adj", "logfoldchanges"]``
        where "names" are signature names.  Mirrors R VISION's
        ``@ClusterComparisons$Signatures``.
        """
        if not hasattr(self, "sig_adata"):
            return
        result = {}
        for c in self.cat_obs_cols:
            key = f"rank_genes_groups_{c}"
            if key not in self.sig_adata.uns:
                continue
            groups = list(self.sig_adata.obs[c].astype("category").cat.categories)
            frames = []
            for g in groups:
                try:
                    df = sc.get.rank_genes_groups_df(self.sig_adata, group=str(g), key=key)
                    df.insert(0, "group", str(g))
                    frames.append(df)
                except Exception:
                    continue
            if frames:
                result[c] = pd.concat(frames, ignore_index=True)
        if result:
            self.adata.uns["vision_signature_differential"] = result

    def persist_meta_differential(self):
        """Persist one-vs-all metadata differential results to adata.uns.

        Writes ``adata.uns["vision_meta_differential"]`` as a nested dict:
        ``{groupby_col: {"numeric": DataFrame, "categorical": DataFrame}}``.

        * ``numeric`` — Wilcoxon test for each numeric obs column grouped by
          the categorical variable.  Columns: ``["group", "names", "scores",
          "pvals", "pvals_adj", "logfoldchanges"]``.
        * ``categorical`` — Cramér's V between each pair of categorical obs
          columns.  Columns: ``["group", "names", "cramers_v", "pval"]``.

        Mirrors R VISION's ``@ClusterComparisons$Meta``.
        """
        if not hasattr(self, "obs_adata"):
            return
        result = {}
        for c in self.cat_obs_cols:
            col_result = {}
            groups = list(self.obs_adata.obs[c].astype("category").cat.categories)

            # ── Numeric metadata (Wilcoxon) ───────────────────────────────────
            key = f"rank_genes_groups_{c}"
            if key in self.obs_adata.uns and self.obs_adata.n_vars > 0:
                frames = []
                for g in groups:
                    try:
                        df = sc.get.rank_genes_groups_df(
                            self.obs_adata, group=str(g), key=key
                        )
                        df.insert(0, "group", str(g))
                        frames.append(df)
                    except Exception:
                        continue
                if frames:
                    col_result["numeric"] = pd.concat(frames, ignore_index=True)

            # ── Categorical pairs (Cramér's V) ────────────────────────────────
            cat_rows = []
            for g in groups:
                for j in self.cat_obs_cols:
                    chi_key = f"chi_sq_{j}_{g}"
                    if chi_key in self.obs_adata.uns:
                        info = self.obs_adata.uns[chi_key]
                        cat_rows.append(
                            {"group": str(g), "names": str(j),
                             "cramers_v": info["stat"], "pval": info["pval"]}
                        )
            if cat_rows:
                col_result["categorical"] = pd.DataFrame(cat_rows)

            if col_result:
                result[c] = col_result
        if result:
            self.adata.uns["vision_meta_differential"] = result

    def persist_gene_importance(self):
        """Persist per-signature gene importance scores to adata.uns.

        Writes ``adata.uns["vision_gene_importance"]`` as a dict keyed by
        signature name.  Each value is a DataFrame indexed by gene name with
        columns ``["importance", "sign"]`` where
        ``importance = sign × Pearson(gene_expr, sig_score)``.
        Mirrors R VISION's ``@SigGeneImportance``.
        """
        if not hasattr(self, "gene_score_sig") or not self.gene_score_sig:
            return
        result = {}
        for sig_name, info in self.gene_score_sig.items():
            if info["genes"]:
                result[sig_name] = pd.DataFrame(
                    {"importance": info["values"], "sign": info["signs"]},
                    index=info["genes"],
                )
        if result:
            self.adata.uns["vision_gene_importance"] = result

    def compute_one_vs_one_de(self, key: str, group1: str, group2: str):
        rank_genes_groups(
            self.adata,
            groupby=key,
            groups=[group1],
            reference=group2,
            key_added=f"rank_genes_groups_{key}",
            method="wilcoxon",
            use_raw=self.norm_data_key == "use_raw",
            layer=self.norm_data_key if self.norm_data_key != "use_raw" else None,
        )
        return sc.get.rank_genes_groups_df(self.adata, group1, key=f"rank_genes_groups_{key}")


def categories(col):
    return col.astype("category").cat.categories
