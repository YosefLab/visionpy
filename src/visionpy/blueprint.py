import json
from collections import OrderedDict
from dataclasses import dataclass
from typing import List

import numpy as np
import pandas as pd
import scanpy as sc
from flask import Blueprint, jsonify, render_template, request
from scipy.stats import entropy

from visionpy import data_accessor

bp = Blueprint("api", __name__)

adata = data_accessor.adata


@dataclass
class ServerSigProjMatrix:
    zscores: List[float]
    pvals: List[float]
    proj_labels: List[str]
    sig_labels: List[str]

    def prepare_json(self):
        return {
            "class": ["ServerSigProjMatrix"],
            "sig_labels": self.sig_labels,
            "proj_labels": self.proj_labels,
            "zscores": self.zscores,
            "pvals": self.pvals,
        }


@bp.route("/")
def app():
    return render_template("Results.html")


@bp.route("/Projections/list", methods=["GET"])
def get_projection_list():
    """Returns dict of col names for each projection."""
    proj_dict = OrderedDict()
    # get only numeric columns
    proj_dict["Obs_metadata"] = data_accessor.numeric_obs_cols
    for k in adata.obsm_keys():
        if k[:2] == "X_":
            name = k.split("_")[1]
            proj_dict[k] = [name.upper() + f"{i + 1}" for i in range(adata.obsm[k].shape[1])]

    return jsonify(proj_dict)


@bp.route("/Projections/<proj_name>/coordinates/<proj_col>", methods=["GET"])
def get_projection_column(proj_name, proj_col):
    if proj_name == "Obs_metadata":
        data = adata.obs[proj_col].values.tolist()
    else:
        name_prefix = proj_name.split("_", 1)[1].upper()
        column_ind = int(proj_col[len(name_prefix):]) - 1
        data = adata.obsm[proj_name]
        if isinstance(data, pd.DataFrame):
            data = data.iloc[:, column_ind].values.tolist()
        else:
            data = data[:, column_ind].tolist()

    # the order matters here
    df = pd.DataFrame()
    df["Coords"] = data
    df["Barcodes"] = adata.obs_names.tolist()

    # need a list of lists, where each internal list is a coord value
    # and then a barcode value
    # coord value needs to be float, barcode str
    return jsonify(df.to_numpy().tolist())


@bp.route("/Tree/Projections/list", methods=["GET"])
def get_tree_projection_list():
    projs = adata.uns.get("vision_trajectory_projections", {})
    return jsonify(list(projs.keys()))


@bp.route("/Tree/Projections/<proj_name>/coordinates", methods=["GET"])
def get_tree_projection_coordinates(proj_name):
    projs = adata.uns.get("vision_trajectory_projections", {})
    if proj_name not in projs:
        return jsonify([]), 404

    proj   = projs[proj_name]
    p_data = np.array(proj["p_data"])     # (n_cells, 2)
    v_data = np.array(proj["v_data"])     # (n_milestones, 2)
    adj_bin = np.array(proj["adj_binary"])  # (n_milestones, n_milestones)
    cell_ids = proj["cell_ids"]

    # Format expected by the VISION frontend:
    # [ [[x, y, cell_id], ...], [[x, y], ...], [[0/1, ...], ...] ]
    coords = [
        [float(p_data[i, 0]), float(p_data[i, 1]), cell_ids[i]]
        for i in range(len(cell_ids))
    ]
    return jsonify([coords, v_data.tolist(), adj_bin.tolist()])


@bp.route("/Tree/SigProjMatrix/Normal", methods=["GET"])
def get_tree_sigprojmatrix_normal():
    traj_ac = adata.uns.get("vision_trajectory_autocorr", {})
    sig_labels, stats, pvals = [], [], []

    if "Signatures" in traj_ac:
        df = traj_ac["Signatures"]
        sig_labels += df.index.tolist()
        stats += [[v] for v in df["c_prime"].tolist()]
        pvals += [[v] for v in df["fdr"].tolist()]

    if "Meta" in traj_ac:
        df = traj_ac["Meta"]
        num_only = df.loc[df.index.isin(data_accessor.numeric_obs_cols)]
        sig_labels += num_only.index.tolist()
        stats += [[v] for v in num_only["c_prime"].tolist()]
        pvals += [[v] for v in num_only["fdr"].tolist()]

    matrix = ServerSigProjMatrix(
        sig_labels=sig_labels, proj_labels=["Score"],
        zscores=stats, pvals=pvals,
    )
    return jsonify(matrix.prepare_json())


@bp.route("/Tree/SigProjMatrix/Meta", methods=["GET"])
def get_tree_sigprojmatrix_meta():
    traj_ac = adata.uns.get("vision_trajectory_autocorr", {})
    sig_labels, stats, pvals = [], [], []

    if "Meta" in traj_ac:
        df = traj_ac["Meta"]
        sig_labels = df.index.tolist()
        stats = [[v] for v in df["c_prime"].tolist()]
        pvals = [[v] for v in df["fdr"].tolist()]

    matrix = ServerSigProjMatrix(
        sig_labels=sig_labels, proj_labels=["Score"],
        zscores=stats, pvals=pvals,
    )
    return jsonify(matrix.prepare_json())


@bp.route("/Tree/ProteinMatrix", methods=["GET"])
def get_tree_protein_matrix():
    traj_ac = adata.uns.get("vision_trajectory_autocorr", {})
    sig_labels, stats, pvals = [], [], []

    if "Proteins" in traj_ac:
        df = traj_ac["Proteins"]
        sig_labels = df.index.tolist()
        stats = [[v] for v in df["c_prime"].tolist()]
        pvals = [[v] for v in df["fdr"].tolist()]

    matrix = ServerSigProjMatrix(
        sig_labels=sig_labels, proj_labels=["Score"],
        zscores=stats, pvals=pvals,
    )
    return jsonify(matrix.prepare_json())


@bp.route("/SessionInfo", methods=["GET"])
def get_session_info():
    info = {}
    info["name"] = adata.uns["vision_session_name"]
    info["has_tree"] = bool(adata.uns.get("vision_trajectory_projections"))
    info["meta_sigs"] = adata.obs.columns.tolist()
    # info[["pooled"]] <- object@params$micropooling$pool
    info["pooled"] = False
    info["ncells"] = adata.n_obs
    info["has_sigs"] = "vision_signatures" in adata.obsm_keys()
    # TODO: Fix
    info["has_proteins"] = "vision_protein_autocorr" in adata.uns
    info["has_lca"] = "vision_lca" in adata.uns
    # TODO: implement hotspot
    info["has_mods"] = False
    info["has_dendrogram"] = "vision_dendrogram" in adata.uns
    if info["has_dendrogram"]:
        info["dendrogram"] = adata.uns["vision_dendrogram"]
    # Do not change this
    info["has_seurat"] = False
    return jsonify(info)


@bp.route("/Clusters/MetaLevels", methods=["GET"])
def get_clusters_metalevels():
    cat_vars = data_accessor.cat_obs_cols

    cat_key_vals = {}
    for c in cat_vars:
        cat_key_vals[c] = list(adata.obs[c].astype("category").cat.categories)

    return jsonify(cat_key_vals)


@bp.route("/Clusters/<cluster_variable>/Cells", methods=["GET"])
def get_clusters_cluster_var_cells(cluster_variable):
    cats = adata.obs[cluster_variable].astype("category")
    result = {
        str(level): adata.obs_names[cats == level].tolist()
        for level in cats.cat.categories
    }
    return jsonify(result)


@bp.route("/Signature/Expression/<sig_name>/<cluster_var>", methods=["GET"])
def get_signature_expression(sig_name, cluster_var):
    sig_df = data_accessor.get_genes_by_signature(sig_name)
    gene_names = sig_df.index.tolist()

    # Expression matrix (n_cells, n_genes)
    expr = data_accessor.get_gene_expression(gene_names, return_list=False)
    if hasattr(expr, "toarray"):
        expr = expr.toarray()
    expr = np.asarray(expr, dtype=float)

    cats = adata.obs[cluster_var].astype("category")
    cluster_levels = cats.cat.categories.tolist()

    # Mean per cluster → (n_genes, n_clusters)
    mean_expr = np.zeros((len(gene_names), len(cluster_levels)))
    for j, level in enumerate(cluster_levels):
        mask = (cats == level).to_numpy()
        mean_expr[:, j] = expr[mask].mean(axis=0)

    # Z-score across clusters per gene
    row_mean = mean_expr.mean(axis=1, keepdims=True)
    row_std = mean_expr.std(axis=1, keepdims=True)
    row_std[row_std == 0] = 1.0
    z_matrix = (mean_expr - row_mean) / row_std

    return jsonify({
        "data": z_matrix.tolist(),
        "sample_labels": cluster_levels,
        "gene_labels": gene_names,
    })


@bp.route("/Cells/Selections", methods=["GET"])
def get_cells_selections():
    return jsonify(list(data_accessor.cells_selections))


@bp.route("/Cells/Selections/<selection_id>", methods=["GET", "POST"])
def get_cells_selections_sel_id(selection_id):
    if request.method == "GET":
        return jsonify(data_accessor.get_cells_selection(selection_id))
    else:
        subset = json.loads(list(dict(request.form.lists()).keys())[0])
        data_accessor.add_cells_selection(selection_id, subset)
        return jsonify([])


@bp.route("/Expression/Genes/List", methods=["GET"])
def get_gene_names():
    return jsonify(adata.var_names.tolist())


@bp.route("/Expression/Gene/<gene_name>", methods=["GET"])
def get_gene_expression(gene_name):
    return jsonify(
        dict(
            cells=adata.obs_names.tolist(),
            values=data_accessor.get_gene_expression(gene_name),
        )
    )


@bp.route("/Expression/Gene/<gene_name>/Raw", methods=["GET"])
def get_gene_expression_raw(gene_name):
    """Return raw (unnormalized) counts for a gene."""
    return jsonify(
        dict(
            cells=adata.obs_names.tolist(),
            values=data_accessor.get_gene_expression_raw(gene_name),
        )
    )


@bp.route("/Signature/Meta/<sig_name>", methods=["GET"])
def get_signature_meta(sig_name):
    return jsonify(dict(cells=adata.obs_names.tolist(), values=adata.obs.loc[:, sig_name].tolist()))


@bp.route("/Signature/Scores/<sig_name>", methods=["GET"])
def get_signature_score(sig_name):
    return jsonify(
        dict(
            cells=adata.obs_names.tolist(),
            values=adata.obsm["vision_signatures"][sig_name].to_numpy().tolist(),
        )
    )


@bp.route("/Signature/Info/<sig_name>", methods=["GET"])
def get_signature_info(sig_name):
    info = data_accessor.gene_score_sig[sig_name]
    return jsonify(
        dict(
            sigDict=info["sigDict"],
            geneImportance=info["geneImportance"],
            name=sig_name,
            source="user-defined",
            metaData="",
        )
    )


@bp.route("/FilterGroup/SigClusters/Normal", methods=["GET"])
def get_sigclusters_normal():
    if "vision_sig_clusters" in adata.uns:
        return jsonify(adata.uns["vision_sig_clusters"])
    cols = adata.obsm["vision_signatures"].columns.to_list()
    return jsonify({c: 1 for c in cols})


@bp.route("/FilterGroup/SigClusters/Meta", methods=["GET"])
def get_sigclusters_meta():
    cols = adata.obs.columns.to_list()
    clusters = {}
    for c in cols:
        clusters[c] = 1
    return jsonify(clusters)


@bp.route("/Clusters/<cluster_variable>/SigProjMatrix/Meta", methods=["GET"])
def get_sigprojmatrix_meta(cluster_variable):
    sig_labels = adata.obs.columns.tolist()
    # TODO: amortize computation in metalevels route
    proj_labels = ["Score"] + list(adata.obs[cluster_variable].astype("category").cat.categories)

    scores = adata.uns["vision_obs_df_scores"]["c_prime"].to_numpy().reshape(-1, 1)
    sigs_by_projs_stats = pd.DataFrame(index=sig_labels, columns=proj_labels[1:], data=0)
    sigs_by_projs_pvals = pd.DataFrame(index=sig_labels, columns=proj_labels[1:], data=0)
    obs_adata = data_accessor.obs_adata
    # TODO: test for categorical data with chi sq
    for p in proj_labels[1:]:
        temp_df = sc.get.rank_genes_groups_df(obs_adata, key=f"rank_genes_groups_{cluster_variable}", group=p)
        temp_df.set_index("names", inplace=True)
        sigs_by_projs_stats.loc[temp_df.index, p] = temp_df["scores"].copy()
        sigs_by_projs_pvals.loc[temp_df.index, p] = temp_df["pvals_adj"].copy()
        for cat_c in data_accessor.cat_obs_cols:
            key = f"chi_sq_{cluster_variable}_{p}"
            sigs_by_projs_stats.loc[cat_c, p] = obs_adata.uns[key]["stat"]
            sigs_by_projs_pvals.loc[cat_c, p] = obs_adata.uns[key]["pval"]

    stats = np.hstack([scores, sigs_by_projs_stats.to_numpy()]).tolist()
    pvals = np.hstack([np.zeros_like(scores), sigs_by_projs_pvals.to_numpy()]).tolist()

    matrix = ServerSigProjMatrix(sig_labels=sig_labels, proj_labels=proj_labels, zscores=stats, pvals=pvals)
    return jsonify(matrix.prepare_json())


@bp.route("/Clusters/<cluster_variable>/SigProjMatrix/Normal", methods=["GET"])
def get_sigprojmatrix_normal(cluster_variable):
    sig_labels = adata.obsm["vision_signatures"].columns.tolist()
    # TODO: amortize computation in metalevels route
    proj_labels = ["Score"] + list(adata.obs[cluster_variable].astype("category").cat.categories)

    # TODO: properly compute all information
    scores = adata.uns["vision_signature_scores"]["c_prime"].to_numpy().reshape(-1, 1)

    sigs_by_projs_stats = pd.DataFrame(index=sig_labels, columns=proj_labels[1:])
    sigs_by_projs_pvals = pd.DataFrame(index=sig_labels, columns=proj_labels[1:])
    sig_adata = data_accessor.sig_adata
    for p in proj_labels[1:]:
        temp_df = sc.get.rank_genes_groups_df(sig_adata, key=f"rank_genes_groups_{cluster_variable}", group=p)
        temp_df.set_index("names", inplace=True)
        temp_df = temp_df.loc[sigs_by_projs_stats.index]
        sigs_by_projs_stats[p] = temp_df["scores"].copy()
        sigs_by_projs_pvals[p] = temp_df["pvals_adj"].copy()

    stats = np.hstack([scores, sigs_by_projs_stats.to_numpy()]).tolist()
    pvals = np.hstack([np.zeros_like(scores), sigs_by_projs_pvals.to_numpy()]).tolist()

    matrix = ServerSigProjMatrix(sig_labels=sig_labels, proj_labels=proj_labels, zscores=stats, pvals=pvals)
    return jsonify(matrix.prepare_json())


@bp.route("/Cell/<cell_id>/Meta", methods=["GET", "POST"])
def get_cell_metadata(cell_id):
    if request.method == "GET":
        cell = adata.obs.loc[cell_id].to_dict()
        new_cell = {}
        for k, v in cell.items():
            if not isinstance(v, str):
                v_ = np.round(v, 4)
            else:
                v_ = v
            # doesn't work if starts with _
            if k[0] != "_":
                new_cell[k] = str(v_)
        return jsonify(new_cell)
    else:
        pass


@bp.route("/Cells/Meta", methods=["POST"])
def send_cells_meta():
    # TODO: Make nicer
    subset = json.loads(list(dict(request.form.lists()).keys())[0])
    df = adata[subset].obs
    numerical_df = df.loc[:, data_accessor.numeric_obs_cols]
    categorical_df = df.loc[:, data_accessor.cat_obs_cols]

    num_percentiles = np.percentile(numerical_df, [0, 50, 100], axis=0)
    numerical_stats = {}
    for i, c in enumerate(numerical_df.columns):
        numerical_stats[c] = {
            "Min": num_percentiles[0, i],
            "Median": num_percentiles[1, i],
            "Max": num_percentiles[2, i],
        }

    categorical_stats = {}
    for c in categorical_df.columns:
        dist = (categorical_df[c].value_counts(normalize=True) * 100).to_dict()
        categorical_stats[c] = {"levels": dist, "entropy": entropy(list(dist.values()), base=2)}

    return_dict = dict(numeric=numerical_stats, factor=categorical_stats)

    return jsonify(return_dict)


@bp.route("/DE", methods=["POST"])
def send_de():
    """Differential expression analysis."""
    # load request with json
    content = json.loads(list(request.form.to_dict().keys())[0])
    type_n = content["type_n"]
    type_d = content["type_d"]
    subtype_n = content["subtype_n"]
    subtype_d = content["subtype_d"]
    group_num = content["group_num"]
    group_denom = content["group_denom"]

    if type_n == "current":
        cell_ids_1 = subtype_n
    elif type_n == "saved_selection":
        cell_ids_1 = data_accessor.get_cells_selection(subtype_n)
    elif type_n == "meta":
        cell_ids_1 = adata.obs[adata.obs[subtype_n] == group_num].index.tolist()
    else:
        raise ValueError("type_n not recognized")

    if type_d == "current":
        cell_ids_2 = subtype_d
    elif type_d == "saved_selection":
        cell_ids_2 = data_accessor.get_cells_selection(subtype_d)
    elif type_d == "meta":
        cell_ids_2 = adata.obs[adata.obs[subtype_d] == group_denom].index.tolist()
    elif type_d == "remainder":
        cell_ids_2 = adata.obs.index.difference(cell_ids_1).tolist()
    else:
        raise ValueError("type_d not recognized")

    adata.obs["_vision_de_cache"] = "None"
    adata.obs.loc[cell_ids_1, "_vision_de_cache"] = "numerator"
    adata.obs.loc[cell_ids_2, "_vision_de_cache"] = "denominator"

    de_res = data_accessor.compute_one_vs_one_de("_vision_de_cache", "numerator", "denominator")

    send_de_res = dict(
        logFC=de_res["logfoldchanges"].tolist(),
        pval=de_res["pvals_adj"].tolist(),
        stat=de_res["scores"].tolist(),
        Feature=de_res["names"].tolist(),
        Type=["Gene"] * len(de_res),
        de_test_type="wilcox",
    )

    return jsonify(send_de_res)


@bp.route("/PearsonCorr/Normal", methods=["GET"])
def get_pearson_corr_normal():
    return jsonify(adata.uns["vision_lca"])


@bp.route("/PearsonCorr/Meta", methods=["GET"])
def get_pearson_corr_meta():
    empty = {"sig_labels": [], "proj_labels": [], "zscores": [], "pvals": []}
    return jsonify(adata.uns.get("vision_lca_meta", empty))


@bp.route("/Proteins/list", methods=["GET"])
def get_protein_list():
    autocorr = adata.uns.get("vision_protein_autocorr")
    if autocorr is None:
        return jsonify([])
    result = []
    for name, row in autocorr.iterrows():
        result.append({
            "name": name,
            "autocorrelation": row["c_prime"],
            "pval": row["pvals"],
            "fdr": row["fdr"],
        })
    return jsonify(result)


@bp.route("/Proteins/Autocorrelation", methods=["GET"])
def get_protein_autocorrelation():
    autocorr = adata.uns.get("vision_protein_autocorr")
    if autocorr is None:
        return jsonify({})
    return jsonify(autocorr.to_dict(orient="index"))


@bp.route("/Proteins/<prot_name>/Scores/<cluster_var>", methods=["GET"])
def get_protein_scores(prot_name, cluster_var):
    autocorr = adata.uns.get("vision_protein_autocorr")
    if autocorr is None or prot_name not in autocorr.index:
        return jsonify({}), 404

    prot_key = adata.uns.get("vision_protein_differential_key")
    if prot_key is None:
        return jsonify({}), 404

    mat = adata.obsm[prot_key]
    if isinstance(mat, pd.DataFrame):
        if prot_name not in mat.columns:
            return jsonify({}), 404
        values = mat[prot_name].values.tolist()
    else:
        prot_names = autocorr.index.tolist()
        if prot_name not in prot_names:
            return jsonify({}), 404
        col_idx = prot_names.index(prot_name)
        values = np.asarray(mat)[:, col_idx].tolist()

    return jsonify(dict(cells=adata.obs_names.tolist(), values=values))


@bp.route("/Proteins/<prot_name>/Meta/<cluster_var>", methods=["GET"])
def get_protein_meta(prot_name, cluster_var):
    prot_adata = data_accessor.protein_adata
    if prot_adata is None or cluster_var not in prot_adata.obs.columns:
        return jsonify({}), 404

    key = f"rank_genes_groups_{cluster_var}"
    if key not in prot_adata.uns:
        return jsonify({}), 404

    proj_labels = list(prot_adata.obs[cluster_var].astype("category").cat.categories)
    stats = []
    pvals = []
    for level in proj_labels:
        temp_df = sc.get.rank_genes_groups_df(prot_adata, key=key, group=str(level))
        temp_df.set_index("names", inplace=True)
        if prot_name in temp_df.index:
            stats.append(float(temp_df.loc[prot_name, "scores"]))
            pvals.append(float(temp_df.loc[prot_name, "pvals_adj"]))
        else:
            stats.append(0.0)
            pvals.append(1.0)

    return jsonify(dict(proj_labels=proj_labels, stats=stats, pvals=pvals))


@bp.route("/Download/Selections", methods=["GET"])
def download_selections():
    selections = {k: data_accessor.get_cells_selection(k) for k in data_accessor.cells_selections}
    return jsonify(selections)


@bp.route("/Download/DE", methods=["GET"])
def download_de():
    return jsonify([])
