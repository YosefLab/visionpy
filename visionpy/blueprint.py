from collections import OrderedDict
import json
import scanpy as sc
import numpy as np
import pandas as pd
from flask import Blueprint, jsonify, render_template, request
from visionpy import data_accessor

bp = Blueprint("api", __name__)

adata = data_accessor.adata


class ServerSigProjMatrix:
    def __init__(self, zscores, pvals, proj_labels, sig_labels) -> None:
        self.zscores = zscores
        self.pvals = pvals
        self.proj_labels = proj_labels
        self.sig_labels = sig_labels

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
    for k in adata.obsm_keys():
        if k[:2] == "X_":
            name = k.split("_")[1]
            proj_dict[k] = [
                name.upper() + "{}".format(i + 1) for i in range(adata.obsm[k].shape[1])
            ]
    # get only numeric columns
    proj_dict["Obs_metadata"] = data_accessor.numeric_obs_cols

    if "X_umap" in adata.obsm_keys():
        proj_dict.move_to_end("X_umap", last=False)

    return jsonify(proj_dict)


@bp.route("/Projections/<proj_name>/coordinates/<proj_col>", methods=["GET"])
def get_projection_column(proj_name, proj_col):
    if proj_name == "Obs_metadata":
        data = adata.obs[proj_col].values.tolist()
    else:
        column_ind = int(proj_col[-1]) - 1
        data = adata.obsm[proj_name]
        if isinstance(data, pd.DataFrame):
            data = data.iloc[:, column_ind].values.tolist()
        else:
            data = data[:, column_ind].tolist()

    # the order matters here
    data = np.vstack([data, adata.obs_names.tolist()]).transpose()

    # need a list of lists, where each internal list is a coord value
    # and then a barcode value
    return jsonify(data.tolist())


@bp.route("/Tree/Projections/list", methods=["GET"])
def get_tree_projection_list():
    # TODO: Implement
    return jsonify([])


@bp.route("/SessionInfo", methods=["GET"])
def get_session_info():
    info = {}
    info["name"] = adata.uns["vision_session_name"]
    info["has_tree"] = False
    info["meta_sigs"] = adata.obs.columns.tolist()
    # info[["pooled"]] <- object@params$micropooling$pool
    info["pooled"] = False
    info["ncells"] = adata.n_obs
    # info[["has_sigs"]] <- length(object@sigData) > 0
    info["has_sigs"] = "vision_signatures" in adata.obsm_keys()
    # TODO: Fix
    info["has_proteins"] = False
    # TODO: Fix
    info["has_lca"] = False
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
    return "No clusters"


@bp.route("/Cells/Selections", methods=["GET"])
def get_cells_selections():

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


@bp.route("/Signature/Meta/<sig_name>", methods=["GET"])
def get_signature_meta(sig_name):
    return jsonify(
        dict(cells=adata.obs_names.tolist(), values=adata.obs.loc[:, sig_name].tolist())
    )


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
    cols = adata.obsm["vision_signatures"].columns.to_list()
    clusters = {}
    for c in cols:
        clusters[c] = 1
    return jsonify(clusters)


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
    proj_labels = ["Score"] + list(
        adata.obs[cluster_variable].astype("category").cat.categories
    )

    scores = adata.uns["vision_obs_df_scores"]["c_prime"].to_numpy().reshape(-1, 1)
    sigs_by_projs_stats = pd.DataFrame(
        index=sig_labels, columns=proj_labels[1:], data=0
    )
    sigs_by_projs_pvals = pd.DataFrame(
        index=sig_labels, columns=proj_labels[1:], data=0
    )
    obs_adata = data_accessor.obs_adata
    # TODO: test for categorical data with chi sq
    for p in proj_labels[1:]:
        temp_df = sc.get.rank_genes_groups_df(
            obs_adata, key="rank_genes_groups_{}".format(cluster_variable), group=p
        )
        temp_df.set_index("names", inplace=True)
        sigs_by_projs_stats.loc[temp_df.index, p] = temp_df["scores"].copy()
        sigs_by_projs_pvals.loc[temp_df.index, p] = temp_df["pvals_adj"].copy()
        for cat_c in data_accessor.cat_obs_cols:
            key = "chi_sq_{}_{}".format(cluster_variable, p)
            sigs_by_projs_stats.loc[cat_c, p] = obs_adata.uns[key]["stat"]
            sigs_by_projs_pvals.loc[cat_c, p] = obs_adata.uns[key]["pval"]

    stats = np.hstack([scores, sigs_by_projs_stats.to_numpy()]).tolist()
    pvals = np.hstack([np.zeros_like(scores), sigs_by_projs_pvals.to_numpy()]).tolist()

    matrix = ServerSigProjMatrix(
        sig_labels=sig_labels, proj_labels=proj_labels, zscores=stats, pvals=pvals
    )
    return jsonify(matrix.prepare_json())


@bp.route("/Clusters/<cluster_variable>/SigProjMatrix/Normal", methods=["GET"])
def get_sigprojmatrix_normal(cluster_variable):

    sig_labels = adata.obsm["vision_signatures"].columns.tolist()
    # TODO: amortize computation in metalevels route
    proj_labels = ["Score"] + list(
        adata.obs[cluster_variable].astype("category").cat.categories
    )

    # TODO: properly compute all information
    scores = adata.uns["vision_signature_scores"]["c_prime"].to_numpy().reshape(-1, 1)

    sigs_by_projs_stats = pd.DataFrame(index=sig_labels, columns=proj_labels[1:])
    sigs_by_projs_pvals = pd.DataFrame(index=sig_labels, columns=proj_labels[1:])
    sig_adata = data_accessor.sig_adata
    for p in proj_labels[1:]:
        temp_df = sc.get.rank_genes_groups_df(
            sig_adata, key="rank_genes_groups_{}".format(cluster_variable), group=p
        )
        temp_df.set_index("names", inplace=True)
        temp_df = temp_df.loc[sigs_by_projs_stats.index]
        sigs_by_projs_stats[p] = temp_df["scores"].copy()
        sigs_by_projs_pvals[p] = temp_df["pvals_adj"].copy()

    stats = np.hstack([scores, sigs_by_projs_stats.to_numpy()]).tolist()
    pvals = np.hstack([np.zeros_like(scores), sigs_by_projs_pvals.to_numpy()]).tolist()

    matrix = ServerSigProjMatrix(
        sig_labels=sig_labels, proj_labels=proj_labels, zscores=stats, pvals=pvals
    )
    return jsonify(matrix.prepare_json())


@bp.route("/Cell/<cell_id>/Meta", methods=["GET", "POST"])
def get_cell_metadata(cell_id):
    if request.method == "GET":
        cell = adata.obs.loc[cell_id].to_dict()
        for k, v in cell.items():
            cell[k] = str(v)
        return jsonify(cell)
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
        categorical_stats[c] = dist

    return_dict = dict(numeric=numerical_stats, factor=categorical_stats)

    return jsonify(return_dict)
