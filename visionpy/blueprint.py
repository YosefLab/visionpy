from collections import OrderedDict

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
    proj_dict["Obs_metadata"] = adata.obs._get_numeric_data().columns.tolist()

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
    info["has_sigs"] = False
    # TODO: Fix
    info["has_proteins"] = False
    # TODO: Fix
    info["has_lca"] = False
    return jsonify(info)


@bp.route("/Clusters/MetaLevels", methods=["GET"])
def get_clusters_metalevels():
    num_cols = adata.obs._get_numeric_data().columns.tolist()
    cols = adata.obs.columns.tolist()
    cat_vars = list(set(cols) - set(num_cols))

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
    # TODO
    return jsonify(
        dict(cells=adata.obs_names.tolist(), values=np.zeros((adata.n_obs).tolist()))
    )


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

    # TODO: properly compute all information
    scores = adata.uns["vision_obs_df_scores"]["c_prime"].to_numpy().reshape(-1, 1)
    zscores = np.hstack(
        [scores, np.zeros((len(sig_labels), len(proj_labels) - 1))]
    ).tolist()
    pvals = np.zeros_like(zscores).tolist()

    matrix = ServerSigProjMatrix(
        sig_labels=sig_labels, proj_labels=proj_labels, zscores=zscores, pvals=pvals
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
        # TODO IMPLEMENT
        # get subset from json postBody
        # adata.obs == METADATA
        # adata.obs._get_numeric_data() is a dataframe subset with numeric columns
        # need dict[column] = dict with keys ("Min", "Median", "Max") for numeric columns
        # for categorical obs
        # need dict[column] = dict with keys (categories of column) value is frequency (e.g., 90 for 90%)
        # take two dictionaries and put in dictionary with key "numeric" for first, "factor" for second
        pass
