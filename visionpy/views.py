from collections import OrderedDict

from flask import Blueprint, jsonify, render_template
from visionpy import data_accessor

bp = Blueprint("api", __name__)

adata = data_accessor.adata


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
    proj_dict["Obs metadata"] = adata.obs._get_numeric_data().columns.tolist()

    if "X_umap" in adata.obsm_keys():
        proj_dict.move_to_end("X_umap", last=False)

    return jsonify(proj_dict)


@bp.route("/Projections/<proj_name>/coordinates/<proj_col>", methods=["GET"])
def get_projection_column(proj_name, proj_col):
    if proj_name == "Obs metadata":
        data = adata.obs[proj_col].values.tolist()
    else:
        column_ind = int(proj_col[-1] - 1)
        data = adata.obsm[proj_name][:, column_ind].values.tolist()

    return jsonify(data)


@bp.route("/Tree/Projections/list", methods=["GET"])
def get_tree_projection_list():
    # TODO: Implement
    return jsonify([])


@bp.route("/SessionInfo", methods=["GET"])
def get_session_info():
    info = {}
    info["name"] = "name"
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
    return jsonify([])


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
        dict(cells=adata.obs_names.tolist(), values=adata[:, gene_name].X.tolist())
    )


@bp.route("/FilterGroup/SigClusters/Meta", methods=["GET"])
def get_sigclusters_meta():
    return jsonify(["1"] * len(adata.obs.columns))


@bp.route("/Clusters/<cluster_variable>/SigProjMatrix/Meta", methods=["GET"])
def get_sigprojmatrix_meta(cluster_variable):
    return jsonify([])
