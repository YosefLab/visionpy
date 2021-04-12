from flask import Blueprint, jsonify, render_template
from visionpy import data_accessor

bp = Blueprint("api", __name__)

adata = data_accessor.adata


@bp.route("/")
def app():
    return render_template("Results.html")


@bp.route("/Projections/list", methods=["GET"])
def get_projection_list():
    return jsonify(["X_umap", "X_pca"])


# @bp.route("/Projections/<proj_name>/coordinates/<proj_col>", methods=("GET"))
# def get_projection_coords(sig_name, cluster_var):
#     pass


# @bp.route("/Expression/Gene/<gene_name>", methods=("GET"))
# def get_gene_expression(sig_name, cluster_var):

#     if request.method == "GET":
#         return None
