from flask import (
    Blueprint,
)

bp = Blueprint("hello", __name__)


@bp.route("/hello")
def hello():
    return "Hello, World!"
