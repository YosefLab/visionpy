"""visionpy."""

import logging
import os
from importlib.metadata import version

from flask import Flask
from flask_compress import Compress
from rich.console import Console
from rich.logging import RichHandler

from .anndata import AnnDataAccessor
from .signature import compute_signatures_anndata, signatures_from_gmt

__version__ = version("visionpy-sc")

logger = logging.getLogger(__name__)
# set the logging level
logger.setLevel(logging.INFO)

# nice logging outputs
console = Console(force_terminal=True)
if console.is_jupyter is True:
    console.is_jupyter = False
ch = RichHandler(show_path=False, console=console, show_time=False)
formatter = logging.Formatter("visionpy: %(message)s")
ch.setFormatter(formatter)
logger.addHandler(ch)

# this prevents double outputs
logger.propagate = False


data_accessor = AnnDataAccessor()
compress = Compress()


# https://flask.palletsprojects.com/en/1.1.x/tutorial/factory/
def create_app(test_config=None):
    # create and configure the app
    app = Flask(__name__, instance_relative_config=True, static_url_path="")
    app.config.from_mapping(
        SECRET_KEY="dev",
        DATABASE=os.path.join(app.instance_path, "flaskr.sqlite"),
    )

    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_pyfile("config.py", silent=True)
    else:
        # load the test config if passed in
        app.config.from_mapping(test_config)

    # ensure the instance folder exists
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    from . import blueprint

    app.register_blueprint(blueprint.bp)
    # https://stackoverflow.com/questions/43263356/prevent-flask-jsonify-from-sorting-the-data/43263483
    # app.config["JSON_SORT_KEYS"] = False
    compress.init_app(app)

    return app
