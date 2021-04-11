from visionpy import create_app, data_accessor

import scanpy

adata = scanpy.datasets.pbmc3k_processed()


def run():
    data_accessor.adata = adata
    app = create_app()
    app.run(threaded=False, processes=1)
