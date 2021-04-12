from visionpy.api import start_vision
import scanpy

adata = scanpy.datasets.pbmc3k_processed()

start_vision(adata, debug=True, norm_data_key="use_raw")
