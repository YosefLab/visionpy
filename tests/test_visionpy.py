from visionpy.api import start_vision
import scanpy
import pandas as pd
import numpy as np

adata = scanpy.datasets.pbmc3k_processed()
N_SIGS = 10
adata.raw.varm["signatures"] = pd.DataFrame(
    index=adata.raw.var_names, data=np.zeros((adata.raw.n_vars, N_SIGS))
)
adata.raw.varm["signatures"].iloc[:20] = 1
adata.raw.varm["signatures"].iloc[-25:] = -1
adata.uns["sig_names"] = ["Sig_{}".format(i) for i in range(N_SIGS)]

start_vision(
    adata,
    name="Test session",
    debug=True,
    norm_data_key="use_raw",
    signature_varm_key="signatures",
    signature_names_uns_key="sig_names",
    use_raw_for_signatures=True,
)
