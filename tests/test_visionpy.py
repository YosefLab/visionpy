import numpy as np
import pandas as pd
import scanpy

from visionpy.api import _prepare_vision


def test_preparation(save_path):
    scanpy.settings.datasetdir = save_path
    adata = scanpy.datasets.pbmc3k_processed()
    N_SIGS = 10
    adata.raw.varm["signatures"] = pd.DataFrame(index=adata.raw.var_names, data=np.zeros((adata.raw.n_vars, N_SIGS)))
    adata.raw.varm["signatures"].iloc[:20] = 1
    adata.raw.varm["signatures"].iloc[-25:] = -1
    adata.uns["sig_names"] = [f"Sig_{i}" for i in range(N_SIGS)]

    _prepare_vision(
        adata,
        name="Test session",
        norm_data_key="use_raw",
        signature_varm_key="signatures",
        signature_names_uns_key="sig_names",
    )
