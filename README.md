<img src="https://raw.githubusercontent.com/YosefLab/VISION/master/man/figures/logo.svg" align="right" width="150" />

# Functional Interpretation for <br/> scRNA-seq Data

[![Tests][badge-tests]][link-tests]
[![Documentation][badge-docs]][link-docs]

[badge-tests]: https://img.shields.io/github/actions/workflow/status/yoseflab/visionpy/test.yaml?branch=main
[link-tests]: https://github.com/yoseflab/visionpy/actions/workflows/test.yml
[badge-docs]: https://img.shields.io/readthedocs/visionpy

NOTE: THIS PACKAGE IS UNDER ACTIVE DEVELOPMENT AND IS NOT YET READY FOR USE.

This is a Python port of the VISION [R package](https://github.com/yoseflab/vision). VISION aids in the interpretation of single-cell RNA-seq (scRNA-seq) data by selecting for gene signatures which describe coordinated variation between cells. While the software only requires an expression matrix and a signature library (available in online databases), it is also designed to integrate into existing scRNA-seq analysis pipelines by taking advantage of precomputed dimensionality reductions, trajectory inferences or clustering results. The results of this analysis are made available through a dynamic web-app which can be shared with collaborators without requiring them to install any additional software.

-   [Nature Communications publication](https://www.nature.com/articles/s41467-019-12235-0)
-   [Full R Documentation](https://yoseflab.github.io/VISION/)
-   [R package](https://github.com/yoseflab/vision)

## Installing visionpy

You need to have Python 3.8 or newer installed on your system. If you don't have
Python installed, we recommend installing [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

There are several alternative options to install visionpy:

<!--
1) Install the latest release of `visionpy` from `PyPI <https://pypi.org/project/visionpy/>`_:

```bash
pip install visionpy
```
-->

1. Install the latest release on PyPI:

```bash
pip install visionpy-sc
```

2. Install the latest development version:

```bash
pip install git+https://github.com/yoseflab/visionpy.git@main
```

## How to run visionpy

### From the command line

```
visionpy --adata ./my_adata.h5ad --norm_data_key use_raw --compute_neighbors_on_key X_scvi --name Test Vision
```

### From Python

```python
from visionpy.api import start_vision
from visionpy import signatures_from_gmt

adata.varm["signatures"] = signatures_from_gmt(["./signatures.gmt"], adata)
start_vision(
    adata=adata,
    name="Test Session",
    norm_data_key="log1pcp10k",
    compute_neighbors_on_key="X_pca",
    signture_varm_key="signatures",
    name="Test Vision",
)
```

[link-docs]: https://visionpy.readthedocs.io
