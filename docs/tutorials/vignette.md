# Introduction to visionpy

## Preliminaries

If you have yet to install visionpy, we recommend installing the package from Github to install this package. Full source code can be found at the visionpy Github repository, available [here](http://www.github.com/YosefLab/visionpypy).

If you encounter errors in the installation, it is likely because a dependency is not installed correctly. Pay special attention to the error message, and in particular whether or not it notifies you that specific dependencies are not found.

After a successful installation, you can proceed by loading visionpy using `import visionpy`.

## Using visionpy

Running an analysis with visionpy consists of three steps:

1. Creating the visionpy object
2. Running the `analyze` function
3. Browsing results

### Creating the visionpy object

Creating the visionpy object requires a gene expression matrix and _either_ a list of Gene Signatures or a data.frame of meta-data.

In this example - both are provided.

**Expression Data**

The provided expression data should be scaled and normalized. The example above shows just a simple UMI-scaling, but it is recommended to apply more advanced normalization procedures such as batch correction or removal of technical confounders.

The expression data should not be log-transformed prior to loading into visionpy.

**Signatures**

Signatures can be provided as a list of paths to signature files (\*.gmt) or Signature objects.

See [the signature vignette](Signatures.html) for more information on finding or creating gene Signatures.

**Meta Data**

An R data.frame with cell-level meta-data. This could be confounding variables (e.g. percentage of mitochondrial RNA, number of genes detected) or experimental covariates (e.g. genotype, donor, batch).

This input is optional if Signatures are provided.

**Other Options**

Other options and inputs can be provided to customize how visionpy runs. For information on this, see the "Customizing visionpy Analysis" section below.

### Running an Analysis

To run an analysis, simply call the analyze function:

### Viewing Results

With the processed visionpy object, a dynamic web report can be generated with the `viewResults()` function.

This will launch a browser running the interactive report.

Other options (port, host, browser) can be provided to control how this occurs. For example, if you are launching a report on a remote server (such as an AWS instance) and want to make it accessible to others, run this with `host="0.0.0.0"`, some selected port number (e.g. `port=8888`), and `browser=FALSE` (so a browser isn't auto-opened). Then the report should be available at "\<your instance IP address\>:8888". (Note: You will also likely need to enable inbound traffic on your selected port for this to work correctly).

Alternately, you can work with the visionpy object directly in R. For example:

### Customizing the Latent Space

visionpy requires a latent space to model the similarity between cells (used to determine a cell's local neighborhood).

### Adding 2d Projections

Two-dimensional projections are used to visualize the data in the output report.

By default, visionpy computes tSNE on the latent space for visualization. However, other options are availabled via the `projection_methods` argument.

Often times, users will have pre-computed projections that they would like to use for visualizing their data (e.g. a pre-computed tSNE, or UMAP projection). In this case, the `addProjection()` method can be used to add this view of the data to the output report.

```{r, collapse=T, message=F, results=F, eval=F}
projection <- read.csv("umap_results.csv")

# projection is a matirx or data.frame of dimension (Cells x 2)

vis <- addProjection(vis, "UMAP", projection)
```

### Hotspot analysis

As of version 3.0.0, we have enabled users to perform de-novo gene module identification [Hotspot](<https://www.cell.com/cell-systems/fulltext/S2405-4712(21)00114-9>) (DeTomaso and Yosef, _Cell Systems_ 2021) from within visionpy.

As described in the original Hotspot vignette, we'll use the `danb` model for standard single-cell RNA-seq datasets. While this works in general, pay special attention to characteristics of your data (e.g., in spatial examples where capture is low, we recommend using the `bernoulli` model).

For more information about the analysis pipeline & Hotspot API, you can refer the documentation website [here](https://yoseflab.github.io/Hotspot/index.html) and [the Phylovisionpy vignette](phylovisionpy.html).
