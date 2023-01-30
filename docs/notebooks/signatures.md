# Gene signatures

## What is a Signature?

A gene signature is a set of genes involved in some biological process.

Signatures come in two flavors:

-   **Unsigned** - A set of genes that have some common annotation. For example, genes involved in a pathway of interest. These signatures commonly derive from manual annotation approaches These signatures commonly derive from manual annotation approaches.
-   **Signed** - A set of genes which describes the _contrast_ between two conditions. Signed signatures compare condition A vs. condition B and have a set of 'positive' genes which are increased in A (compared to B) and a set of 'negative' genes which are decreased in A (compared to B). Signed signatures are often defined computationally as the result of a differential expression test in some previous experiment.

## Where to find Signatures

A great resource for gene signatures is **MSigDB**, curated by the Broad institute. Signatures can be browsed, searched and downloaded from [here](http://software.broadinstitute.org/gsea/msigdb/) as .gmt files, then provided to VISION to be included in the analysis.

Alternately, many signature libraries in .gmt format are available through **[Enrichr](http://amp.pharm.mssm.edu/Enrichr)** on their ["Libraries"](http://amp.pharm.mssm.edu/Enrichr/#stats) page.

Finally, if you have your own lists of genes, you can define you own signatures following the instructions below.

## Creating Signatures Manually

If there is a set of proprietary genes of interest, a user-defined signature can be created in two ways:

### 1. Creating a Signature in Python

Once a set of genes that are up or down regulated in the process or cell type of interest are selected, creating a Signature object from them is relatively straightforward:

For the sigData vector, the names represent gene names and the values (1 or -1) represent the 'sign' of the gene. For an unsigned signature, just use 1 for the value of every gene.

A list of these sig objects can then be passed into the VISION object constructor instead of paths to .gmt files.

Additionally, you can mix and match user-created signatures and paths to signature files in the `signatures` argument. When doing this, signatures in each library file are read and combined with the user-created signatures.

### 2. Create your own .gmt files

Signature files are supported in the .gmt format - a textual format which is easy to create and view

The file is tab-delimited and contains one signature per line, (however 'positive' and 'negative' genes are split into two lines)

Each line should look like this:

```
<Signature Name> TAB <Signature Description> TAB <Gene1> TAB <Gene2> … (etc)
```

To denote signed signatures, use two lines to show the signature, with the "positive" genes in one line and the "negative" genes in the other. Add "\_plus" to the signature name on the line with the positive genes and "\_minus" to the signature name on the line with the negative genes.

For example:

```
MEMORY_VS_NAIVE_CD8_TCELL_plus    GSE16522   RHOC    OFD1     MLF1   ...
MEMORY_VS_NAIVE_CD8_TCELL_minus   GSE16522   PTPRK   S100A5   IL1A   ...
BCELL_VS_LUPUS_BCELL_plus         GSE10325   OXT     KCNH2    BTBD7  ...
BCELL_VS_LUPUS_BCELL_minus        GSE10325   VAMP5   WSB2     CCR2   ...
```

## Signature Scores

For each signature an overall measure of expression is evaluated for each cell as the _Signature Score_. This is computed as:

$$ s*j = \frac{1}{|G*{pos}| + |G*{neg}} \Bigg(\sum*{g \in G*{pos}}{X*{gj}} - \sum*{g \in G*{neg}}{X\_{gj}}\Bigg) $$

Where:

-   $s_j$ is the signature score in cell _j_
-   $X_{gj}$ is the _log_, _scaled_ expression of gene _g_ in cell _j_
-   $G_{pos}$ and $G_{neg}$ are the set of 'positive' and 'negative' signature genes respectively

VISION expects the input expression data to already be scaled (i.e., the gene counts in each cell are divided by the total number of counts in the cell and multiplied by some constant such as 10,000 or the median number of UMI across cells) and normalized (significant technical confounders regressed-out). Internally, VISION additionally log-transforms ($log_2(x+1)$) the expression data to compress the dynamic range. This is done so that genes which are naturally expressed at higher magnitudes don't dominate the signature score as heavily.

However, with just the above formulation, we have observed that depending on the scaling or normalization method, signature scores can still be highly correlated with cell-level metrics such as the number of UMI per cell. To account for this, we note that the expected value of a randomly-drawn signature in cell _j_ is:

$$ \frac{|G*{pos}| - |G*{neg}|}{|G*{pos}| + |G*{neg}|} \bar{X_j} $$

And the expected variance of a random signature score is:

$$ \frac{var(X*j)}{|G*{pos}| + |G\_{neg}|} $$

Therefore we can remove global cell-specific distributional effects from the signature scores by first Z-normalizing the expression data, $X_j$, within each cell (so that $\bar{X_j}$ is 0 and $var(X_j)$ is 1) prior to calculating signature scores as defined above. This is the default operation performed by VISION (controlled through the `sig_norm_method` parameter).
