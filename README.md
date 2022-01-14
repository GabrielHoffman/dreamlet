
<img src="man/figures/logo.png" align="right" alt="" width="160" />

### Scalable differential expression analysis of single cell transcriptomics datasets with complex study designs

The `dreamlet` package enables differential expression analysis on multi-sample single cell datasets using linear (mixed) models with precision weights.

Major functionality of `dreamlet` package using the [Bioconductor](https://www.bioconductor.org) [`SingleCellExperiment`](https://www.bioconductor.org/packages/SingleCellExperiment/) interface:

+ [`aggregateToPseudoBulk()`](https://gabrielhoffman.github.io/dreamlet/reference/aggregateToPseudoBulk.html)      Computationally scalable evaluation of pseudobulk from raw counts
+ [`processAssays()`](https://gabrielhoffman.github.io/dreamlet/reference/processAssays.html)                      Normalize counts, compute precision weights
+ [`fitVarPart()`](https://gabrielhoffman.github.io/dreamlet/reference/fitVarPart.html)                            Variance partitioning analysis
+ [`dreamlet()`](https://gabrielhoffman.github.io/dreamlet/reference/dreamlet.html)                                Differential expression analysis across samples accounting for study design
+ [`dreamletCompareClusters()`](https://gabrielhoffman.github.io/dreamlet/reference/dreamletCompareClusters.html)  Differential expression analysis across cell clusters  for study design
+ [`zenith_gsa()`](https://gabrielhoffman.github.io/dreamlet/reference/zenith_gsa.html)                            Gene set analysis using full spectrum of tests statistics


### Install
```r
# repo is currently private, so need to include your userid and password
devtools::install_github("GabrielHoffman/dreamlet", auth_token=XXXXX)
```

#### Dependencies
In case code above doesn't install these automatically
```r
devtools::install_github("GabrielHoffman/variancePartition")
devtools::install_github("GabrielHoffman/zenith")
devtools::install_github("stephenslab/mashr")
```


