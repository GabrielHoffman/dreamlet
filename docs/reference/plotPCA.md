# Plot PCA of gene expression for an assay

Compute PCA of gene expression for an assay, and plot samples coloring
by outlier score

## Usage

``` r
# S4 method for class 'list'
plotPCA(
  object,
  assays = names(object),
  nPC = 2,
  robust = FALSE,
  ...,
  maxOutlierZ = 20,
  nrow = 2,
  size = 2,
  fdr.cutoff = 0.05
)
```

## Arguments

- object:

  `dreamletProcessedData` from
  [`processAssays()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/processAssays.md)
  or a `list` from
  [`residuals()`](https://rdrr.io/r/stats/residuals.html)

- assays:

  assays / cell types to analyze

- nPC:

  number of PCs to uses for outlier score with
  [`outlier()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/outlier.md)

- robust:

  use robust covariance method, defaults to `FALSE`

- ...:

  arguments passed to
  [`MASS::cov.rob()`](https://rdrr.io/pkg/MASS/man/cov.rob.html)

- maxOutlierZ:

  cap outlier z-scores at this value for plotting to maintain consistent
  color scale

- nrow:

  number of rows in plot

- size:

  size passed to `geom_point()`

- fdr.cutoff:

  FDR cutoff to determine outlier

## See also

[`outlierByAssay()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/outlierByAssay.md)

## Examples

``` r
library(muscat)
library(SingleCellExperiment)

data(example_sce)

# create pseudobulk for each sample and cell cluster
pb <- aggregateToPseudoBulk(example_sce,
  assay = "counts",
  cluster_id = "cluster_id",
  sample_id = "sample_id",
  verbose = FALSE
)

# voom-style normalization
res.proc <- processAssays(pb, ~group_id)
#>   B cells...
#> 0.055 secs
#>   CD14+ Monocytes...
#> 0.086 secs
#>   CD4 T cells...
#> 0.062 secs
#>   CD8 T cells...
#> 0.046 secs
#>   FCGR3A+ Monocytes...
#> 0.081 secs

# PCA to identify outliers
# from normalized expression
plotPCA( res.proc, c("B cells", "CD14+ Monocytes"))


# Run on regression residuals
#-----------------------------

# Regression analysis
fit = dreamlet(res.proc, ~ group_id)
#>   B cells...
#> 0.05 secs
#>   CD14+ Monocytes...
#> 0.064 secs
#>   CD4 T cells...
#> 0.061 secs
#>   CD8 T cells...
#> 0.033 secs
#>   FCGR3A+ Monocytes...
#> 0.061 secs

# Extract regression residuals
residsObj = residuals(fit)

# PCA on residuals
plotPCA( residsObj, c("B cells", "CD14+ Monocytes"))
```
