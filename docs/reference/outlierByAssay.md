# Outlier analysis for each assay

Compute outlier score for each sample in each assay using
[`outlier()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/outlier.md)
run on the top principal components. Mahalanobis distance is used for
outlier detect and multivariate normal assumption is used to compute
p-values

## Usage

``` r
outlierByAssay(object, assays = names(object), nPC = 2, robust = FALSE, ...)
```

## Arguments

- object:

  `dreamletProcessedData` from
  [`processAssays()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/processAssays.md)

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

## Value

- `ID`::

  sample identifier

- `assay`::

  specify assay

- `PCs`::

  principal components

- `chisq`::

  mahalanobis distance that is distributed as chisq(k) k = nPC if the
  data is multivariate gaussian

- `z`::

  z-score corresponding to the chisq distance

## See also

[`outlier()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/outlier.md)

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
#> 0.054 secs
#>   CD14+ Monocytes...
#> 0.087 secs
#>   CD4 T cells...
#> 0.064 secs
#>   CD8 T cells...
#> 0.039 secs
#>   FCGR3A+ Monocytes...
#> 0.1 secs

# Compute PCs and outlier scores
outlierByAssay( res.proc, c("B cells", "CD14+ Monocytes")) 
#> # A tibble: 8 × 7
#>   ID      assay              PC1    PC2 chisq     z pValue
#>   <chr>   <chr>            <dbl>  <dbl> <dbl> <dbl>  <dbl>
#> 1 ctrl101 B cells         -0.457  1.09   1.40 0.681  0.496
#> 2 ctrl107 B cells         -1.19  -0.626  1.81 0.834  0.404
#> 3 stim101 B cells          0.909  0.575  1.16 0.582  0.561
#> 4 stim107 B cells          0.738 -1.04   1.63 0.768  0.443
#> 5 ctrl101 CD14+ Monocytes -0.976  0.709  1.45 0.701  0.483
#> 6 ctrl107 CD14+ Monocytes -0.736 -0.971  1.48 0.713  0.476
#> 7 stim101 CD14+ Monocytes  0.715  1.00   1.52 0.725  0.469
#> 8 stim107 CD14+ Monocytes  0.998 -0.741  1.54 0.736  0.462
```
