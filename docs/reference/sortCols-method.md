# Sort variance partition statistics

Sort variance partition statistics

## Usage

``` r
# S4 method for class 'vpDF'
sortCols(
  x,
  FUN = sum,
  decreasing = TRUE,
  last = c("Residuals", "Measurement.error"),
  ...
)
```

## Arguments

- x:

  object returned by
  [`fitVarPart()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/fitVarPart.md)

- FUN:

  function giving summary statistic to sort by. Defaults to sum

- decreasing:

  logical. Should the sorting be increasing or decreasing?

- last:

  columns to be placed on the right, regardless of values in these
  columns

- ...:

  other arguments to sort

## Value

`data.frame` with columns sorted by mean value, with Residuals in last
column

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
#> 0.066 secs
#>   CD14+ Monocytes...
#> 0.087 secs
#>   CD4 T cells...
#> 0.063 secs
#>   CD8 T cells...
#> 0.039 secs
#>   FCGR3A+ Monocytes...
#> 0.11 secs

# variance partitioning analysis
vp <- fitVarPart(res.proc, ~group_id)
#>   B cells...
#> 0.61 secs
#>   CD14+ Monocytes...
#> 0.83 secs
#>   CD4 T cells...
#> 0.64 secs
#>   CD8 T cells...
#> 0.39 secs
#>   FCGR3A+ Monocytes...
#> 0.77 secs
#> 

# Summarize variance fractions genome-wide for each cell type
plotVarPart(sortCols(vp))

```
