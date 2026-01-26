# Compute log normalized counts

Compute normalized counts as log2 counts per million

## Usage

``` r
computeLogCPM(
  sce,
  lib.size = colSums2(counts(sce)),
  prior.count = 2,
  scaledByLib = FALSE
)
```

## Arguments

- sce:

  `SingleCellExperiment` with counts stored as `counts(sce)`

- lib.size:

  library size for each cell

- prior.count:

  average count to be added to each observation to avoid taking log of
  zero

- scaledByLib:

  if `TRUE`, scale pseudocount by `lib.size`. Else do standard constant
  pseudocount addition

## Value

matrix of log CPM values

## Details

This function gives same result as `edgeR::cpm(counts(sce), log=TRUE)`

## See also

also [`edgeR::cpm()`](https://rdrr.io/pkg/edgeR/man/cpm.html)

## Examples

``` r
library(muscat)
library(SingleCellExperiment)

data(example_sce)

logcounts(example_sce) <- computeLogCPM(example_sce)
```
