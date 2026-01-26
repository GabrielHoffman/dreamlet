# Compute normalized counts

Compute normalized counts as counts per million

## Usage

``` r
computeNormCounts(sce)
```

## Arguments

- sce:

  `SingleCellExperiment` with counts stored as `counts(sce)`

## Value

matrix of CPM values

## Details

This function gives same result as `edgeR::cpm(counts(sce), log=FALSE)`

## See also

also [`edgeR::cpm()`](https://rdrr.io/pkg/edgeR/man/cpm.html)

## Examples

``` r
library(muscat)
library(SingleCellExperiment)

data(example_sce)

normcounts(example_sce) <- computeNormCounts(example_sce)
```
