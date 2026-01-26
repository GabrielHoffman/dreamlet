# Extract cell counts

Extract matrix of cell counts from `SingleCellExperiment`

## Usage

``` r
cellCounts(x)
```

## Arguments

- x:

  a `SingleCellExperiment`

## Value

matrix of cell counts with samples as rows and cell types as columns

## See also

[`computeCellCounts()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/computeCellCounts.md)

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

# get matrix of cell counts for each sample
cellCounts(pb)
#>         B cells CD14+ Monocytes CD4 T cells CD8 T cells FCGR3A+ Monocytes
#> ctrl101     100             100         100          74                80
#> ctrl107      44             100         100          20                32
#> stim101     100             100         100         100               100
#> stim107      54             100         100          15                37
```
