# Get cell counts with metadata

Get cell counts with metadata for each sample

## Usage

``` r
computeCellCounts(sce, annotation, sampleIDs)
```

## Arguments

- sce:

  `SingleCellExperiment`

- annotation:

  string indicating column in `colData(sce)` storing cell type
  annotations

- sampleIDs:

  string indicating column in `colData(sce)` storing sample identifers

## Value

`matrix` storing cell counts

## Examples

``` r
library(muscat)
library(SingleCellExperiment)

data(example_sce)

counts <- computeCellCounts(example_sce, "cluster_id", "sample_id")

counts[1:4, 1:4]
#>         B cells CD14+ Monocytes CD4 T cells CD8 T cells
#> ctrl101     100             100         100          74
#> ctrl107      44             100         100          20
#> stim101     100             100         100         100
#> stim107      54             100         100          15
```
