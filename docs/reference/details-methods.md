# Extract details from dreamletProcessedData

Extract details from `dreamletProcessedData`

## Usage

``` r
details(object)

# S4 method for class 'dreamletProcessedData'
details(object)

# S4 method for class 'dreamletResult'
details(object)

# S4 method for class 'vpDF'
details(object)
```

## Arguments

- object:

  A `dreamletProcessedData` object

## Value

Extract detailed information from some classes

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
#> 0.086 secs
#>   CD4 T cells...
#> 0.062 secs
#>   CD8 T cells...
#> 0.038 secs
#>   FCGR3A+ Monocytes...
#> 0.082 secs

# For each cell type, number of samples retained,
# and variables retained
details(res.proc)
#>               assay n_retain   formula formDropsTerms n_genes n_errors
#> 1           B cells        4 ~group_id          FALSE     847        0
#> 2   CD14+ Monocytes        4 ~group_id          FALSE    1130        0
#> 3       CD4 T cells        4 ~group_id          FALSE     897        0
#> 4       CD8 T cells        4 ~group_id          FALSE     531        0
#> 5 FCGR3A+ Monocytes        4 ~group_id          FALSE    1086        0
#>   error_initial
#> 1         FALSE
#> 2         FALSE
#> 3         FALSE
#> 4         FALSE
#> 5         FALSE
```
