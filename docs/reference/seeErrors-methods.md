# Get error text

Get error text

## Usage

``` r
seeErrors(obj)

# S4 method for class 'dreamletResult'
seeErrors(obj)

# S4 method for class 'dreamletProcessedData'
seeErrors(obj)

# S4 method for class 'vpDF'
seeErrors(obj)
```

## Arguments

- obj:

  A `dreamletResult` object

## Value

`tibble` storing error text

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
#> 0.11 secs
#>   CD4 T cells...
#> 0.063 secs
#>   CD8 T cells...
#> 0.047 secs
#>   FCGR3A+ Monocytes...
#> 0.081 secs

# Differential expression analysis within each assay,
# evaluated on the voom normalized data
res.dl <- dreamlet(res.proc, ~group_id)
#>   B cells...
#> 0.05 secs
#>   CD14+ Monocytes...
#> 0.064 secs
#>   CD4 T cells...
#> 0.053 secs
#>   CD8 T cells...
#> 0.034 secs
#>   FCGR3A+ Monocytes...
#> 0.072 secs

# show errors
# but none are reported
res.err = seeErrors(res.dl)
#>    Assay-level errors: 0
#>    Gene-level errors: 0
```
