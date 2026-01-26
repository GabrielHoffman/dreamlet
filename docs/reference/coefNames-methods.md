# Get coefficient names

Get coefficient names

## Usage

``` r
coefNames(obj)

# S4 method for class 'dreamletResult'
coefNames(obj)
```

## Arguments

- obj:

  A `dreamletResult` object

## Value

array storing names of coefficients

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
#> 0.065 secs
#>   CD8 T cells...
#> 0.042 secs
#>   FCGR3A+ Monocytes...
#> 0.087 secs

# Differential expression analysis within each assay,
# evaluated on the voom normalized data
res.dl <- dreamlet(res.proc, ~group_id)
#>   B cells...
#> 0.054 secs
#>   CD14+ Monocytes...
#> 0.076 secs
#>   CD4 T cells...
#> 0.052 secs
#>   CD8 T cells...
#> 0.033 secs
#>   FCGR3A+ Monocytes...
#> 0.062 secs

# show coefficients estimated for each cell type
coefNames(res.dl)
#> [1] "(Intercept)"  "group_idstim"
```
