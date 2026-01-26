# Table of Top Genes from dreamlet fit

Extract a table of the top-ranked genes from a dreamlet fit.

## Usage

``` r
# S4 method for class 'dreamletResult'
topTable(
  fit,
  coef = NULL,
  number = 10,
  genelist = NULL,
  adjust.method = "BH",
  sort.by = "P",
  resort.by = NULL,
  p.value = 1,
  lfc = 0,
  confint = FALSE
)
```

## Arguments

- fit:

  dreamletResult object

- coef:

  coef

- number:

  number

- genelist:

  genelist

- adjust.method:

  adjust.method

- sort.by:

  sort.by

- resort.by:

  resort.by

- p.value:

  p.value

- lfc:

  lfc

- confint:

  confint

## Value

`data.frame` storing hypothesis test for each gene and cell type

## See also

[`limma::topTable()`](https://rdrr.io/pkg/limma/man/toptable.html),
`variancePartition::topTable()`

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
#> 0.066 secs
#>   CD8 T cells...
#> 0.048 secs
#>   FCGR3A+ Monocytes...
#> 0.081 secs

# Differential expression analysis within each assay,
# evaluated on the voom normalized data
res.dl <- dreamlet(res.proc, ~group_id)
#>   B cells...
#> 0.049 secs
#>   CD14+ Monocytes...
#> 0.064 secs
#>   CD4 T cells...
#> 0.053 secs
#>   CD8 T cells...
#> 0.034 secs
#>   FCGR3A+ Monocytes...
#> 0.07 secs

# show coefficients estimated for each cell type
coefNames(res.dl)
#> [1] "(Intercept)"  "group_idstim"

# extract results using limma-style syntax
# combines all cell types together
# adj.P.Val gives study-wide FDR
topTable(res.dl, coef = "group_idstim", number = 3)
#> DataFrame with 3 rows and 9 columns
#>               assay          ID     logFC   AveExpr         t     P.Value
#>         <character> <character> <numeric> <numeric> <numeric>   <numeric>
#> 1           B cells       ISG15   6.17666   10.2306   19.1083 1.26348e-14
#> 2 FCGR3A+ Monocytes      CXCL10   5.27610   11.9149   27.3747 7.09630e-14
#> 3           B cells       ISG20   3.60083   11.5794   16.0129 3.92209e-13
#>     adj.P.Val         B     z.std
#>     <numeric> <numeric> <numeric>
#> 1 5.67431e-11   23.2102   19.1083
#> 2 1.59347e-10   22.0773   27.3747
#> 3 5.87136e-10   20.2656   16.0129
```
