# Test if coefficient is different from a specified value

Test if coefficient is different from a specified value

## Usage

``` r
# S4 method for class 'dreamletResult'
getTreat(fit, lfc = log2(1.2), coef = NULL, number = 10, sort.by = "p")
```

## Arguments

- fit:

  dreamletResult object

- lfc:

  a minimum log2-fold-change below which changes not considered
  scientifically meaningful

- coef:

  which coefficient to test

- number:

  number of genes to return

- sort.by:

  column to sort by

## Value

`DataFrame` storing hypothesis test for each gene and cell type

## See also

[`limma::topTreat()`](https://rdrr.io/pkg/limma/man/toptable.html),
[`variancePartition::getTreat()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/getTreat-method.md)

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
#> 0.055 secs
#>   CD14+ Monocytes...
#> 0.098 secs
#>   CD4 T cells...
#> 0.062 secs
#>   CD8 T cells...
#> 0.037 secs
#>   FCGR3A+ Monocytes...
#> 0.081 secs

# Differential expression analysis within each assay,
# evaluated on the voom normalized data
res.dl <- dreamlet(res.proc, ~group_id)
#>   B cells...
#> 0.049 secs
#>   CD14+ Monocytes...
#> 0.066 secs
#>   CD4 T cells...
#> 0.062 secs
#>   CD8 T cells...
#> 0.034 secs
#>   FCGR3A+ Monocytes...
#> 0.062 secs

# show coefficients estimated for each cell type
coefNames(res.dl)
#> [1] "(Intercept)"  "group_idstim"

# extract results using limma-style syntax
# combines all cell types together
# adj.P.Val gives study-wide FDR
getTreat(res.dl, coef = "group_idstim", number = 3)
#> DataFrame with 3 rows and 7 columns
#>               assay     logFC   AveExpr         t     P.Value   adj.P.Val
#>         <character> <numeric> <numeric> <numeric>   <numeric>   <numeric>
#> 1           B cells   6.17666   10.2306   18.2945 1.75997e-14 7.90401e-11
#> 2 FCGR3A+ Monocytes   5.27610   11.9149   26.0100 9.11216e-14 2.04613e-10
#> 3 FCGR3A+ Monocytes   7.43146   11.1948   22.7202 6.79966e-13 9.94753e-10
#>           B
#>   <numeric>
#> 1   23.2102
#> 2   22.0773
#> 3   19.7485
```
