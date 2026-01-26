# Forest plot

Forest plot

## Usage

``` r
plotForest(x, gene, coef, ...)

# S4 method for class 'dreamletResult'
plotForest(x, gene, coef, assays = names(x), ylim = NULL)

# S4 method for class 'dreamlet_mash_result'
plotForest(x, gene, coef, assays = colnames(x$logFC.original), ylim = NULL)
```

## Arguments

- x:

  result from `dreamlet`

- gene:

  gene to show results for

- coef:

  coefficient to test with `topTable`

- ...:

  other arguments

- assays:

  array of assays to plot

- ylim:

  limits for the y axis

## Value

Plot showing effect sizes

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
#> 0.087 secs
#>   CD4 T cells...
#> 0.064 secs
#>   CD8 T cells...
#> 0.039 secs
#>   FCGR3A+ Monocytes...
#> 0.49 secs

# Differential expression analysis within each assay,
# evaluated on the voom normalized data
res.dl <- dreamlet(res.proc, ~group_id)
#>   B cells...
#> 0.049 secs
#>   CD14+ Monocytes...
#> 0.064 secs
#>   CD4 T cells...
#> 0.051 secs
#>   CD8 T cells...
#> 0.032 secs
#>   FCGR3A+ Monocytes...
#> 0.061 secs

# show coefficients estimated for each cell type
coefNames(res.dl)
#> [1] "(Intercept)"  "group_idstim"

# Show estimated log fold change with in each cell type
plotForest(res.dl, gene = "ISG20", coef = "group_idstim")

```
