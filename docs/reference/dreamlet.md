# Differential expression for each assay

Perform differential expression for each assay using linear (mixed)
models

## Usage

``` r
dreamlet(
  x,
  formula,
  data = colData(x),
  assays = assayNames(x),
  contrasts = NULL,
  min.cells = 10,
  robust = FALSE,
  quiet = FALSE,
  BPPARAM = SerialParam(),
  use.eBayes = TRUE,
  ...
)

# S4 method for class 'dreamletProcessedData'
dreamlet(
  x,
  formula,
  data = colData(x),
  assays = assayNames(x),
  contrasts = NULL,
  min.cells = 10,
  robust = FALSE,
  quiet = FALSE,
  BPPARAM = SerialParam(),
  use.eBayes = TRUE,
  ...
)
```

## Arguments

- x:

  SingleCellExperiment or dreamletProcessedData object

- formula:

  regression formula for differential expression analysis

- data:

  metadata used in regression formula

- assays:

  array of assay names to include in analysis. Defaults to
  `assayNames(x)`

- contrasts:

  character vector specifying contrasts specifying linear combinations
  of fixed effects to test. This is fed into
  `makeContrastsDream( formula, data, contrasts=contrasts)`

- min.cells:

  minimum number of observed cells for a sample to be included in the
  analysis

- robust:

  logical, use eBayes method that is robust to outlier genes

- quiet:

  show messages

- BPPARAM:

  parameters for parallel evaluation

- use.eBayes:

  should `eBayes` be used on result? (defualt: TRUE)

- ...:

  other arguments passed to `dream`

## Value

Object of class `dreamletResult` storing results for each cell type

## Details

Fit linear (mixed) model on each cell type separately. For advanced use
of contrasts see
[`variancePartition::makeContrastsDream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/makeContrastsDream.md)
and vignette
<https://gabrielhoffman.github.io/variancePartition/articles/dream.html#advanced-hypothesis-testing-1>.

## See also

[`variancePartition::dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md),
[`variancePartition::makeContrastsDream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/makeContrastsDream.md)

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
#> 0.053 secs
#>   CD14+ Monocytes...
#> 0.086 secs
#>   CD4 T cells...
#> 0.063 secs
#>   CD8 T cells...
#> 0.038 secs
#>   FCGR3A+ Monocytes...
#> 0.083 secs

# Differential expression analysis within each assay,
# evaluated on the voom normalized data
res.dl <- dreamlet(res.proc, ~group_id)
#>   B cells...
#> 0.059 secs
#>   CD14+ Monocytes...
#> 0.064 secs
#>   CD4 T cells...
#> 0.053 secs
#>   CD8 T cells...
#> 0.033 secs
#>   FCGR3A+ Monocytes...
#> 0.063 secs

# Examine results
res.dl
#> class: dreamletResult 
#> assays(5): B cells CD14+ Monocytes CD4 T cells CD8 T cells FCGR3A+
#>   Monocytes
#> Genes:
#>  min: 531 
#>  max: 1130 
#> details(7): assay n_retain ... n_errors error_initial
#> coefNames(2): (Intercept) group_idstim

# Examine details for each assay
details(res.dl)
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
