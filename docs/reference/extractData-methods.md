# Extract normalized expression and `colData`

Extract normalized expression and `colData`

Extract normalized (i.e. log2 CPM) expression and `colData` from
`dreamletProcessedData`

## Usage

``` r
extractData(x, assay, cols = colnames(colData(x)), genes = rownames(x))

# S4 method for class 'dreamletProcessedData,character'
extractData(
  x,
  assay,
  cols = colnames(colData(x)),
  genes = rownames(assay(x, assay))
)
```

## Arguments

- x:

  `dreamletProcessedData` object

- assay:

  assay to extract

- cols:

  columns in `colData(x)` to extract. defaults to all columns as
  `colnames(colData(x))`

- genes:

  genes to extract from `assay(x, assay)$E`. defaults to all genes as
  `rownames(x)`

## Value

`data.frame` or `DataFrame` of merged expression and colData

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
#> 0.084 secs

# Extract all:
# Extract tibble of colData merged with expression.
# variables and genes are stored as columns, samples as rows
df_merge <- extractData(res.proc, "B cells")

# first few columns
df_merge[, 1:6]
#> # A tibble: 4 × 6
#>   Row.names group_id ISG15 AURKAIP1 MRPL20 SSU72
#>   <chr>     <fct>    <dbl>    <dbl>  <dbl> <dbl>
#> 1 ctrl101   ctrl      7.02     7.80   7.02  8.16
#> 2 ctrl107   ctrl      7.29     8.35   7.29  8.55
#> 3 stim101   stim     13.4      8.28   7.58  7.94
#> 4 stim107   stim     13.2      7.79   6.85  8.43

# Extract subset:
df_merge <- extractData(res.proc, "B cells", cols = "group_id", genes = c("SSU72", "U2AF1"))

df_merge
#> # A tibble: 4 × 4
#>   Row.names group_id SSU72 U2AF1
#>   <chr>     <fct>    <dbl> <dbl>
#> 1 ctrl101   ctrl      8.16  9.19
#> 2 ctrl107   ctrl      8.55  9.32
#> 3 stim101   stim      7.94  9.02
#> 4 stim107   stim      8.43  9.08

# Boxplot of expression
boxplot(SSU72 ~ group_id, df_merge)

#
```
