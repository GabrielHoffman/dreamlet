# Bar plot of cell compositions

Bar plot of cell compositions

## Usage

``` r
plotCellComposition(obj, col, width = NULL)

# S4 method for class 'SingleCellExperiment'
plotCellComposition(obj, col, width = NULL)

# S4 method for class 'matrix'
plotCellComposition(obj, col, width = NULL)

# S4 method for class 'data.frame'
plotCellComposition(obj, col, width = NULL)
```

## Arguments

- obj:

  matrix of \[cells\] x \[samples\] or `SingleCellExperiment` from
  `aggregateToPseudoBulk`

- col:

  array of colors. If missing, use default colors. If `names(col)` is
  the same as `arrayNames(obj)`, then colors will be assigned by assay
  name#'

- width:

  specify width of bars

## Value

Barplot showing cell fractions

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

# show cell composition bar plots
plotCellComposition(pb)
#> Warning: Ignoring empty aesthetic: `width`.


# extract cell counts
df_cellCounts <- cellCounts(pb)

# show cell composition bar plots
plotCellComposition(df_cellCounts)
#> Warning: Ignoring empty aesthetic: `width`.

```
