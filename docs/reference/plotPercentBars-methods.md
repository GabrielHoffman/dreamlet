# Bar plot of variance fractions

Bar plot of variance fractions for a subset of genes

## Usage

``` r
# S4 method for class 'vpDF'
plotPercentBars(
  x,
  col = c(ggColorHue(ncol(x) - 3), "grey85"),
  genes = unique(x$gene),
  width = NULL,
  ncol = 3,
  ...
)

# S4 method for class 'cellSpecificityValues'
plotPercentBars(
  x,
  col = ggColorHue(ncol(x)),
  genes = rownames(x),
  width = NULL,
  ...
)
```

## Arguments

- x:

  `vpDF` object returned by
  [`fitVarPart()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/fitVarPart.md)

- col:

  color of bars for each variable

- genes:

  name of genes to plot

- width:

  specify width of bars

- ncol:

  number of columns in the plot

- ...:

  other arguments

## Value

Bar plot showing variance fractions for each gene

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
#> 0.063 secs
#>   CD8 T cells...
#> 0.038 secs
#>   FCGR3A+ Monocytes...
#> 0.083 secs

# variance partitioning analysis
vp <- fitVarPart(res.proc, ~group_id)
#>   B cells...
#> 0.6 secs
#>   CD14+ Monocytes...
#> 0.84 secs
#>   CD4 T cells...
#> 0.64 secs
#>   CD8 T cells...
#> 0.4 secs
#>   FCGR3A+ Monocytes...
#> 0.77 secs
#> 

# Show variance fractions at the gene-level for each cell type
plotPercentBars(vp, genes = vp$gene[2:4], ncol = 2)
#> Warning: Ignoring empty aesthetic: `width`.

```
