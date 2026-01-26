# Violin plot of variance fractions

Violin plot of variance fraction for each gene and each variable

## Usage

``` r
# S4 method for class 'DataFrame'
plotVarPart(
  obj,
  col = c(ggColorHue(base::ncol(obj) - 3), "grey85"),
  label.angle = 20,
  main = "",
  ylab = "",
  convertToPercent = TRUE,
  ncol = 3,
  ...
)
```

## Arguments

- obj:

  `varParFrac` object returned by `fitExtractVarPart` or
  `extractVarPart`

- col:

  vector of colors

- label.angle:

  angle of labels on x-axis

- main:

  title of plot

- ylab:

  text on y-axis

- convertToPercent:

  multiply fractions by 100 to convert to percent values

- ncol:

  number of columns in the plot

- ...:

  additional arguments

## Value

Violin plot showing variance fractions

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
#> 0.089 secs
#>   CD4 T cells...
#> 0.064 secs
#>   CD8 T cells...
#> 0.038 secs
#>   FCGR3A+ Monocytes...
#> 0.093 secs

# variance partitioning analysis
vp <- fitVarPart(res.proc, ~group_id)
#>   B cells...
#> 0.65 secs
#>   CD14+ Monocytes...
#> 0.84 secs
#>   CD4 T cells...
#> 0.66 secs
#>   CD8 T cells...
#> 0.4 secs
#>   FCGR3A+ Monocytes...
#> 0.86 secs
#> 

# Summarize variance fractions genome-wide for each cell type
plotVarPart(vp)

```
