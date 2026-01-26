# Extract residuals from `dreamletResult`

Extract residuals from `dreamletResult`

## Usage

``` r
# S4 method for class 'dreamletResult'
residuals(object, y, ..., type = c("response", "pearson"))
```

## Arguments

- object:

  `dreamletResult` object

- y:

  `dreamletProcessedData` object

- ...:

  other arguments

- type:

  compute either `"response"` residuals or `"pearson"` residuals.

## Value

residuals from model fit

## Details

`"response"` residuals are the typical residuals returned from
[`lm()`](https://rdrr.io/r/stats/lm.html). `"pearson"` residuals divides
each residual value by its estimated standard error. This requires
specifying `y`

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
#> 0.063 secs
#>   CD8 T cells...
#> 0.038 secs
#>   FCGR3A+ Monocytes...
#> 0.09 secs

# Differential expression analysis within each assay,
# evaluated on the voom normalized data
res.dl <- dreamlet(res.proc, ~group_id)
#>   B cells...
#> 0.049 secs
#>   CD14+ Monocytes...
#> 0.064 secs
#>   CD4 T cells...
#> 0.052 secs
#>   CD8 T cells...
#> 0.033 secs
#>   FCGR3A+ Monocytes...
#> 0.069 secs

# extract typical residuals for each assay (i.e. cell type)
# Return list with entry for each assay with for retained samples and genes
resid.lst <- residuals(res.dl)

# Get Pearson residuals:
# typical residuals scaled by the standard deviation
residPearson.lst <- residuals(res.dl, res.proc, type = "pearson")
```
