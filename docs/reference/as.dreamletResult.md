# Convert list of regression fits to `dreamletResult`

Convert list of regression fits to `dreamletResult` for downstream
analysis

## Usage

``` r
as.dreamletResult(fitList, df_details = NULL)
```

## Arguments

- fitList:

  list of regression fit with `dream()`

- df_details:

  `data.frame` storing assay details

## Value

object of class `dreamletResult`

## Details

Useful for combining multiple runs of
[`dreamletCompareClusters()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/dreamletCompareClusters.md)
into a single `dreamletResult` for downstream analysis

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

# first comparison
ct.pairs <- c("B cells", "CD14+ Monocytes")
fit <- dreamletCompareClusters(pb, ct.pairs, method = "fixed")
#> Initial filtering...
#> Filtering for paired samples...

# second comparison
ct.pairs2 <- c("B cells", "CD8 T cells")
fit2 <- dreamletCompareClusters(pb, ct.pairs2, method = "fixed")
#> Initial filtering...
#> Filtering for paired samples...

# Make a list storing each result with a meaningful name
fitList <- list()

id <- paste0("[", ct.pairs[1], "]_vs_[", ct.pairs[2], "]")
fitList[[id]] <- fit

id <- paste0("[", ct.pairs2[1], "]_vs_[", ct.pairs2[2], "]")
fitList[[id]] <- fit2

# create a dreamletResult form this list
res.compare <- as.dreamletResult(fitList)
res.compare
#> class: dreamletResult 
#> assays(2): [B cells]_vs_[CD14+ Monocytes] [B cells]_vs_[CD8 T cells]
#> Genes:
#>  min: 351 
#>  max: 633 
#> details(0):
#> coefNames(6): compare cellClusterbaseline ... Samplestim101
#>   Samplestim107
```
