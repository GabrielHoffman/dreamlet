# Hierarchical clustering on cell types from pseudobulk

Perform hierarchical clustering on cell types from pseudobulk by
aggregating read counts from each cell type.

## Usage

``` r
buildClusterTreeFromPB(
  pb,
  method = c("complete", "ward.D", "single", "average", "mcquitty", "median", "centroid",
    "ward.D2"),
  dist.method = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
  assays = assayNames(pb)
)
```

## Arguments

- pb:

  `SingleCellObject` storing pseudobulk for each cell type in in
  [`assay()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  field

- method:

  clustering method for
  [`hclust()`](https://rdrr.io/r/stats/hclust.html)

- dist.method:

  distance metric

- assays:

  which assays to include

## Value

hierarchical clustering object of class `hclust`

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

# Hierarchical clustering of cell types
hcl <- buildClusterTreeFromPB(pb)

plot(hcl)

```
