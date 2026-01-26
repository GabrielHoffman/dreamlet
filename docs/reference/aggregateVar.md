# Per-sample variance of single-cell counts

Aggregation function for single-cell log-normalized counts to calculate
per-sample variance for dreamlet.

## Usage

``` r
aggregateVar(
  sce,
  assay = NULL,
  cluster_id = NULL,
  sample_id = NULL,
  min.cells = 10,
  min.var = 0.01,
  min.samples = 4,
  min.prop = 0.4,
  verbose = TRUE,
  BPPARAM = SerialParam(progressbar = verbose)
)
```

## Arguments

- sce:

  a
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html).

- assay:

  character string specifying the assay slot to use as input data.
  Defaults to the 1st available (`assayNames(x)[1]`).

- cluster_id:

  character string specifying which variable to use as cluster id

- sample_id:

  character string specifying which variable to use as sample id

- min.cells:

  minimum number of observed cells for a sample to be included in the
  analysis

- min.var:

  minimum variance for a gene to be considered expressed in a sample

- min.samples:

  minimum number of samples passing cutoffs for cell cluster to be
  retained

- min.prop:

  minimum proportion of retained samples with non-zero counts for a gene
  to be

- verbose:

  logical. Should information on progress be reported?

- BPPARAM:

  a
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object specifying how aggregation should be parallelized.

## Value

a `dreamletProcessedData` object

## Details

The `dreamlet` workflow can also be applied to model gene expression
variance. In this case, a per-sample per-gene variance is calculated
across all cells from a given sample and cell type. Here
`aggregateVar()` performs the roles of
[`aggregateToPseudoBulk()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/aggregateToPseudoBulk.md)
followed by
[`processAssays()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/processAssays.md)
but using log-normalized count data.

For each cell cluster, samples with at least min.cells are retained.
Only clusters with at least min.samples retained samples are kept.
Features are retained if they have at least min.var in at least min.prop
fraction of the samples.

The precision of a measurement is the inverse of its sampling variance.
The precision weights are computed as `1/sem^2`, where
`sem = sd / sqrt(n)` and `n` is the number of cells.

## Examples

``` r
library(muscat)
library(SingleCellExperiment)

data(example_sce)

# Compute variance for each sample and cell cluster
pbVar <- aggregateVar(example_sce,
  assay = "counts",
  cluster_id = "cluster_id",
  sample_id = "sample_id",
  verbose = FALSE
)
```
