# Aggregation of single-cell signals

Aggregation of single-cell to pseudobulk data for non-count data.

## Usage

``` r
aggregateNonCountSignal(
  sce,
  assay = NULL,
  sample_id = NULL,
  cluster_id = NULL,
  min.cells = 10,
  min.signal = 0.01,
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

- sample_id:

  character string specifying which variable to use as sample id

- cluster_id:

  character string specifying which variable to use as cluster id

- min.cells:

  minimum number of observed cells for a sample to be included in the
  analysis

- min.signal:

  minimum signal value for a gene to be considered expressed in a
  sample. Proper value for this cutoff depends on the type of signal
  value

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

The `dreamlet` workflow can also be applied to non-count data. In this
case, a signal is averaged across all cells from a given sample and cell
type. Here `aggregateNonCountSignal()` performs the roles of
[`aggregateToPseudoBulk()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/aggregateToPseudoBulk.md)
followed by
[`processAssays()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/processAssays.md)
but using non-count data.

For each cell cluster, samples with at least `min.cells` are retained.
Only clusters with at least `min.samples` retained samples are kept.
Features are retained if they have at least `min.signal` in at least
min.prop fraction of the samples.

The precision of a measurement is the inverse of its sampling variance.
The precision weights are computed as `1/sem^2`, where
`sem = sd(signal) / sqrt(n)`, `signal` stores the values averaged across
cells, and `n` is the number of cells.

## Examples

``` r
library(muscat)
library(SingleCellExperiment)

data(example_sce)

# create pseudobulk for each sample and cell cluster
# using non-count signal
pb.signal <- aggregateNonCountSignal(example_sce,
  assay = "logcounts",
  cluster_id = "cluster_id",
  sample_id = "sample_id",
  verbose = FALSE
)

# Differential expression analysis within each assay,
# evaluated on the voom normalized data
res.dl <- dreamlet(pb.signal, ~group_id)
#>   B cells...
#> 0.4 secs
#>   CD14+ Monocytes...
#> 0.11 secs
#>   CD4 T cells...
#> 0.075 secs
#>   CD8 T cells...
#> 0.075 secs
#>   FCGR3A+ Monocytes...
#> 0.08 secs
```
