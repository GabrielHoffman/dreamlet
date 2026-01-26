# Aggregation of single-cell to pseudobulk data

Aggregation of single-cell to pseudobulk data. Adapted from
[`muscat::aggregateData`](https://rdrr.io/pkg/muscat/man/aggregateData.html)
and has same syntax and results. But can be much faster for
`SingleCellExperiment` backed by H5AD files using on-disk storage.

## Usage

``` r
aggregateToPseudoBulk(
  x,
  assay = NULL,
  sample_id = NULL,
  cluster_id = NULL,
  fun = c("sum", "mean", "median", "prop.detected", "num.detected", "sem", "number"),
  scale = FALSE,
  verbose = TRUE,
  BPPARAM = SerialParam(progressbar = verbose),
  checkValues = TRUE,
  h5adBlockSizes = 1e+09
)
```

## Arguments

- x:

  a
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html).

- assay:

  character string specifying the assay slot to use as input data.
  Defaults to the 1st available (`assayNames(x)[1]`).

- sample_id:

  character string specifying which variable to use as sample id

- cluster_id:

  character string specifying which variable to use as cluster id

- fun:

  a character string. Specifies the function to use as summary
  statistic. Passed to `summarizeAssayByGroup2`.

- scale:

  logical. Should pseudo-bulks be scaled with the effective library size
  & multiplied by 1M?

- verbose:

  logical. Should information on progress be reported?

- BPPARAM:

  a
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object specifying how aggregation should be parallelized.

- checkValues:

  logical. Should we check that signal values are positive integers?

- h5adBlockSizes:

  set the automatic block size block size (in bytes) for DelayedArray to
  read an H5AD file. Larger values use more memory but are faster.

## Value

a
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html).

Aggregation parameters (`assay, by, fun, scaled`) are stored in
`metadata()$agg_pars`, where `by = c(cluster_id, sample_id)`. The number
of cells that were aggregated are accessible in `int_colData()$n_cells`.

## Details

Adapted from
[`muscat::aggregateData`](https://rdrr.io/pkg/muscat/man/aggregateData.html)
and has similar syntax and same results. This is much faster for
`SingleCellExperiment` backed by H5AD files using `DelayedMatrix`
because this summarizes counts using
[`DelayedMatrixStats`](https://rdrr.io/pkg/DelayedMatrixStats/man/DelayedMatrixStats.html).
But this function also includes optmizations for `sparseMatrix` used by
[`Seurat`](https://satijalab.org/seurat/reference/Seurat-package.html)
by using `sparseMatrixStats`.

Keeps variables from
[`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
that are constant within `sample_id`. For example, sex will be constant
for all cells from the same `sample_id`, so it is retained as a variable
in the pseudobulk result. But number of expressed genes varies across
cells within each `sample_id`, so it is dropped from
[`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html).
Instead the mean value per cell type is stored in
`metadata(pb)$aggr_means`, and these can be included in regression
formulas downstream. In that case, the value of the covariates used per
sample will depend on the cell type analyzed.

## References

Crowell, HL, Soneson, C, Germain, P-L, Calini, D, Collin, L, Raposo, C,
Malhotra, D & Robinson, MD: Muscat detects subpopulation-specific state
transitions from multi-sample multi-condition single-cell
transcriptomics data. *Nature Communications* **11(1):6077** (2020).
doi: <https://doi.org/10.1038/s41467-020-19894-4>

## Author

Gabriel Hoffman, Helena L Crowell & Mark D Robinson

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

# pseudobulk data from each cell type
# is stored as its own assay
pb
#> class: SingleCellExperiment 
#> dim: 1267 4 
#> metadata(3): experiment_info agg_pars aggr_means
#> assays(5): B cells CD14+ Monocytes CD4 T cells CD8 T cells FCGR3A+
#>   Monocytes
#> rownames(1267): HES4 ISG15 ... CSTB PRMT2
#> rowData names(2): ENSEMBL SYMBOL
#> colnames(4): ctrl101 ctrl107 stim101 stim107
#> colData names(1): group_id
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):

# aggregate by cluster only,
# collapsing all samples into the same pseudobulk
pb2 <- aggregateToPseudoBulk(example_sce, 
 cluster_id = "cluster_id", 
 verbose = FALSE)

pb2
#> class: SingleCellExperiment 
#> dim: 1267 5 
#> metadata(3): experiment_info agg_pars aggr_means
#> assays(1): ''
#> rownames(1267): HES4 ISG15 ... CSTB PRMT2
#> rowData names(2): ENSEMBL SYMBOL
#> colnames(5): B cells CD14+ Monocytes CD4 T cells CD8 T cells FCGR3A+
#>   Monocytes
#> colData names(0):
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#
```
