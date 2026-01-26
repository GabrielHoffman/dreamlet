# Processing SingleCellExperiment to dreamletProcessedData

For raw counts, estimate precision weights using linear mixed model
weighting by number of cells observed for each sample. For normalized
data, only weight by number of cells.

## Usage

``` r
processAssays(
  sceObj,
  formula,
  assays = assayNames(sceObj),
  min.cells = 5,
  min.count = 5,
  min.samples = 4,
  min.prop = 0.4,
  isCounts = TRUE,
  normalize.method = "TMM",
  span = "auto",
  quiet = FALSE,
  weightsList = NULL,
  BPPARAM = SerialParam(),
  ...
)
```

## Arguments

- sceObj:

  SingleCellExperiment object

- formula:

  regression formula for differential expression analysis

- assays:

  array of assay names to include in analysis. Defaults to
  `assayNames(sceObj)`

- min.cells:

  minimum number of observed cells for a sample to be included in the
  analysis

- min.count:

  used to compute a CPM threshold of
  `CPM.cutoff = min.count/median(lib.size)*1e6`. Passed to
  [`edgeR::filterByExpr()`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html)

- min.samples:

  minimum number of samples passing cutoffs for cell cluster to be
  retained

- min.prop:

  minimum proportion of retained samples with `CPM > CPM.cutoff`

- isCounts:

  logical, indicating if data is raw counts

- normalize.method:

  normalization method to be used by `calcNormFactors`

- span:

  Lowess smoothing parameter using by
  [`variancePartition::voomWithDreamWeights()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/voomWithDreamWeights.md)

- quiet:

  show messages

- weightsList:

  list storing matrix of precision weights for each cell type. If `NULL`
  precision weights are set to 1

- BPPARAM:

  parameters for parallel evaluation

- ...:

  other arguments passed to `dream`

## Value

Object of class `dreamletProcessedData` storing voom-style normalized
expression data

## Details

For each cell cluster, samples with at least `min.cells` are retained.
Only clusters with at least `min.samples` retained samples are kept.
Genes are retained if they have at least `min.count` reads in at least
`min.prop` fraction of the samples. Current values are reasonable
defaults, since genes that don't pass these cutoffs are very
underpowered for differential expression analysis and only increase the
multiple testing burden. But values of `min.cells = 2` and
`min.count = 2` are also reasonable to include more genes in the
analysis.

The precision weights are estimated using the residuals fit from the
specified formula. These weights are robust to changes in the formula as
long as the major variables explaining the highest fraction of the
variance are included.

If `weightsList` is `NULL`, precision weights are set to 1 internally.

## See also

`voomWithDreamWeights()`

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
#> 0.062 secs
#>   CD8 T cells...
#> 0.038 secs
#>   FCGR3A+ Monocytes...
#> 0.083 secs

# Differential expression analysis within each assay,
# evaluated on the voom normalized data
res.dl <- dreamlet(res.proc, ~group_id)
#>   B cells...
#> 0.049 secs
#>   CD14+ Monocytes...
#> 0.071 secs
#>   CD4 T cells...
#> 0.052 secs
#>   CD8 T cells...
#> 0.033 secs
#>   FCGR3A+ Monocytes...
#> 0.062 secs
#
```
