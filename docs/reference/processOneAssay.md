# Processing expression data from assay

For raw counts, filter genes and samples, then estimate precision
weights using linear mixed model weighting by number of cells observed
for each sample. For normalized data, only weight by number of cells

## Usage

``` r
processOneAssay(
  y,
  formula,
  data,
  n.cells,
  min.cells = 5,
  min.count = 2,
  min.samples = 4,
  min.prop = 0.4,
  min.total.count = 15,
  isCounts = TRUE,
  normalize.method = "TMM",
  span = "auto",
  quiet = TRUE,
  weights = NULL,
  rescaleWeightsAfter = FALSE,
  BPPARAM = SerialParam(),
  ...
)
```

## Arguments

- y:

  matrix of counts or log2 CPM

- formula:

  regression formula for differential expression analysis

- data:

  metadata used in regression formula

- n.cells:

  array of cell count for each sample

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

- min.total.count:

  minimum total count required per gene for inclusion

- isCounts:

  logical, indicating if data is raw counts

- normalize.method:

  normalization method to be used by `calcNormFactors`

- span:

  Lowess smoothing parameter using by
  [`variancePartition::voomWithDreamWeights()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/voomWithDreamWeights.md)

- quiet:

  show messages

- weights:

  matrix of precision weights

- rescaleWeightsAfter:

  default = FALSE, should the output weights be scaled by the input
  weights

- BPPARAM:

  parameters for parallel evaluation

- ...:

  other arguments passed to `dream`

## Value

`EList` object storing log2 CPM and precision weights

## See also

[`processAssays()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/processAssays.md)
