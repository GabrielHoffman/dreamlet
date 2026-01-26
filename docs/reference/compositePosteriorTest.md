# Perform composite test on results from mashr

The posterior probabilities for all genes and conditions is obtained as
`1-lFSR`. Let `prob` be an array storing results for one gene. The
probability that \_no\_ conditions in the exclusion set are non-zero is
`prod(1 - prob[exclude])`. The probability that \_all\_ conditions in
the inclusion set are non-zero is `prod(prob[include])`. The probability
that \_at least one\_ condition in the inclusion set is non-zero is
`1 - prod(1 - prob[include])`. The composite test is the product of the
probabilties computed from the inclusion and exclusion sets.

## Usage

``` r
compositePosteriorTest(
  x,
  include,
  exclude = NULL,
  test = c("at least 1", "all")
)
```

## Arguments

- x:

  `"dreamlet_mash_result"` from
  [`run_mash()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/run_mash.md)

- include:

  array of conditions in the inclusion set

- exclude:

  array of conditions in the exclusion set. Defaults to `NULL` for no
  exclusion

- test:

  evaluate the posterior probability of a non-zero effect in
  `"at least 1"` or `"all"` conditions

## Details

Perform composite test evaluating the specificity of an effect. Evalute
the posterior probability that an a non-zero effect present in \_all\_
or \_at least one\_ condition in the inclusion set, but \_no
conditions\_ in the exclusion set.

## See also

[`run_mash()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/run_mash.md)

## Examples

``` r
library(muscat)
library(mashr)
#> Loading required package: ashr
#> 
#> Attaching package: ‘ashr’
#> The following object is masked from ‘package:generics’:
#> 
#>     prune
library(SingleCellExperiment)

data(example_sce)

# create pseudobulk for each sample and cell cluster
pb <- aggregateToPseudoBulk(example_sce[1:100, ],
  assay = "counts",
  cluster_id = "cluster_id",
  sample_id = "sample_id",
  verbose = FALSE
)

# voom-style normalization
res.proc <- processAssays(pb, ~group_id)
#>   B cells...
#> 0.015 secs
#>   CD14+ Monocytes...
#> 0.017 secs
#>   CD4 T cells...
#> 0.018 secs
#>   CD8 T cells...
#> 0.016 secs
#>   FCGR3A+ Monocytes...
#> 0.016 secs

# Differential expression analysis within each assay,
# evaluated on the voom normalized data
res.dl <- dreamlet(res.proc, ~group_id)
#>   B cells...
#> 0.0089 secs
#>   CD14+ Monocytes...
#> 0.0097 secs
#>   CD4 T cells...
#> 0.0087 secs
#>   CD8 T cells...
#> 0.015 secs
#>   FCGR3A+ Monocytes...
#> 0.0094 secs

# run MASH model
# This can take 10s of minutes on real data
# This small datasets should take ~30s
res_mash <- run_mash(res.dl, "group_idstim")

# Composite test based on posterior probabilities
# to identify effect present in *at least 1* monocyte type
# and *NO* T-cell type.
include <- c("CD14+ Monocytes", "FCGR3A+ Monocytes")
exclude <- c("CD4 T cells", "CD8 T cells")

# Perform composite test
prob <- compositePosteriorTest(res_mash, include, exclude)

# examine the lFSR for top gene
get_lfsr(res_mash$model)[which.max(prob), , drop = FALSE]
#>       B cells CD14+ Monocytes CD4 T cells CD8 T cells FCGR3A+ Monocytes
#> FBXO6      NA              NA          NA          NA      1.981561e-17

# Test if *all* cell types have non-zero effect
prob <- compositePosteriorTest(res_mash, assayNames(res.dl))
```
