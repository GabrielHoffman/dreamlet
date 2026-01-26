# Error handling

[`dreamlet()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/dreamlet.md)
evaluates precision-weighted linear (mixed) models on each gene that
passes standard filters. The linear mixed model used by `dream()` can be
a little fragile for small sample sizes and correlated covariates.
[`dreamlet()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/dreamlet.md)
runs
[`variancePartition::dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)
in the backend for each cell cluster. `dream()` reports model failures
for each cell cluster and
[`dreamlet()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/dreamlet.md)
reports these failures to the user.
[`dreamlet()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/dreamlet.md)
returns all **successful** model fits to be used for downstream
analysis.

See details from `variancePartition` [error
page](https://diseaseneurogenomics.github.io/variancePartition/articles/errors.html).

## Errors with random effects

Due to a recent [bug](https://github.com/lme4/lme4/issues/763) in the
dependency `Matrix` package, all random effects models may fail for
technical reasons. If your random effects analysis is failing for all
genes in cases with no good explanation, this bug may be responsible.
This case can be detected and resolved as follows:

``` r

library(lme4)

# Fit simple mixed model
lmer(Reaction ~ (1 | Subject), sleepstudy)
# Error in initializePtr() : 
#  function 'chm_factor_ldetL2' not provided by package 'Matrix'
```

This error indicates incompatible installs of `Matrix` and `lme4`. This
can be solved with

``` r

install.packages("lme4", type = "source") 
```

followed by **restarting** R.

## Errors at the assay- and gene-level

The most common issue is when
[`dreamlet()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/dreamlet.md)
analysis succeeds for most genes, but a handful of genes fail in each
cell cluster. These genes can fail if the iterative process of fitting
the linear mixed model does not converge, or if the estimated covariance
matrix that is supposed be positive definite has an eigen-value that is
negative or too close to zero due to rounding errors in floating point
arithmetic.

In these cases,
[`dreamlet()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/dreamlet.md)
stores a summary of these failures for all cell clusters that is
accessible with
[`details()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/details-methods.md).
Specific failure messages for each cell cluster and gene can be
extracted using
[`seeErrors()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/seeErrors-methods.md)

Here we demonstrate how
[`dreamlet()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/dreamlet.md)
handles model failures:

``` r

library(dreamlet)
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

# voom-style normalization for each cell cluster
res.proc <- processAssays(
  pb[1:300, ],
  ~group_id
)

# Redundant formula
# This example is an extreme example of redundancy
# but more subtle cases often show up in real data
form <- ~ group_id + (1 | group_id)

# fit dreamlet model
res.dl <- dreamlet(res.proc, form)
##  B cells...7.9 secs
##  CD14+ Monocytes...10 secs
##  CD4 T cells...9 secs
##  CD8 T cells...4.4 secs
##  FCGR3A+ Monocytes...11 secs
##
## Of 1,062 models fit across all assays, 96.2% failed

# summary of models
res.dl
## class: dreamletResult
## assays(5): B cells CD14+ Monocytes CD4 T cells CD8 T cells FCGR3A+ Monocytes
## Genes:
##  min: 3
##  max: 11
## details(7): assay n_retain ... n_errors error_initial
## coefNames(2): (Intercept) group_idstim
##
## Of 1,062 models fit across all assays, 96.2% failed

# summary of models for each cell cluster
details(res.dl)
##               assay n_retain                    formula formDropsTerms n_genes n_errors error_initial
## 1           B cells        4 ~group_id + (1 | group_id)          FALSE     201      190         FALSE
## 2   CD14+ Monocytes        4 ~group_id + (1 | group_id)          FALSE     269      263         FALSE
## 3       CD4 T cells        4 ~group_id + (1 | group_id)          FALSE     216      207         FALSE
## 4       CD8 T cells        4 ~group_id + (1 | group_id)          FALSE     118      115         FALSE
## 5 FCGR3A+ Monocytes        4 ~group_id + (1 | group_id)          FALSE     258      247         FALSE
```

- `assay`: cell type
- `n_retain`: number of samples retained
- `formula`: regression formula used after variable filtering
- `formDropsTerms`: whether a variable was dropped from the formula for
  having zero variance following filtering
- `n_genes`: number of genes analyzed
- `n_errors`: number of genes with errors
- `error_initial`: indicator for assay-level error

### Assay-level errors

Before the full dataset is analyzed,
[`dreamlet()`](http://DiseaseNeurogenomics.github.io/dreamlet/reference/dreamlet.md)
runs a test for each assay to see if the model succeeds. If the model
fails, its does not continue analysis for that assay. These assay-level
errors are reported above in the `error_initial` column, and details are
returned here.

``` r

# Extract errors as a tibble
res.err = seeErrors(res.dl)
##   Assay-level errors: 0
##   Gene-level errors: 1038

# No errors at the assay level
res.err$assayLevel

# the most common error is:
"Some predictor variables are on very different scales: consider rescaling"
```

This indicates that the scale of the predictor variables are very
different and can affect the numerical stability of the iterative
algorithm. This can be solved by running
[`scale()`](https://rdrr.io/r/base/scale.html) on each variable in the
formula:

``` r

form = ~ scale(x) + scale(y) + ...
```

### Gene-level errors

A model can fail for a single gene if covariates are too correlated, or
for other numerical issues. Failed models are reported here and are not
included in downstream analysis.

``` r

# See gene-level errors for each assay
res.err$geneLevel[1:2,]
## # A tibble: 2 × 3
##   assay   feature  errorText
##   <chr>   <chr>    <chr>                               
## B cells ISG15    "Error in lmerTest:::as_lmerModLT(model, devfun, tol = tol):…
## B cells AURKAIP1 "Error in lmerTest:::as_lmerModLT(model, devfun, tol = tol):…

# See full error message text
res.err$geneLevel$errorText[1]
"Error in lmerTest:::as_lmerModLT(model, devfun, tol = tol): (converted from warning) 
Model may not have converged with 1 eigenvalue close to zero: 1.4e-09\n"
```

This message indicates that the model was numerically unstable likely
because the variables are closely correlated.

## Session Info

    ## R version 4.5.1 (2025-06-13)
    ## Platform: aarch64-apple-darwin23.6.0
    ## Running under: macOS Sonoma 14.7.1
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /opt/homebrew/Cellar/openblas/0.3.30/lib/libopenblasp-r0.3.30.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
    ##  [5] xfun_0.56         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
    ##  [9] rmarkdown_2.30    lifecycle_1.0.5   cli_3.6.5         sass_0.4.10      
    ## [13] pkgdown_2.2.0     textshaping_1.0.4 jquerylib_0.1.4   systemfonts_1.3.1
    ## [17] compiler_4.5.1    tools_4.5.1       ragg_1.5.0        bslib_0.9.0      
    ## [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0   
    ## [25] rlang_1.1.7       fs_1.6.6          htmlwidgets_1.6.4
