# Meta-analysis across multiple studies

Meta-analysis across multiple studies

## Usage

``` r
meta_analysis(
  x,
  method = "FE",
  group = c("ID", "assay"),
  control = list(maxiter = 2000)
)
```

## Arguments

- x:

  `data.frame` rbind'ing results across genes, cell types and datasets

- method:

  meta-analysis method. Values are fed into
  [`metafor::rma()`](https://wviechtb.github.io/metafor/reference/rma.uni.html),
  except for `'RE2C'` which calls
  [`remaCor::RE2C()`](http://gabrielhoffman.github.io/remaCor/reference/RE2C.md).

- group:

  colums in `x` to group by. For results from `dreamlet::topTable()`,
  results are aggregrated by gene and cell type (i.e. `'ID'` and
  `'assay'`). If `x` is not from this function, this argument allows the
  function to group results properly

- control:

  passed to `rma(..,control)`

## Details

- `'FE'`: :

  fixed effects meta-analysis

- `'REML'`: :

  random effects meta-analysis

- `'RE2C'`: :

  joint testing of fixed and random effects

## Examples

``` r
library(dreamlet)
library(muscat)

data(example_sce)

# create pseudobulk for each sample and cell cluster
pb <- aggregateToPseudoBulk(example_sce,
  assay = "counts",
  cluster_id = "cluster_id",
  sample_id = "sample_id",
  verbose = FALSE
)

# voom-style normalization
# just 'CD14+ Monocytes' for speed
res.proc <- processAssays(pb, ~group_id, assays = "CD14+ Monocytes")
#>   CD14+ Monocytes...
#> 0.12 secs

# dreamlet
res.dl <- dreamlet(res.proc, ~group_id)
#>   CD14+ Monocytes...
#> 0.065 secs

tab1 <- topTable(res.dl, coef = "group_idstim", number = Inf)
tab1$Dataset <- "1"

# Results from a second cohort
# Here, just a copy of the same results for simplicity
tab2 <- tab1
tab2$Dataset <- "2"

# rbind
tab_combined <- rbind(tab1, tab2)

# Perform fixed effects meta-analysis
res <- meta_analysis(tab_combined, method = "FE")

res[1:3, ]
#> # A tibble: 3 × 8
#> # Groups:   ID, assay [3]
#>   ID     assay           estimate std.error statistic  p.value n.studies method
#>   <chr>  <chr>              <dbl>     <dbl>     <dbl>    <dbl>     <int> <chr> 
#> 1 ABRACL CD14+ Monocytes    0.631     0.236      2.68 7.35e- 3         2 FE    
#> 2 ACOT9  CD14+ Monocytes    2.65      0.199     13.3  2.17e-40         2 FE    
#> 3 ACP5   CD14+ Monocytes   -0.220     0.140     -1.57 1.17e- 1         2 FE    
```
