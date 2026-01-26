# Drop redundant terms from the model

Detect co-linear fixed effects and drop the last one

## Usage

``` r
dropRedundantTerms(formula, data, tol = 0.001)
```

## Arguments

- formula:

  original formula

- data:

  data.frame

- tol:

  tolerance to test difference of correlation from 1 or -1

## Value

a formula, possibly with terms omitted.

## Examples

``` r

# Valid formula
dropRedundantTerms(~ group + extra, sleep)
#> ~group + extra
#> <environment: 0x357907cf0>
```
