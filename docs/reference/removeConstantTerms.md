# Remove constant terms from formula

Remove constant terms from formula. Also remove categorical variables
with a max of one example per category

## Usage

``` r
removeConstantTerms(formula, data)
```

## Arguments

- formula:

  original formula

- data:

  data.frame

## Value

a formula, possibly with terms omitted.

## Details

Adapted from `MoEClust::drop_constants`

## Examples

``` r

# Valid formula
removeConstantTerms(~ group + extra, sleep)
#> ~group + extra
#> <environment: 0x3021a7eb0>

# there is no variation in 'group' in this dataset
removeConstantTerms(~ group + extra, sleep[1:3, ])
#> ~1 + extra
#> <environment: 0x3021a7eb0>
```
