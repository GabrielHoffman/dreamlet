# Check variables in a formula

Check that variables in formula are present in the data

## Usage

``` r
checkFormula(formula, data)
```

## Arguments

- formula:

  formula of variables to check

- data:

  data.frame storing variables in the formula

## Value

If formula is valid, return TRUE. Else throw error

## Examples

``` r

# Valid formula
dreamlet:::checkFormula(~speed, cars)

# Not valid formula
# dreamlet:::checkFormula( ~ speed + a, cars)
```
