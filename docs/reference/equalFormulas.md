# Check if two formulas are equal

Check if two formulas are equal by evaluating the formulas and
extracting terms

## Usage

``` r
equalFormulas(formula1, formula2)
```

## Arguments

- formula1:

  first formula

- formula2:

  second formula

## Value

boolean value indciating of formulas are equivalent

## Examples

``` r

# These formulas are equivalent
formula1 <- ~ Size + 1
formula2 <- ~ 1 + Size

dreamlet:::equalFormulas(formula1, formula2)
#> [1] TRUE
```
