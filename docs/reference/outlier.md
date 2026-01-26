# Multivariate outlier detection

Detect multivariante outliers using Mahalanobis distance using mean and
covariance estimated either with standard or robust methods.

## Usage

``` r
outlier(data, robust = FALSE, ...)
```

## Arguments

- data:

  matrix of data

- robust:

  use robust covariance method, defaults to `FALSE`

- ...:

  arguments passed to
  [`MASS::cov.rob()`](https://rdrr.io/pkg/MASS/man/cov.rob.html)

## Value

`data.frame` storing chisq and z-score for each entry indicating
deviation from the mean. The z-score is computed by evaluating the
p-value of chisq statistic and converting it into a z-score

## Details

The distance follow a chisq distrubtion under the null with standard
method for mean and covariance. It is approximate if the robust method
is used. So use `qchisq(p = 0.999 , df = k)` to get cutoff to keep 99.9%
of samples under the null for data with `k=2` columns.

## Examples

``` r
data <- matrix(rnorm(200), 100, 2)

res <- outlier(data)

res[1:4,]
#>       chisq         z     pValue
#> 1 0.6375144 0.3490496 0.72705207
#> 2 1.6240517 0.7655271 0.44395776
#> 3 8.0789883 2.3737939 0.01760638
#> 4 5.8063084 1.9200630 0.05484994
```
