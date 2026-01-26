# Convert results table to matrix

Convert results table to matrix

## Usage

``` r
tabToMatrix(tab, col, rn = "ID", cn = "assay")
```

## Arguments

- tab:

  results table from `topTable()`

- col:

  which column to extract

- rn:

  column id storing rownames

- cn:

  column id storing colnames

## Value

matrix storing values of column `col` in rows defind by `rn` and columns
defined by `cn`
