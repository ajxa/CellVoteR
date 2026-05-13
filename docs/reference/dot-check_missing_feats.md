# Check for missing features in a set

Verifies that all features in `set1` are present in `set2`.

## Usage

``` r
.check_missing_feats(set1, set2)
```

## Arguments

- set1:

  A character vector of query features (e.g., markers to look for).

- set2:

  A character vector of target features (e.g., available genes in data).

## Value

NULL. Throws an error if any features in `set1` are missing from `set2`.
