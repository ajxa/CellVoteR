# Unpack marker list to Tibble

Converts a named list of markers (either simple or nested) into a
long-format tibble.

## Usage

``` r
.unpack_markers(marker_list, val_col = "gene", ind_col = "type")
```

## Arguments

- marker_list:

  A list of markers. Can be a simple named list of vectors or a nested
  list.

- val_col:

  Name of the output column for values (default: "gene").

- ind_col:

  Name of the output column for indices/names (default: "type").

## Value

A tibble with columns defined by `val_col` and `ind_col`.
