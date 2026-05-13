# Validate broad marker structure

Checks if a broad marker list is named and contains the required
structural fields.

## Usage

``` r
.valid_broad_marker_struct(
  broad_marker_list,
  struct_fields = c("markers", "expr_threshold", "coexp_min", "priority")
)
```

## Arguments

- broad_marker_list:

  A list of broad markers.

- struct_fields:

  A character vector of fields required inside each list element.
  Defaults to `c("markers", "expr_threshold", "coexp_min", "priority")`.

## Value

`TRUE` if the structure is valid, `FALSE` otherwise. Throws an error if
elements are unnamed.
