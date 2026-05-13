# Define a complete feature-space

Extracts a unique vector of all features present in both the broad and
fine marker lists. Handles both structured (complex) and simple broad
marker lists.

## Usage

``` r
.feature_space(broad_marker_list, fine_marker_list)
```

## Arguments

- broad_marker_list:

  A list of broad markers (simple or structured).

- fine_marker_list:

  A list of fine markers.

## Value

A unique character vector of all features found in the inputs.
