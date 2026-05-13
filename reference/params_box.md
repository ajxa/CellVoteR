# Print a formatted box of parameters

prints a stylized CLI box containing a list of parameters and their
values to the console. It is designed for logging configuration settings
at the start of a script or process.

## Usage

``` r
params_box(param_list)
```

## Arguments

- param_list:

  A named list. The names of the list elements are used as labels, and
  the values are displayed in blue.

## Value

No return value, called for side effects (printing to console).

## Examples

``` r
if (FALSE) { # \dontrun{
  config <- list(
    date = "2023-10-27",
    model = "XGBoost",
    iterations = 100,
    verbose = TRUE
  )
  params_box(config)
} # }
```
