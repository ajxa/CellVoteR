# Check for missing arguments in parent call

Checks if any specific arguments are missing in the function calling
this helper. This looks into the
[`parent.frame()`](https://rdrr.io/r/base/sys.parent.html) to inspect
the arguments of the caller.

## Usage

``` r
.check_missing_args(required_args)
```

## Arguments

- required_args:

  A character vector of argument names to check for presence.

## Value

NULL. Throws an error if any required arguments are missing.
