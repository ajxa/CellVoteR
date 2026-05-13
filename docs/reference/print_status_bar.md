# Print a status bar with icons

Prints a summary line with two statuses separated by a vertical bar.
Success items get a Green Tick; Failure/Unknown items get a Red Cross.

## Usage

``` r
print_status_bar(
  n_pass,
  n_fail,
  label_text = "Labelled",
  unknown_text = "Unknown",
  use_props = FALSE
)
```

## Arguments

- n_pass:

  Number (or percentage) for the "good" category (e.g., labelled cells).

- n_fail:

  Number (or percentage) for the "bad" category (e.g., unknown cells).

- label_text:

  Character. Label for the good category. Defaults to "labelled".

- unknown_text:

  Character. Label for the bad category. Defaults to "unknown".

- use_props:

  Logical. If TRUE, appends a "%" sign to the numbers. Defaults to
  FALSE.
