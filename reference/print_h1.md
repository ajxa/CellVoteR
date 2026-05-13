# Print a styled H1 header with optional contextual information

This function prints a double-line horizontal rule to the console,
acting as a major section header. It supports a main text label and an
optional, italicized "info" label (e.g., for Sample IDs or timestamps),
allowing for distinct coloring of the line, the main text, and the info
text.

## Usage

``` r
print_h1(
  text,
  info = NULL,
  head_color = "cyan",
  text_color = "cyan",
  info_color = "grey60"
)
```

## Arguments

- text:

  Character string. The primary header text.

- info:

  Character string (optional). Secondary text to display in parentheses
  next to the main text. This text is automatically italicized. Defaults
  to `NULL`.

- head_color:

  Character string. Color of the horizontal rule lines. Defaults to
  `"cyan"`.

- text_color:

  Character string. Color of the main text. Defaults to `"white"`.

- info_color:

  Character string. Color of the secondary info text. Defaults to
  `"grey60"`.

## Value

No return value; prints a formatted message to the console.

## Details

The function combines
[`cli::cli_rule`](https://cli.r-lib.org/reference/cli_rule.html) for the
structure with inline ANSI styling for the labels. The `info` argument,
if provided, is formatted as `(info)`, colored according to
`info_color`, and italicized using
[`cli::style_italic`](https://cli.r-lib.org/reference/ansi-styles.html).

## Examples

``` r
if (FALSE) { # \dontrun{
  # Standard usage
  print_h1("Data Ingestion")

  # With context info (e.g. Sample ID)
  # Output: == CellVoterR (Sample_001) ==
  print_h1("CellVoterR", info = "Sample_001")

  # Custom colors
  print_h1("Analysis", info = "Batch A", head_color = "red", info_color = "yellow")
} # }
```
