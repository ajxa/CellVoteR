# Format text for the CLI progress steps

A helper function that applies ANSI color and font styles (bold/italic)
to a text string. This is designed to be used inside
[`cli::cli_progress_step()`](https://cli.r-lib.org/reference/cli_progress_step.html)
to style the status message while preserving the spinner's
functionality.

## Usage

``` r
step_fmt(text, color = "grey60", face = "n")
```

## Arguments

- text:

  Character string. The message text to format.

- color:

  Character string. The text color. Accepts standard names (e.g., "red",
  "blue") or hex codes. Defaults to "grey60".

- face:

  Character string. A shorthand code for the font face:

  - `"n"`: Normal (Default)

  - `"b"`: Bold

  - `"i"`: Italic

  - `"bi"`: Bold Italic

## Value

A character string containing the ANSI-styled text, ready to be passed
to CLI printing functions.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Standard usage inside a progress step
  cli::cli_progress_step(step_fmt("Loading data...", color = "blue"))

  # Bold and Red for critical steps
  cli::cli_progress_step(step_fmt("Critical Check", color = "red", face = "b"))

  # Italic Grey (Subtle)
  cli::cli_progress_step(step_fmt("Calculating metrics", face = "i"))
} # }
```
