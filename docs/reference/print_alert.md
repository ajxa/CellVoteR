# Print a styled alert message

A flexible wrapper for `cli` alert functions that uses `cli`'s theming
system for consistent styling. Supports glue-style interpolation
including `cli` in-line text formatting for the following: `{.val ...}`;
`{.arg ...}`; `{.field ...}` and `{.path ...}`. Moreover, each of these
in-line values is rendered in a distinct highlight color (defaulting to
`"grey98"`), making it easy to draw attention to specific variables or
values within the message. The alert type controls the leading icon,
while the base text color and font face can be customized for further
emphasis.

## Usage

``` r
print_alert(
  text,
  type = "n",
  face = "n",
  color = "grey80",
  highlight_color = "grey98",
  .envir = parent.frame()
)
```

## Arguments

- text:

  Character string. The message to display. Supports `cli` glue-style
  interpolation and inline markup (e.g. `{.val x}`, `{.emph x}`,
  `{?singular/plural}`).

- type:

  Character string. Alert type controlling the leading icon:

  `"s"`

  :   Success (green tick)

  `"d"`

  :   Danger (red cross)

  `"w"`

  :   Warning (orange exclamation)

  `"i"`

  :   Info (blue/cyan 'i')

  `"n"`

  :   Plain text, no icon (default)

  Any other value falls back to a generic alert (arrow).

- face:

  Character string. Font face applied to the base text:

  `"n"`

  :   Normal / plain (default)

  `"b"`

  :   Bold

  `"i"`

  :   Italic

  `"bi"`

  :   Bold italic

- color:

  Character string. Base text color applied to the message body. Accepts
  standard color names (e.g. `"red"`) or hex codes (e.g. `"#FF5733"`).
  Defaults to `"grey80"`.

- highlight_color:

  Character string. Color applied to `{.val ...}`; `{.arg ...}`;
  `{.field ...}` and `{.path ...}` spans. Each of the aforementioned
  spans are rendered using their default themes: only the colour is
  changed according to this argument which accepts the same values as
  `color` and defaults to `"grey98"`.

- .envir:

  Environment used for glue interpolation. Defaults to the caller's
  environment so that variables referenced in `text` are resolved from
  the calling scope.

## Value

Invisible `NULL`; called for its side effect of printing a formatted
message to the console.

## Details

Styling is applied via
[`start_app`](https://cli.r-lib.org/reference/start_app.html) themes
rather than pre-wrapping the text in ANSI codes. This ensures that
`cli`'s own inline markup (`{.val}`, `{.emph}`, etc.) composes correctly
with the user-specified base color and face.

The `"n"` (plain) alert type uses
[`cli_text`](https://cli.r-lib.org/reference/cli_text.html) under the
`body` theme class, producing output with no icon and no leading
whitespace — visually identical to the other alert types minus the
bullet character.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Plain message — no icon, default grey
  print_alert("Starting process...")

  # Highlight interpolated values
  n <- 5
  print_alert("Loaded {.val {n}} sample{?s}", type = "s", color = "green")

  # Warning with italic text and custom highlight
  print_alert("Column {.val {col}} has missing values",
              type = "w", color = "orange", face = "i",
              highlight_color = "yellow")

  # Bold italic danger alert
  print_alert("Fatal: {.val {msg}}", type = "d", color = "red", face = "bi")
} # }
```
