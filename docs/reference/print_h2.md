# Print a styled H2 header

Prints a single-line horizontal rule with a split layout. The main text
appears on the left (bold, colored), and optional info appears on the
right (italic, grey).

## Usage

``` r
print_h2(
  text,
  info = NULL,
  line_color = "cyan",
  text_color = "cyan",
  info_color = "grey80"
)
```

## Arguments

- text:

  Character string. The text for the left side of the header.

- info:

  Character string (Optional). Text for the right side, displayed in
  brackets.

- line_color:

  Character string. Color of the horizontal line. Defaults to "cyan".

- text_color:

  Character string. Color of the left-side text. Defaults to "cyan".

- info_color:

  Character string. Color of the right-side info text. Defaults to
  "grey60".
