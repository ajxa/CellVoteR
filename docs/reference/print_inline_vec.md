# Print an inline vector summary

Prints a list of labeled counts separated by vertical pipes. Format:
"label 1 (count 1) \| label 2 (count 2) \| ..."

## Usage

``` r
print_inline_vec(counts, labels, text_color = "grey50", num_color = "grey90")
```

## Arguments

- counts:

  Numeric vector. The count values.

- labels:

  Character vector. The names corresponding to the counts.

- text_color:

  Character. Color for the label text.

- num_color:

  Character. Color for the count numbers inside the brackets.
