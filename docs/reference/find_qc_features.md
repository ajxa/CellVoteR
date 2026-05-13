# Auto-detect QC feature groups

Scans the row names of a
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object to identify genes matching specific patterns (e.g. mitochondrial
or ribosomal) and formats them for use in downstream quality control
filtering.

## Usage

``` r
find_qc_features(
  sce,
  group_configs = list(mito = list(pattern = "^MT-", max_pct = 20), ribo = list(pattern =
    "^RP[SL]", max_pct = 50))
)
```

## Arguments

- sce:

  A
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  object.

- group_configs:

  A named list of configurations. Each element must be a list
  containing:

  pattern

  :   A regex string to match gene names (e.g. `"^MT-"`).

  max_pct

  :   Numeric scalar. The maximum allowed percentage for this group,
      carried forward for use in QC filtering.

  Defaults to mitochondrial (`^MT-`, 20%) and ribosomal (`^RP[SL]`, 50%)
  genes.

## Value

A named list where each element contains:

- features:

  Character vector of matching gene names found in the object. Empty
  character vector if no genes matched.

- max_pct:

  Numeric threshold carried over from the config.

## Examples

``` r
if (FALSE) { # \dontrun{
# Default usage
qc_feats <- find_qc_features(sce)

# Custom patterns (e.g. hemoglobin genes)
qc_feats <- find_qc_features(sce, group_configs = list(
  hb = list(pattern = "^HB[AB]", max_pct = 5)
))
} # }
```
