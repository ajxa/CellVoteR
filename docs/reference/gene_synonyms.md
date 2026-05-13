# Retrieve gene synonyms from Ensembl

Connects to the Ensembl BioMart database to retrieve gene synonyms for a
provided list of gene names. It queries both the 'external_gene_name'
and 'external_synonym' filters to maximize match rates.

## Usage

``` r
gene_synonyms(lookup_genes)
```

## Arguments

- lookup_genes:

  Character vector. A list of gene names (or synonyms) to search for.

## Value

A named list where:

- **Names** are the standard Ensembl external gene names.

- **Values** are character vectors of unique synonyms for that gene.

Returns `NULL` if the connection to Ensembl fails or if no genes are
found.

## Details

The function performs the following steps:

1.  Checks if `biomaRt` is installed.

2.  safely connects to the `hsapiens_gene_ensembl` dataset.

3.  Queries Ensembl twice: once matching against gene names and once
    matching against existing synonyms.

4.  Cleans the results by converting empty strings `""` to `NA` and
    removing incomplete records.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Basic usage
  syns <- gene_synonyms(c("TP53", "BRCA1"))

  # Access synonyms for TP53
  print(syns$TP53)
} # }
```
