# Extract and save Seurat counts and metadata to RDS

Extract and save Seurat counts and metadata to RDS

## Usage

``` r
save_seurat_components(
  seurat_obj,
  metadata_cols = NULL,
  out_dir = getwd(),
  base_name = "data",
  assay = SeuratObject::DefaultAssay(seurat_obj),
  layer = "counts"
)
```

## Arguments

- seurat_obj:

  A Seurat object.

- metadata_cols:

  A character vector of column names expected in metadata. If NULL
  (default), all metadata columns are extracted.

- out_dir:

  Character. Directory where the files will be saved.

- base_name:

  Character. Prefix for the saved filenames.

- assay:

  Character. The specific assay to get data from or set data for;
  defaults to
  [`DefaultAssay`](https://satijalab.github.io/seurat-object/reference/DefaultAssay.html)

- layer:

  The name of the layer to get (defaults to "count" for raw counts).

## Value

Logical TRUE if successful.
