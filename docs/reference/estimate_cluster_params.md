# Estimate clustering parameters using the cell count

Computes sensible defaults for the number of principal components, SNN
graph neighbourhood size (K), and Leiden clustering resolution based on
the number of cells. All three parameters scale dynamically with dataset
size using bounded transformations and are intended to allow for slight
over-clustering.

## Usage

``` r
estimate_cluster_params(
  n_cells,
  min_cluster_cells = 50L,
  res_min = 0.6,
  res_max = 2,
  res_saturation_n = 100000L,
  k_min = 10L,
  k_max = 50L,
  k_saturation_n = 200000L,
  npc_min = 20L,
  npc_max = 50L,
  npc_slope = 4
)
```

## Arguments

- n_cells:

  Integer scalar. Number of cells in the dataset.

- min_cluster_cells:

  Integer scalar. Minimum number of cells required to attempt
  clustering. If `n_cells` is below this threshold the function returns
  a parameter set with `skip = TRUE`. Defaults to `50`.

- res_min, res_max:

  Numeric scalars. Bounds for Leiden resolution. Defaults to `0.6` and
  `2.0`.

- res_saturation_n:

  Integer scalar. Cell count at which resolution reaches `res_max`.
  Defaults to `100000`.

- k_min, k_max:

  Integer scalars. Bounds for SNN nearest-neighbour K. Defaults to `10`
  and `50`.

- k_saturation_n:

  Integer scalar. Cell count at which K reaches `k_max`. Defaults to
  `200000`.

- npc_min, npc_max:

  Integer scalars. Bounds for number of principal components. Defaults
  to `20` and `50`.

- npc_slope:

  Numeric scalar. Controls the log2-based scaling rate for PCs. Defaults
  to `4`.

## Value

A named list with components:

- n_pcs:

  Integer. Number of principal components.

- dims:

  Integer vector `seq_len(n_pcs)`.

- k:

  Integer. Nearest-neighbour count for SNN graph.

- resolution:

  Numeric. Leiden clustering resolution.

- skip:

  Logical. `TRUE` if `n_cells` is below the minimum threshold.

## Parameter scaling logic

- **PCs**:

  Log2-scaled between `npc_min` and `npc_max`. Small datasets get fewer
  components to avoid overfitting noise; large datasets saturate at
  `npc_max`.

- **K (SNN neighbours)**:

  Square-root-scaled between `k_min` and `k_max`. Grows slowly because
  neighbourhood size affects graph topology - too large and distinct
  populations merge.

- **Resolution (Leiden)**:

  Square-root-scaled between `res_min` and `res_max`. Grows faster than
  K because larger datasets typically contain more distinct communities
  and require finer partitioning.

## Examples

``` r
if (FALSE) { # \dontrun{
params <- estimate_cluster_params(n_cells = ncol(sce))
if (!params$skip) {
  ctx <- prepare_context(sce, n_pcs = params$n_pcs, k = params$k)
}
} # }
```
