# Resolve final cell labels via majority vote

Combines multiple annotation inputs and two tie-breakers to determine a
single consensus label for each cell.

## Usage

``` r
resolve_consensus_labels(
  label_list,
  method_names,
  tie_breaker_names,
  unassigned_label = "Unknown",
  ordered_tiebreak = TRUE,
  allow_even_split = FALSE
)
```

## Arguments

- label_list:

  Named list of factors or character vectors, one entry per cell. Must
  contain all names in `method_names` and `tie_breaker_names`, and all
  elements must be the same length.

- method_names:

  Character vector. Names within `label_list` representing the primary
  voting methods.

- tie_breaker_names:

  Character vector of length 2. Names within `label_list` to use as
  ordered fallbacks.

- unassigned_label:

  Character scalar. Label assigned when no consensus is reached.
  Defaults to `"Unknown"`.

- ordered_tiebreak:

  Logical scalar. When `TRUE` (default), tie-breakers are consulted in
  order: tie-breaker 1 is tried before tie-breaker 2. When `FALSE`,
  either tie-breaker agreeing with the leading candidate is sufficient,
  with no priority between them.

- allow_even_split:

  Logical scalar. When `FALSE` (default), a majority requires strictly
  more than 50\\ 4). When `TRUE`, a label with more votes than any other
  is accepted as a majority even if it does not exceed 50\\ others have
  1).

## Value

A named list with two equal-length factors:

- `label`:

  The resolved cell type label.

- `method`:

  The decision rule used to reach that label.

## Details

**Decision hierarchy (evaluated per cell):**

1.  **Strong majority:** If any label receives strictly more votes than
    all others, and that count exceeds the majority threshold, it is
    selected immediately. The threshold behaviour is controlled by
    `allow_even_split`.

2.  **Split — tie-breaker agreement:** If there is a genuine leading
    candidate (at least one label has more votes than the others), and
    both tie-breakers agree with each other and match that candidate,
    the label is assigned.

3.  **Split — priority order** (when `ordered_tiebreak = TRUE`, the
    default): Tie-breaker 1 is tried first against the leading
    candidate. If it does not match, tie-breaker 2 is tried.

4.  **Split — either agrees** (when `ordered_tiebreak = FALSE`): Either
    tie-breaker agreeing with the leading candidate is sufficient, with
    no priority between them.

5.  **Unresolved:** If every method disagrees (no leading candidate), or
    no tie-breaker can resolve the split, `unassigned_label` is
    assigned.
