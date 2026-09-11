# Bivariate Nearest Neighbor G(r)

`bi_NN_G()` computes the cross-type nearest-neighbour distribution
function: for each ordered pair of markers, the proportion of
anchor-positive cells whose nearest *counted*-positive neighbour lies
within r. Cells positive for both markers of a pair are excluded from
that pair, so the anchor and counted sets are always disjoint.

## Usage

``` r
bi_NN_G(
  mif,
  mnames,
  r_range = 0:100,
  num_permutations = 50,
  edge_correction = "rs",
  keep_permutation_distribution = FALSE,
  workers = 1,
  overwrite = FALSE,
  xloc = NULL,
  yloc = NULL,
  ...
)
```

## Arguments

- mif:

  object of class `mif` created by function
  [`create_mif()`](https://fridleylab.github.io/spatialTIME/reference/create_mif.md)

- mnames:

  character vector of column names within the spatial files indicating
  whether a cell is positive for a phenotype, or a two-column data frame
  of specific anchor/counted combinations

- r_range:

  numeric vector of radii at which to evaluate G(r)

- num_permutations:

  integer number of permutations used to estimate sample-specific
  complete spatial randomness (CSR)

- edge_correction:

  edge correction method: one of "rs", "km", "han" or "none"

- keep_permutation_distribution:

  boolean; keep each permutation's result or average them to one row per
  marker pair and radius

- workers:

  integer number of CPU cores used to process samples in parallel

- overwrite:

  boolean; replace an existing `bivariate_NN` slot rather than appending
  it as a new `Run`

- xloc, yloc:

  the x and y columns giving cell centres. If left `NULL`, `XMin`,
  `XMax`, `YMin` and `YMax` must be present.

- ...:

  support for deprecated argument names. `keep_perm_dis` is accepted as
  an alias for `keep_permutation_distribution`. Anything else is an
  error.

## Value

object of class `mif` with a `bivariate_NN` table in the `derived` slot

## Accuracy

Estimates come from
[`spatstat.explore::Gcross()`](https://rdrr.io/pkg/spatstat.explore/man/Gcross.html).
Before 2.0.0 this function hand-rolled the `rs` and `han` estimators on
a full `as.matrix(dist(...))`, needing memory proportional to the square
of the cell count. Those hand-rolled results agreed with `Gcross()`
exactly on well-populated markers, but returned `NaN` where `Gcross()`
correctly returns `0` for sparse markers, and the `rs` branch referenced
an undefined variable. Delegating fixes the sparse-marker case and drops
memory to O(n).

The observation window is the convex hull of **every** cell in the
sample and is held fixed across all marker pairs and permutations.

## Why there is no exact CSR for G

See the corresponding section of
[`NN_G()`](https://fridleylab.github.io/spatialTIME/reference/NN_G.md).
The `Exact CSR` column exists for schema consistency with
[`ripleys_k()`](https://fridleylab.github.io/spatialTIME/reference/ripleys_k.md)
and is always `NA`.

## Examples

``` r
x <- spatialTIME::create_mif(clinical_data = spatialTIME::example_clinical %>%
  dplyr::mutate(deidentified_id = as.character(deidentified_id)),
  sample_data = spatialTIME::example_summary %>%
  dplyr::mutate(deidentified_id = as.character(deidentified_id)),
  spatial_list = spatialTIME::example_spatial[1],
  patient_id = "deidentified_id",
  sample_id = "deidentified_sample")

x2 = bi_NN_G(mif = x,
      mnames = c("CD3..Opal.570..Positive", "CD8..Opal.520..Positive"),
      r_range = seq(0, 100, 10), num_permutations = 10,
      edge_correction = "rs", workers = 1, overwrite = TRUE)
```
