# Univariate Nearest Neighbor G(r)

`NN_G()` computes the nearest-neighbour distance distribution function
G(r) for each marker: the proportion of marker-positive cells whose
nearest marker-positive neighbour lies within r. CSR is estimated by
permutation – relabelling the same cell locations – so the null accounts
for the sample's own geometry rather than assuming a homogeneous Poisson
process.

## Usage

``` r
NN_G(
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
  whether a cell is positive for a phenotype

- r_range:

  numeric vector of radii at which to evaluate G(r)

- num_permutations:

  integer number of permutations used to estimate sample-specific
  complete spatial randomness (CSR)

- edge_correction:

  edge correction method: one of "rs", "km", "han" or "none"

- keep_permutation_distribution:

  boolean; keep each permutation's result or average them to one row per
  marker and radius

- workers:

  integer number of CPU cores used to process samples in parallel

- overwrite:

  boolean; replace an existing `univariate_NN` slot rather than
  appending it as a new `Run`

- xloc, yloc:

  the x and y columns giving cell centres. If left `NULL`, `XMin`,
  `XMax`, `YMin` and `YMax` must be present.

- ...:

  support for deprecated argument names. `keep_perm_dis` is accepted as
  an alias for `keep_permutation_distribution`. Anything else is an
  error.

## Value

object of class `mif` with a `univariate_NN` table in the `derived` slot

## Accuracy

Estimates come from
[`spatstat.explore::Gest()`](https://rdrr.io/pkg/spatstat.explore/man/Gest.html),
so they agree with spatstat exactly. The observation window is the
convex hull of **every** cell in the sample and is held fixed across all
markers and permutations, which is what makes G comparable between
markers within a sample.

## Why there is no exact CSR for G

[`ripleys_k()`](https://fridleylab.github.io/spatialTIME/reference/ripleys_k.md)
can skip permutations because the K of all cells is the closed form for
the expected K of a random subset of them. No such shortcut exists for
G, because G depends on the *intensity* of the point set, not just its
geometry: a random subset of the cells has a lower intensity and
therefore larger nearest-neighbour distances. Concretely, on one of the
example samples `Gest()` over all cells gives 0.52 at a radius where the
mean permuted G is 0.13. The `Exact CSR` column is therefore present for
schema consistency with the other metrics but is always `NA`; use
`Permuted CSR` or `Theoretical CSR`.

## Examples

``` r
x <- spatialTIME::create_mif(clinical_data = spatialTIME::example_clinical %>%
  dplyr::mutate(deidentified_id = as.character(deidentified_id)),
  sample_data = spatialTIME::example_summary %>%
  dplyr::mutate(deidentified_id = as.character(deidentified_id)),
  spatial_list = spatialTIME::example_spatial[1],
  patient_id = "deidentified_id",
  sample_id = "deidentified_sample")

x2 = NN_G(mif = x, mnames = "CD3..Opal.570..Positive",
  r_range = seq(0, 100, 10), num_permutations = 10,
  edge_correction = "rs", workers = 1, overwrite = TRUE)
```
