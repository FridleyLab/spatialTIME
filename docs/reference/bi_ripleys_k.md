# Bivariate Ripley's K

`bi_ripleys_k()` takes a `mIF` object plus marker names and a range of
radii, and measures bivariate clustering (co-localization) between each
ordered pair of markers. Cells positive for both markers of a pair are
excluded from that pair, so the anchor and counted sets are always
disjoint.

Either estimate CSR by permutation (`permute = TRUE`) or use the exact
CSR estimate (`permute = FALSE`). The exact estimate is the univariate K
of *all* cells in the sample, which is the closed form for the expected
cross-K under random labelling – verified against a 500-permutation
Monte Carlo to within Monte Carlo error.

## Usage

``` r
bi_ripleys_k(
  mif,
  mnames,
  r_range = 0:100,
  edge_correction = "translation",
  num_permutations = 50,
  permute = FALSE,
  keep_permutation_distribution = FALSE,
  overwrite = FALSE,
  workers = 1,
  xloc = NULL,
  yloc = NULL,
  big = 10000,
  ...
)
```

## Arguments

- mif:

  mIF object with spatial data frames, clinical, and per-sample summary
  information

- mnames:

  vector of column names for phenotypes, or a two-column data frame of
  specific anchor/counted marker combinations to run

- r_range:

  vector range of radii at which to calculate co-localization *K*

- edge_correction:

  edge correction method: one of "translation", "isotropic", "border" or
  "none"

- num_permutations:

  integer number of permutations used to estimate CSR. Ignored when
  `permute = FALSE`.

- permute:

  whether to estimate CSR by permutation (`TRUE`) or to use the exact
  closed-form CSR estimate (`FALSE`, the default)

- keep_permutation_distribution:

  boolean; keep each permutation's result or average them to one row per
  marker pair and radius

- overwrite:

  boolean; replace an existing `bivariate_Count` slot rather than
  appending it as a new `Run`

- workers:

  integer number of CPU workers to use

- xloc, yloc:

  the x and y columns giving cell centres. If left `NULL`, `XMin`,
  `XMax`, `YMin` and `YMax` must be present.

- big:

  cell count above which per-pair edge weights are computed in chunks to
  bound peak memory. Memory and speed only – results are identical
  either way and the requested `edge_correction` is always honoured.

- ...:

  support for deprecated argument names. `keep_perm_dis` aliases
  `keep_permutation_distribution` and `nlarge` aliases `big`; `force` is
  accepted and ignored (there is no longer a cell-count limit). Anything
  else is an error.

## Value

mif object with bivariate Ripley's K calculated

## Whole slide images

As of 2.0.0 this function handles whole-slide-scale data directly and
[`bi_ripleys_k_WSI()`](https://fridleylab.github.io/spatialTIME/reference/bi_ripleys_k_WSI.md)
is gone. Only the cell pairs closer than `max(r_range)` are ever
materialised, so memory scales with the number of nearby pairs rather
than with n squared – at 200,000 cells roughly 75 MB instead of 319 GB.
The `big` argument bounds peak memory further by chunking, and unlike
the old `big`/`nlarge` arguments it never changes the statistic. In
particular the requested edge correction is no longer silently replaced
with `"none"` above a cell-count threshold.

## Accuracy

Values agree with
[`spatstat.explore::Kcross()`](https://rdrr.io/pkg/spatstat.explore/man/Kcross.html)
to floating-point precision. The observation window is the convex hull
of **every** cell in the sample and is held fixed across all marker
pairs and permutations.

## Examples

``` r
x <- spatialTIME::create_mif(clinical_data = spatialTIME::example_clinical %>%
                               dplyr::mutate(deidentified_id = as.character(deidentified_id)),
                             sample_data = spatialTIME::example_summary %>%
                               dplyr::mutate(deidentified_id = as.character(deidentified_id)),
                             spatial_list = spatialTIME::example_spatial[1],
                             patient_id = "deidentified_id",
                             sample_id = "deidentified_sample")
x2 = bi_ripleys_k(mif = x,
                  mnames = c("CD3..Opal.570..Positive", "CD8..Opal.520..Positive"),
                  r_range = seq(0, 100, 10),
                  edge_correction = "translation",
                  permute = FALSE,
                  workers = 1)
```
