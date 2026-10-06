# Calculate Ripley's K

`ripleys_k()` calculates the empirical Ripley's K for the cell types
given in `mnames`. This is useful for exploring the spatial clustering
of single cell types on TMA cores, ROI spots, or whole slide images
following phenotyping with a program such as HALO.

Either estimate CSR by permutation (`permute = TRUE`) or use the exact
CSR estimate (`permute = FALSE`). The exact estimate is the K of *all*
cells in the sample, which is the closed form for the expected K of a
randomly chosen subset of them, so it is both faster and free of Monte
Carlo error. Permutations are still useful if you want the full null
distribution rather than its mean – run 1000 and treat an observed value
outside the 95th percentile as significant.

`Exact CSR` is filled in either way, so with `permute = TRUE` you can
read `Permuted CSR` against it and see directly whether your permutation
count was large enough to converge. They agree exactly in expectation.

## Usage

``` r
ripleys_k(
  mif,
  mnames,
  r_range = seq(0, 100, 1),
  num_permutations = 50,
  edge_correction = "translation",
  permute = FALSE,
  keep_permutation_distribution = FALSE,
  workers = 1,
  overwrite = FALSE,
  xloc = NULL,
  yloc = NULL,
  big = 10000,
  ...
)
```

## Arguments

- mif:

  object of class `mif` created with `create_mif`

- mnames:

  cell phenotype markers to calculate Ripley's K for

- r_range:

  radius range (including 0)

- num_permutations:

  number of permutations to use to estimate CSR. Ignored when
  `permute = FALSE`.

- edge_correction:

  edge correction method: one of "translation", "isotropic", "border" or
  "none". Unlike previous versions this is never silently downgraded for
  large samples. All four agree with
  [`spatstat.explore::Kest()`](https://rdrr.io/pkg/spatstat.explore/man/Kest.html)
  to floating-point precision; note that `"border"` reports no
  `Exact CSR` (see below) and that `"none"` applies no correction at all
  and so is biased downward near the window edge.

- permute:

  whether to estimate CSR by permutation (`TRUE`) or to use the exact
  closed-form CSR estimate (`FALSE`, the default and much faster)

- keep_permutation_distribution:

  whether to keep each permutation's result or average them into a
  single row per marker and radius

- workers:

  number of cores to use for calculations

- overwrite:

  whether to overwrite the `univariate_Count` slot within `mif$derived`

- xloc, yloc:

  columns giving the cell centre. If left `NULL`, `XMin`, `XMax`, `YMin`
  and `YMax` must be present and the centre is their midpoint.

- big:

  cell count above which the per-pair edge-correction weights are
  computed in chunks to bound peak memory. This affects memory and speed
  only: results are identical either way, and the requested
  `edge_correction` is always honoured.

- ...:

  support for deprecated argument names. `keep_perm_dis` is accepted as
  an alias for `keep_permutation_distribution`; `method` is accepted and
  ignored (it was never used). Anything else is an error.

## Value

object of class `mif`

## The CSR columns, and which to trust

Three references are reported for each radius, and
`Degree of Clustering X` is always `Observed K - X CSR`:

- `Theoretical CSR` is \\\pi r^2\\, the K of a homogeneous Poisson
  process. It ignores the shape of the window and the cell density
  actually present, so it is the weakest of the three on real tissue.

- `Exact CSR` is the K of every cell in the sample. Under random
  labelling of a fixed set of cell locations this is *exactly* the
  expected K of a marker-positive subset – not an approximation – for
  translation, isotropic and none. It is the reference to prefer.

- `Permuted CSR` is the mean over `num_permutations` random
  relabellings.

`Exact CSR` is `NA` for `edge_correction = "border"`. Border is a
reduced-sample estimator whose denominator counts only the cells still
further than `r` from the window edge, and that count changes with which
cells are marker-positive. So the cancellation that makes the other
three exact does not apply, and the K of all cells is *not* the expected
K of a subset: measured on 600 cells with 120 positive, it sits about 1%
below the mean permuted K at larger radii. Rather than report a number
that is quietly wrong in a column called "Exact", it is omitted – as it
is for
[`NN_G()`](https://fridleylab.github.io/spatialTIME/reference/NN_G.md)
and
[`pair_correlation()`](https://fridleylab.github.io/spatialTIME/reference/pair_correlation.md),
for the same reason. Use `permute = TRUE` with border.

## Permutation p-values

`Permutations Larger than Observed` counts the permutations at least as
extreme as the observation (`>=`), so ties count toward the count and
therefore toward a larger p-value. `Permutation p-value` is the standard
Monte Carlo p-value, \\(1 + \\\\perm \ge obs\\) / (B + 1)\\, where `B`
is the number of permutations that actually produced a value at that
radius. The `+1`s mean it can never be exactly 0, which a raw count
divided by `B` can be, and which would read as infinite significance
rather than "nothing in this sample was more extreme".

Both are `NA` at radii where `Observed K` is `NA`.

## Accuracy

Values agree with
[`spatstat.explore::Kest()`](https://rdrr.io/pkg/spatstat.explore/man/Kest.html)
to floating-point precision. The observation window is the convex hull
of **every** cell in the sample, and it is held fixed across all markers
and permutations, so K values for different markers within a sample are
directly comparable.

Large samples are handled by only ever materialising the cell pairs
closer than `max(r_range)`, rather than a full n-by-n distance matrix.
At 200,000 cells that is the difference between roughly 75 MB and 319
GB. Because of this, `big` no longer changes the statistic – in versions
before 2.0.0 exceeding it silently replaced your `edge_correction` with
`"none"`.

## Examples

``` r
x <- spatialTIME::create_mif(clinical_data = spatialTIME::example_clinical %>%
  dplyr::mutate(deidentified_id = as.character(deidentified_id)),
  sample_data = spatialTIME::example_summary %>%
  dplyr::mutate(deidentified_id = as.character(deidentified_id)),
  spatial_list = spatialTIME::example_spatial[1],
  patient_id = "deidentified_id",
  sample_id = "deidentified_sample")
x2 = ripleys_k(mif = x,
  mnames = "CD3..Opal.570..Positive",
  r_range = seq(0, 100, 10),
  edge_correction = "translation",
  permute = FALSE,
  workers = 1,
  overwrite = TRUE)
```
