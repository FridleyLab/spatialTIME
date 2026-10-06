# Split a sample into tissue compartments by density difference

For each spatial sample, `split_tissue()` computes kernel density
estimates of `class1` and `class2`, takes their difference, and extracts
the *exact* zero level set of that difference with
[`grDevices::contourLines()`](https://rdrr.io/r/grDevices/contourLines.html)
– not a thresholded band around zero. Every cell (not just
`class1`/`class2` cells) is then labelled by which side of that boundary
it falls on.

## Usage

``` r
split_tissue(
  mif,
  classifier,
  class1,
  class2,
  sigma,
  interface_width,
  rescale = TRUE,
  min_density = NULL,
  hard_threshold = FALSE,
  workers = 1,
  overwrite = FALSE,
  xloc = NULL,
  yloc = NULL,
  ...
)
```

## Arguments

- mif:

  object of class `mif` created with
  [`create_mif()`](https://fridleylab.github.io/spatialTIME/reference/create_mif.md)

- classifier:

  column in each spatial file giving the per-cell tissue classification
  (e.g. Tumor/Stroma/Lymph/Necrotic). May have more than two levels;
  only `class1` and `class2` are compared.

- class1, class2:

  the two classifier levels to compare. The density difference is
  `class1` minus `class2`, so a cell on the positive side of the
  boundary is labelled `class1` and a cell on the negative side is
  labelled `class2`. Neither may be named `"Interface"` (reserved) and
  they must differ.

- sigma:

  kernel density bandwidth, in the same coordinate units as the spatial
  data. There is no default: the boundary's shape depends on it and
  values that make sense differ by imaging platform, so a silent default
  would make cores or cohorts processed with different (undocumented)
  defaults incomparable. Memory is quadratic in `1/sigma` – see
  `check_density_budget()` in `R/utils-density-boundary.R` and the
  pixel-budget error below.

- interface_width:

  full width (diameter, not radius) of the interface band, in the same
  units as `sigma`. Cells within `interface_width / 2` of the boundary
  are labelled `"Interface"` in `refined_density_compartment`.

- rescale:

  rescale each class's density image to `[0, 1]` before differencing
  (default `TRUE`). The boundary then sits where the two *rescaled*
  densities are equal, which is not the same curve as where the absolute
  intensities are equal: each image gets its own affine map, so the zero
  set moves rather than being reparameterised. On
  `example_spatial[["TMA3_[9,K].tif"]]` at `sigma = 40` the unrescaled
  field gives 9 contour pieces / length 6756.5 and the rescaled one 7 /
  6723.8. The practical consequence is that under `rescale = TRUE` the
  boundary – and so `Boundary Length` – depends on each sample's own
  density range, so it is not an absolute criterion you can compare
  across a cohort on its own terms. Use `rescale = FALSE` for that. The
  upside is that a sparse class is not swamped by an abundant one, which
  matters when compartment proportions vary a lot between samples.

- min_density:

  drop pixels where there is essentially no tissue, given as a fraction
  of the sample's own mean `class1 + class2` intensity; `NULL` (default)
  disables it. Where both classes are near zero the difference is at the
  floating-point noise floor and its sign is meaningless, so the zero
  contour fragments into noise there. Measured on a whole-slide sample
  with a large tissue hole: unmasked, the boundary had 227 pieces and
  74% of its length lay in space with essentially no cells; at
  `min_density = 0.05` it had 23 pieces, and only 5 of 1,000,977 cells
  changed compartment. The threshold is insensitive – anything from 0.01
  to 0.10 gave the same answer on that sample, because the gap it
  straddles is several orders of magnitude – but above roughly 0.2 it
  starts clipping real boundary, so check how many cells change label if
  you raise it. Applied to the raw intensities, before `rescale`.

- hard_threshold:

  collapse the density difference to its sign before contouring (default
  `FALSE`). Since the zero level set is already invariant to any
  monotone transform, this cannot change the contour's topology – and
  measurably does not: piece counts are unchanged. What it does change
  is that
  [`contourLines()`](https://rdrr.io/r/grDevices/contourLines.html) can
  then only place a crossing at the midpoint between two pixel centres,
  quantising the boundary to the grid and adding a staircase that
  inflates `Boundary Length` by 5-7%. Provided for comparison with
  implementations that threshold first; leave it off unless you need to
  reproduce one.

- workers:

  number of cores to use for calculations

- overwrite:

  whether to replace `density_compartment`,
  `refined_density_compartment`, `density_score`,
  `mif$derived$density_boundary` and `mif$sample`'s `Boundary Length`
  column if any already exist. There is no `Run`-based append here
  (unlike the metric functions): a cell can carry only one compartment
  label, so run `split_tissue()` twice under different
  `sigma`/`interface_width` on two separate mifs if you want to compare
  them.

- xloc, yloc:

  columns giving the cell centre. If left `NULL`, `XMin`, `XMax`, `YMin`
  and `YMax` must be present and the centre is their midpoint.

- ...:

  `filter_density`, an optional `function(im) im` applied to each
  class's density image before differencing, to suppress the boundary
  inflating/bouncing around near-zero-density holes in the tissue (a
  hole is close to zero in *both* compartments, not just one). It
  affects only the boundary geometry (contour, length, interface
  distances) – the sign that drives
  `density_compartment`/`refined_density_compartment` always comes from
  the *unfiltered* difference, so a filtered-out hole never leaves a
  cell `NA`. Must return an `im` on the same pixel grid it was given.
  Anything else in `...` is an error naming the offender.

## Value

object of class `mif`, with:

- spatial:

  each sample gains `density_compartment` (factor, levels
  `c(class1, class2)`), `refined_density_compartment` (factor, levels
  `c(class1, "Interface", class2)`) and `density_score`, the signed
  `class1 - class2` difference at the cell's own location.
  `density_score` is the field the two factors are derived from, so its
  sign always agrees with `density_compartment`; keeping the magnitude
  gives you "how far into this compartment" as a covariate rather than
  only "which side". Its units follow `rescale` – roughly `[-1, 1]` when
  `TRUE`, intensity difference when `FALSE` – so it is comparable across
  samples only when `rescale = FALSE`.

- derived\$density_boundary:

  a named list (by `names(mif$spatial)`), one data frame per sample with
  columns `<sample_id>`, `piece`, `x`, `y` (0 rows if no boundary
  exists), carrying a `"call_info"` attribute recording the settings
  used, consumed internally by
  [`plot_tissue_split()`](https://fridleylab.github.io/spatialTIME/reference/plot_tissue_split.md)

- sample:

  gains a `Boundary Length` column (`NA` for samples with no matching
  spatial frame; `0`, not `NA`, for a sample with data but no contour)

## What is and is not kept

The density images and point patterns used to compute the boundary are
discarded once the boundary and per-cell labels are derived – keeping
them would multiply the size of the mif by the number of samples. Only
the boundary polyline (`mif$derived$density_boundary`) and the two new
spatial columns survive. To visualise the density alongside the
boundary, call
[`plot_tissue_split()`](https://fridleylab.github.io/spatialTIME/reference/plot_tissue_split.md),
which recomputes the density on demand at plot time.

## Choosing sigma

`sigma` has no default. Its units are the same as your spatial
coordinates, which differ by imaging platform, so there is no value that
is reasonable for every dataset. Memory is quadratic in `1/sigma` (pixel
size is `sigma / 8` internally, so halving `sigma` quadruples the
density grid) – too small a `sigma` errors naming a minimum viable value
for that sample rather than exhausting memory.

## Boundary length

Pixel resolution (`dimyx` in the underlying
[`spatstat.explore::density.ppp()`](https://rdrr.io/pkg/spatstat.explore/man/density.ppp.html)
call) is not exposed. Measured on a real core (Tumor vs Stroma, sigma
40, 1803 cells, `rescale = FALSE`): the number of contour pieces is 9 at
every resolution from `eps = sigma` down to `eps = sigma/32` – topology
is set by `sigma`, not resolution – while boundary length converges from
6016.8 (`eps = sigma`, -11.5%) to 6799.1 (`eps = sigma/32`, reference).
`eps = sigma/8` (6756.5, -0.6% bias) is used internally as a
resolution/cost tradeoff; there is nothing left for the user to tune.
Those figures are from the unrescaled field, but the argument they
support is about the `eps/sigma` ratio, which `rescale` does not affect.

`Boundary Length` is not scale-free – compare
`Boundary Length / sqrt(area)` across cores of different size – and
under the default `rescale = TRUE` it is additionally relative to each
sample's own density range. Two things worth checking before using it as
a covariate: that it is not dominated by contour in near-empty tissue
(see `min_density`), and that it is not inflated by grid quantisation
(see `hard_threshold`).

## Examples

``` r
library(dplyr)
x <- create_mif(clinical_data = example_clinical %>%
  mutate(deidentified_id = as.character(deidentified_id)),
  sample_data = example_summary %>%
  mutate(deidentified_id = as.character(deidentified_id)),
  spatial_list = example_spatial[1],
  patient_id = "deidentified_id",
  sample_id = "deidentified_sample")
x <- split_tissue(x, classifier = "Classifier.Label",
  class1 = "Tumor", class2 = "Stroma",
  sigma = 40, interface_width = 100)
```
