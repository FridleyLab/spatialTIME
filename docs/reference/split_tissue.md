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

- workers:

  number of cores to use for calculations

- overwrite:

  whether to replace `density_compartment`,
  `refined_density_compartment`, `mif$derived$density_boundary` and
  `mif$sample`'s `Boundary Length` column if any already exist. There is
  no `Run`-based append here (unlike the metric functions): a cell can
  carry only one compartment label, so run `split_tissue()` twice under
  different `sigma`/`interface_width` on two separate mifs if you want
  to compare them.

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
  `c(class1, class2)`) and `refined_density_compartment` (factor, levels
  `c(class1, "Interface", class2)`)

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
40, 1803 cells): the number of contour pieces is 9 at every resolution
from `eps = sigma` down to `eps = sigma/32` – topology is set by
`sigma`, not resolution – while boundary length converges from 6016.8
(`eps = sigma`, -11.5%) to 6799.1 (`eps = sigma/32`, reference).
`eps = sigma/8` (6756.5, -0.6% bias) is used internally as a
resolution/cost tradeoff; there is nothing left for the user to tune.
`Boundary Length` is not scale-free – compare
`Boundary Length / sqrt(area)` across cores of different size.

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
