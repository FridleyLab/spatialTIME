# Plot density-based tissue compartments alongside the density difference

[`split_tissue()`](https://fridleylab.github.io/spatialTIME/reference/split_tissue.md)
discards the density images and point patterns it computes the boundary
from, to avoid inflating the mif. `plot_tissue_split()` recomputes them
on demand, at plot time, from the settings recorded by
[`split_tissue()`](https://fridleylab.github.io/spatialTIME/reference/split_tissue.md),
and draws the *stored* boundary polyline (never a recontoured one) on
top of both a raster of the density difference and a scatter of the
per-cell compartment label.

## Usage

``` r
plot_tissue_split(
  mif,
  which = NULL,
  compartment = c("refined_density_compartment", "density_compartment"),
  panels = c("both", "density", "compartment"),
  colors = NULL,
  point_size = 0.4,
  raster_max_pixels = 5e+05,
  workers = 1,
  filename = NULL,
  path = NULL,
  settings = NULL,
  ...
)
```

## Arguments

- mif:

  object of class `mif` that has been through
  [`split_tissue()`](https://fridleylab.github.io/spatialTIME/reference/split_tissue.md)

- which:

  samples to plot: `NULL` for all, or a numeric/character vector
  indexing/naming elements of `mif$spatial`

- compartment:

  which compartment column to colour cells by: the 3-level
  `"refined_density_compartment"` (default) or the 2-level
  `"density_compartment"`

- panels:

  which panel(s) to draw: `"both"` (default, one faceted plot with the
  density-difference raster and the compartment scatter side by side),
  `"density"` alone, or `"compartment"` alone

- colors:

  named character vector of colours for `class1`, `class2` and (for
  `compartment = "refined_density_compartment"`) `"Interface"`. If
  `NULL`, defaults to `RColorBrewer::brewer.pal(3, "Set1")[1:2]` for the
  two classes and `"grey20"` for `"Interface"`; `NA`/unclassified cells
  are always drawn `"grey70"`.

- point_size:

  size passed to
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)
  for the compartment panel

- raster_max_pixels:

  the density raster is recomputed at plot time (see Details) at up to
  this many pixels; above the cap resolution is capped for display only
  – the drawn boundary line is unaffected and stays exact.

- workers:

  number of cores to use, one sample per core

- filename, path:

  if `filename` is given, every requested plot is also written to a
  single PDF at `file.path(path, filename)` (`path` defaults to the
  working directory; a trailing `.pdf` in `filename` is not doubled).

- settings:

  provenance normally read from
  `attr(mif$derived$density_boundary, "call_info")`, which
  [`split_tissue()`](https://fridleylab.github.io/spatialTIME/reference/split_tissue.md)
  attaches. That attribute does not survive `[` subsetting or
  [`dplyr::bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html)
  (e.g. after
  [`merge_mifs()`](https://fridleylab.github.io/spatialTIME/reference/merge_mifs.md)
  disagrees across mifs), so pass the settings list explicitly to
  recover from that: a list with `classifier`, `class1`, `class2`,
  `sigma`, `eps`, `interface_width`, `xloc`, `yloc`, `filter_density`
  and `sample_id`.

- ...:

  accepts no arguments; present only so that a mistyped named argument
  above produces an informative "Unknown argument" error instead of
  being silently absorbed.

## Value

a named list of `ggplot` objects, one per requested sample, named like
`mif$spatial` – **not** the `mif`. Unlike
[`plot_immunoflo()`](https://fridleylab.github.io/spatialTIME/reference/plot_immunoflo.md),
this deliberately does not attach the plots to `mif$derived`: each plot
captures a raster data frame that can be as large as the sample itself,
and
[`merge_mifs()`](https://fridleylab.github.io/spatialTIME/reference/merge_mifs.md)'s
handling of list-valued derived slots is fragile enough already (see
`NEWS.md`) without doubling the exposure.

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
plots <- plot_tissue_split(x)
```
