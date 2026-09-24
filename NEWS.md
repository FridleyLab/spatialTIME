# spatialTIME 2.0.0

A cleanup and correctness release, plus one new capability. Roughly 1,100 lines of
unreachable code are gone, the count-based measures now agree with **spatstat** to
floating-point precision on samples of any size, every metric returns the same
columns, density-based tissue segmentation is new, and test coverage went from 4
assertions to 732.

## New: density-based tissue segmentation

* **`split_tissue()`** segments a sample into tissue compartments from the
  difference of two classes' kernel density estimates (e.g. Tumor vs Stroma).
  Every cell gains `density_compartment` (2 levels, `class1`/`class2`, by the sign
  of the `class1 - class2` density), `refined_density_compartment` (those two plus
  `"Interface"` for cells within `interface_width / 2` of the boundary) and
  `density_score`, the signed difference itself — the field the two factors are
  derived from, kept because "how far into this compartment" is more useful as a
  covariate than "which side". `mif$sample` gains a `Boundary Length` column.
  * `sigma` (the KDE bandwidth) has **no default** — units differ by imaging
    platform, so a silent default would make cores processed with different
    undocumented defaults incomparable.
  * The boundary is the *exact* zero level set of the density difference,
    extracted with `grDevices::contourLines()`, not a thresholded band around
    zero — so there is no `boundary_threshold` to tune.
  * Pixel resolution (`dimyx`) is not exposed either. It is derived as
    `sigma / 8`: contour *topology* is set by `sigma` at every resolution tested
    (9 pieces from `eps = sigma` down to `eps = sigma/32` on a real core, with
    `rescale = FALSE`), while boundary length converges by `sigma/8` (−0.6%
    bias, against −11.5% at `eps = sigma`). Hiding it also avoids the non-square
    pixels `dimyx` gives on a non-square window (4x anisotropy in testing).
  * `rescale` (default `TRUE`) rescales each class's density image to `[0, 1]`
    before differencing, so a sparse class is not swamped by an abundant one. The
    cost: the boundary then sits where the *rescaled* densities are equal, which
    depends on each sample's own density range, so `Boundary Length` is not an
    absolute criterion comparable across a cohort. `rescale = FALSE` gives the
    absolute one — the boundary where the intensities themselves are equal. On
    `TMA3_[9,K].tif` at `sigma = 40`: 7 pieces / length 6723.8 rescaled, 9 /
    6756.5 not.
  * `min_density` (default off) drops pixels where both classes are near zero.
    There the difference sits at the floating-point noise floor and its sign is
    meaningless, so the contour fragments into noise. On a whole slide with a
    large tissue hole, unmasked: 227 contour pieces, with 74% of the boundary
    length lying in space containing essentially no cells (median total density
    1.46e-12 there against 1.24e-3 elsewhere). At `min_density = 0.05`: 23
    pieces, and 5 of 1,000,977 cells changed compartment. Given as a fraction of
    the sample's own mean intensity, so the same value means the same thing across
    samples and coordinate units, and insensitive — anything from 0.01 to 0.10
    gave the same answer, because the gap it straddles is orders of magnitude.
    Above roughly 0.2 it starts clipping real boundary, so count how many labels
    move if you raise it.
  * `filter_density`, through `...`, is an arbitrary `function(im) im` applied to
    each class image, for anything the above does not cover.
  * `min_density` and `filter_density` affect the **boundary geometry only**. The
    per-cell sign and `density_score` always come from the unmasked, unfiltered
    difference, so masking — even masking an entire sample — cannot orphan a cell
    to `NA`.
  * The field is read at each cell by **bilinear interpolation** between the four
    surrounding pixel centres, with a nearest-valid-pixel fallback for cells whose
    neighbourhood is incomplete (23–39 per example core). `contourLines()` places
    the boundary where the linearly interpolated field crosses zero, so
    interpolating the same way makes `density_score` consistent with the drawn
    polyline: evaluated at the polyline's own vertices the field is zero to machine
    precision (`<= 4.3e-15`), where nearest-pixel lookup is off by up to 0.041,
    roughly 4% of the rescaled field's range. The gain is in the *score*, not in
    relabelling — it shifts `density_score` by up to 0.3 for cells near the
    boundary, but flips the 2-level compartment for only 6–12 cells per example
    core, and since essentially all of those lie within half a pixel of the contour
    (and so inside any sensible interface band) the 3-level label changes for 0–1
    cells per core. The fallback also means a lookup never yields `NA`, so an `NA`
    compartment can only mean an exactly-zero field.
  * `hard_threshold` (default off) collapses the field to its sign before
    contouring, for comparison with implementations that threshold first. It
    cannot change the contour's topology, since the zero level set is invariant to
    any monotone transform, and measurably does not: piece counts unchanged on all
    five example cores and the line displaced by at most 0.7 of a pixel. What it
    does do is force every crossing onto a pixel midpoint, inflating
    `Boundary Length` by 5–6%.
  * `overwrite = FALSE` (the default) errors, naming every existing clash, rather
    than appending a new `Run` — a cell can carry only one compartment label, so
    there is nowhere for a second run to go. Keep two mifs to compare two
    settings.
  * The density images and point patterns are discarded once the boundary and
    per-cell labels are derived, to avoid inflating the mif. Only the boundary
    polyline survives, in `mif$derived$density_boundary` (a named list, one data
    frame per sample), with the settings used attached as `"call_info"`.
  * `subset_mif()` rebuilds the mif from scratch and drops `derived` entirely, so
    the boundary slot and `Boundary Length` do not survive it, even though the
    three spatial columns ride along with the row filter. Split *after*
    subsetting, not before.
* **`plot_tissue_split()`** recomputes the density difference on demand, at plot
  time, from the settings `split_tissue()` recorded, and draws the *stored*
  boundary polyline (never a recontoured one) over a raster of the difference and
  a scatter of the compartment label. Unlike `plot_immunoflo()`, it returns a
  **named list of `ggplot` objects**, not the `mif` — each plot's raster can carry
  as much data as the sample itself, and attaching several to `mif$derived` would
  multiply the mif's size for no benefit.

## Read this first: results that were wrong

Two functions were producing incorrect output. If you have used either, re-run it.

* **`marker_freq_diff()` p-values were all wrong.** The Fisher contingency table
  was built with the compartment *total* as its second row instead of the count of
  marker-*negative* cells, so the margin double-counted the positives. On the
  shipped `example_spatial[[1]]` with `CD3..Opal.570..Positive` it returned
  `7.949897e-07` where the correct answer is `6.001431e-07`. Every p-value the
  function has ever produced is affected.

  Separately, marker selection used a substring match, so a marker whose name is a
  prefix of another silently absorbed the other's counts and `fisher.test()` was
  handed a 3x2 table whose result was stored as that marker's 2x2 p-value. Both
  `CD3..CD8.` and `CD3..CD8..FOXP3.` are real columns of the shipped data, so this
  was reachable with the documented examples.

* **`plot_immunoflo()`'s `cell_type` argument never had any effect.** It mapped the
  *string* you passed rather than the column, so a single shape was drawn for every
  cell and the legend had one entry named after the column. Any plot made with
  `cell_type` set showed Tumor and Stroma as identical.

* **`NN_G(edge_correction = "km")` returned a malformed table.** Columns were
  reordered by position, which assumed `Gest()` returns three columns; for `"km"` it
  returns five, so the output had no sample-id column, no `Marker` column, a leaked
  `theohaz` column, and roughly four times too many rows.

## Other changes that alter numbers

* `edge_correction` is **no longer silently downgraded**. Previously, exceeding the
  `big` cell-count threshold replaced whatever correction you asked for with
  `"none"`, and the tiled code path hard-coded translation regardless of the
  argument. Results for large samples change, for the better.
* `edge_correction = "none"` now uses **spatstat's `whist` binning**. `Kest()` has
  two internal binning paths that disagree on tied distances, and the one it takes
  when `"none"` is the only correction requested is not consistent with its own
  translation output. Values from this package are now bit-identical to
  `Kest(correction = c("none", "translation"))`. Only visible when pair distances
  land exactly on a radius, which the half-integer cell centres typical of HALO and
  Vectra exports produce constantly.
* **Permutations are now reproducible.** `set.seed()` previously had no effect:
  24 nested `mclapply()` calls omitted `mc.cores`, so they forked with their own
  RNG streams. Two runs at the same seed could differ by thousands. Results are now
  reproducible *and* independent of `workers`, which also means a fixed seed gives
  different numbers than 1.4.0 did.
* `workers` now means what it says. Only the outermost loop honoured it before, so
  `workers = 1` was not serial and the real process count could be `workers x 2 x 2`
  — over CRAN's two-core limit.
* `subset_mif()`'s `% ` columns are now true **percentages**. They held proportions
  (`sum/nrow`) while `marker_freq_diff()` put percentages in its `%` columns, so the
  two disagreed on what `%` meant. Values change by 100x.
* `dixons_s()` and `marker_freq_diff()` gain a correct `Run` column; see Bug fixes.

## Breaking changes

* **`compute_metrics()` removed.** It was never exported — no `@export`, absent from
  `NAMESPACE`, zero call sites — but it did have a man page and a pkgdown entry.
  Use `ripleys_k()`, `bi_ripleys_k()`, `NN_G()` or `bi_NN_G()` directly.
* **`bi_ripleys_k_WSI()` now raises an error** pointing at `bi_ripleys_k()`, which
  handles whole-slide images natively. Kept exported so the name still resolves and
  the error can explain the migration. `big`/`nlarge` have no equivalent because
  they no longer affect the statistic.
* **Output columns are unified across all seven metrics**, which will break code
  that reads the old names:

  | Old | New |
  |---|---|
  | `Theoretical G`, `Theoretical g` | `Theoretical CSR` |
  | `Permuted G`, `Permuted g`, `Permuted Interaction` | `Permuted CSR` |
  | `From`, `To` | `Anchor`, `Counted` |
  | `Degree of Correlation *`, `Degree of Interaction Permuted` | `Degree of Clustering *` |
  | `Permuted_larger_than_Observed` | `Permutations Larger than Observed` |

  Only the observed column keeps a statistic-specific name: `Observed K`,
  `Observed G`, `Observed g`, `Observed Interaction`. Columns a metric cannot fill
  are present and `NA` rather than absent, so `dplyr::bind_rows()` across metrics
  now lines up. `Exact CSR` is `NA` everywhere except Ripley's K — see below.
* `dixons_s()` output gains `Anchor`/`Counted` columns and its `tablaC` column names
  are no longer padded with spaces. Its `@return` documentation described three
  columns that never existed.
* `create_mif()` no longer computes the unused `sample_string`, and its `@return`
  no longer documents it. A consequence is welcome: mismatched integer/character
  ids are now accepted, so the `as.character()` call that every example carried is
  unnecessary.
* Defaults unified: `workers = 1` (was `6` in `bi_ripleys_k()`) and
  `overwrite = FALSE` (was `TRUE` in `bi_ripleys_k()`). `num_permutations` defaults
  are deliberately unchanged, since altering them would silently change results.
* `Depends: R (>= 4.1)`, up from an implausible `2.10`.

## Deprecations

Old argument names still work for the 2.x line, with a warning naming the
replacement. They will be removed in 3.0.0.

| Old | New |
|---|---|
| `keep_perm_dis` | `keep_permutation_distribution` |
| `nlarge` | `big` |
| `method` (`ripleys_k()`) | *removed* — the body never read it |
| `force` (`bi_ripleys_k()`) | *removed* — there is no longer a size limit to force past |

Unrecognised arguments are now an **error** rather than being silently ignored, so
`workerss = 4` fails loudly.

## Ripley's K rewritten

`ripleys_k()` and `bi_ripleys_k()` share one exact, memory-bounded engine. The old
implementation avoided an n-by-n distance matrix by tiling it, and paid for that
with the accuracy problems listed above. The engine instead materialises only the
cell pairs closer than `max(r_range)` and computes one edge weight per pair, so:

* **200,000 cells with the translation correction now run in ~1 s using ~740 MB**,
  where the matrix the tiling existed to avoid would have been 320 GB.
* Agreement with `Kest()`/`Kcross()` is 0 in most configurations and at worst
  9.1e-13 absolute on K values around 11,000 — floating point only.
* Because translation and isotropic edge weights depend only on a pair's
  displacement and the window, the pair list is computed once per sample and reused
  for every marker and permutation: about 26x faster than recomputing.
* `big` now bounds peak memory only. `big = 100` and `big = 1e9` give identical
  numbers.

One documented exception: bivariate **isotropic** can differ from `Kcross()` by
~1e-6 relative at a single radius. `Kest()` uses `closepairs()` and `Kmulti()` uses
`crosspairs()`, which compute the same distance with up to 1 ulp of difference, and
the isotropic weight jumps discontinuously when the circle of radius `d` passes
exactly through a window vertex — and since the window is the convex hull of the
cells, its vertices *are* cells. Translation, the default, is exact.

`bi_NN_G()`'s hand-rolled `rs` and `han` estimators, built on a full
`as.matrix(dist(...))`, are replaced by `spatstat.explore::Gcross()`. This was
numerically free — the values are identical — and drops memory from O(n²) to O(n).

`Exact CSR` is present for every metric but populated only for Ripley's K. There is
no closed form for `G` under random labelling, because `G` depends on the intensity
of the point set and not only its geometry: on the example data `Gest()` over all
cells gives 0.52 at a radius where the mean permuted G is 0.13. See `?NN_G`.

## Bug fixes

* `dixons_s(overwrite = FALSE)` — the default — always errored, taking `max()` of a
  column that did not exist yet.
* `marker_freq_diff()`'s `overwrite` branches were inverted: `overwrite = FALSE`
  destroyed the previous run on a populated object, and on a fresh one shipped
  `Run = -Inf`.
* `subset_mif()` errored when the first sample had two or fewer cells at the
  requested level, and — worse — silently duplicated the *previous* sample's summary
  row when a later one did, leaving a row matching no retained spatial frame. Its
  summary table also came back entirely character.
* `plot_immunoflo()` called `dev.off()` twice, so it errored on exit when no other
  device was open and silently closed the caller's device when one was. Its `path`
  argument was documented but never used, and a classifier with more than two levels
  errored at draw time.
* `pair_correlation()`, `bi_pair_correlation()` and `interaction_variable()` wrote
  appended runs to misspelled slots (`univaraite_…`, `bivaraite_…`,
  `derived_intraction_variable`), so `overwrite = FALSE` silently discarded them.
* The marker-existence guard `FALSE %in% mnames %in% colnames(spat)` could never
  fire, because `%in%` is left-associative.
* `spatial_exp_to_mif()` called `colData()`, `spatialCoords()` and `metadata()`
  unqualified from three packages that were not declared anywhere, so it failed
  unless you had attached them yourself. They are now in `Suggests`, checked for up
  front, and namespace-qualified. Its examples no longer download
  **VectraPolarisData** during `R CMD check`.
* `print.mif()` could not print a partially built object — `NULL[["id"]]` is an
  error in R — and returned `NULL` instead of `invisible(x)`.
* `merge_mifs()` accepted an empty list and then failed obscurely inside `seq(0)`.
  Its "No variables have been derived yet" message was unreachable.
* `merge_mifs()` called `dplyr::bind_rows()` on every `derived` slot regardless of
  type, so a **list**-valued slot — `spatial_plots`, and now `density_boundary` —
  was silently collapsed into a nameless data frame instead of being merged or
  erroring. List-valued slots are now concatenated, with a warning when the
  inputs' recorded settings disagree. A second bug in the same path prefixed every
  sample name with its source mif's index (`"1.S1"` rather than `"S1"`), because
  `do.call(c, parts)` read the list's own names as argument names.
* `create_mif()`'s "each item must be named" check could never fail.
* `is()` was used in three files without **methods** being declared; replaced with
  `inherits()`.
* `inst/CITATION` had `pages = 4584-4586` as unquoted arithmetic evaluating to `-2`,
  an empty `textVersion`, and a typo in the title.
* `example_spatial`'s documented names matched nothing in the data.
* Removed console noise: per-sample and per-permutation `cat()` output, unsuppressed
  dplyr join messages, a duplicate ggplot y-scale that warned once per sample, and
  **dixon**'s permutation counter.

## Internal

* `R/utils-helpers.R` is 998 lines shorter. The first-generation permutation engine
  (`uni_Rip_K`, `bi_Rip_K`, `uni_NN_G`, `bi_NN_G_sample` and their workers) was
  reachable only from `compute_metrics()`; `get_exactK()` could never run because it
  referenced an undefined variable; `dix_s_z()`/`dix_s_c()` were superseded by
  `dixons_s()`. A duplicate `get_bi_rows()` definition is gone.
* Four dependencies dropped, having become unused: **purrr**, **furrr**, **future**,
  **tidyselect**. `parallel`, `stats` and `utils` are now correctly declared.
* New internals: `k_pairs()`/`k_from_pairs()` (the K engine), `g_univariate()`/
  `g_bivariate()`, `standard_metric_cols()`, `write_derived()`,
  `add_cell_centres()`, `apply_deprecated_args()`.
* The observation window is now the convex hull of *all* cells in a sample by
  construction rather than by convention: the pair list is built from the full
  pattern and marker selection is a mask over it, so a marker subset cannot reach
  the window. `tests/testthat/test-window-invariant.R` enforces this with a fixture
  whose markers sit in one corner, making the otherwise-silent failure an ~8x error.
* `Roxygen: list(markdown = TRUE)` added; 11 man pages were shipping literal
  backticks and dead cross-references.
* `vignettes/spatialexperiment.Rmd` is now excluded from the CRAN tarball while
  remaining published on the pkgdown site.

## Note for maintainers

* The published site under `docs/` still lists `compute_metrics()` and does not link
  the SpatialExperiment article. It needs a `pkgdown::build_site()` refresh, which
  this release deliberately does not perform.
* `interaction_variable()`'s percentage divides an anchor-cell count by *anchor plus
  counted*, so `Observed Interaction` tops out near
  `100 * n_anchor / (n_anchor + n_counted)` rather than 100. Preserved exactly as it
  was and now documented in `?interaction_variable`, including how to rescale, since
  redefining a published metric is not a refactoring decision. Worth confirming
  against Steinhart et al.
