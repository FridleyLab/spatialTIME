# spatialTIME 2.0.0

A cleanup and correctness release. Roughly 1,100 lines of unreachable code are
gone, the count-based measures now agree with **spatstat** to floating-point
precision on samples of any size, every metric returns the same columns, and test
coverage went from 4 assertions to 465 (91%).

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
