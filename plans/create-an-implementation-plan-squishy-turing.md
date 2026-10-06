# 2.1.1 — density-field options for `split_tissue()`

## Context

2.1.0 shipped `split_tissue()`/`plot_tissue_split()` computing the boundary as the exact zero level
set of the **raw** intensity difference `λ_class1 − λ_class2`. Running it against five ovarian WSIs
(294k–1M cells, frames 16k–40k × 21k–59k units) alongside the author's own prototype surfaced three
things:

1. **On real slides the raw difference produces contour fragments in tissue holes.** On
   `110361 A7` the zero set has 227 pieces, and **74% of the reported `Boundary Length` comes from
   pieces sitting in near-empty space** (median total density 1.46e-12 there vs 1.24e-3 elsewhere —
   nine orders of magnitude). `Boundary Length` is the quantity that goes into a survival model, and
   this contamination scales with how much empty space falls inside each slide's convex hull, i.e.
   with tissue fragmentation and scan area. That is a prep artifact acting as a confounder.
2. **The author's independent rewrite (`SplittingTissue/splitting_2.R`) converged on this package's
   design** — full-sample convex hull, sign-of-difference compartments, exact zero contour, psp +
   `nncross` interface — but rescales each class density to [0,1] before differencing and applies a
   ±1 hard threshold before contouring. They want the rescale available and defaulted on.
3. **The halfplane test fixture asserts exact contour geometry on a field that is zero to machine
   epsilon along the contour**, and it broke on a toolchain change with no code change (see Step 0).

Outcome: three new options on `split_tissue()`, a per-cell score column, and a portable fixture.

## Decisions (settled with the author)

| Decision | Value | Rationale |
|---|---|---|
| `rescale` | **`TRUE`** (new default) | Author's call. Changes `Boundary Length` for everyone; documented. |
| `min_density` | **`NULL`** (off) | Author's call — opt-in, so no current number moves. |
| `hard_threshold` | **`FALSE`** (option only) | Must never default on: staircase quantisation inflates length ~6%. |
| Per-cell score | **add `density_score`** | Already computed internally and discarded; a continuous score beats a 3-level factor for modelling. |
| Boundary storage | **unchanged** (data frame + `piece`) | Retains piece identity, plots with `geom_path`, survives `merge_mifs()`, no redundant window. `boundary_psp()` rebuilds the psp. |
| `eps` | **unchanged**, internal `sigma/8` | Author endorsed it; mask is grid-independent (20.9/20.9/21.0% masked across a 16× pixel range). No `resolution` argument. |
| Per-cell lookup | **unchanged** `im_value_at()` | `lookup.im` + `nearest.valid.pixel`; leaves 0 cells NA where `interp.im` leaves 291. |

## Step 0 — make the suite green first (separate commit, ships before anything else)

`tests/testthat/test-split-tissue.R` currently has **4 failures / 109 passes** with a clean working
tree. All four are in `"halfplane_mif recovers the exact half-plane boundary"`. Cause: the fixture's
two halves are exact mirror images, so at `x = 502.5` both densities are `1.2500000e-03` — equal to
16 significant figures, difference `6.5e-19`, ratio `2.6e-16`. That column of machine-epsilon zero
runs the full image height, so `contourLines()` decides "does the contour cross here?" 196 times on
the sign of noise and the line shatters into 47 fragments. The *geometry* is stable (x spans
[502.500, 502.642] — 3% of one 5-unit pixel; length 964.78 vs the ideal 975, −1.05%); only the piece
count and exact equalities are not.

Moving the class split does **not** fix this — verified: a 26-vs-24 split still leaves an
exactly-equal column (x = 522.5, difference −2.2e-19) and still fragments, 36 pieces. Any half-plane
split of a regular lattice preserves the symmetry. So assert what is robust:

- `:79` `pieces == 1` → drop
- `:80` `length(unique(bd$x)) == 1` → drop
- `:81` `unique(bd$x) == 500 + eps/2` → `expect_true(all(abs(bd$x - (500 + eps/2)) < eps))`
- `:83` length `tolerance = 0.01` → `tolerance = 0.02` (relative; −1.05% passes, and so does the
  exact 975 the old toolchain produced)
- `:86-88` (0 NAs, ground truth for all 2500 cells) — unchanged, they pass

Also reword `tests/testthat/helper-mif.R:124-132`, which documents "**1** contour piece … length
exactly `975`" as a guarantee. State instead that the field is antisymmetric about the split, so the
contour's *representation* is platform-dependent while its geometry is stable.

Kept separate because Step 5 re-measures the five real-data targets; a genuine regression must not be
able to hide inside an already-red suite.

## Step 1 — `compartment_diff()` gains the three options

`R/utils-density-boundary.R:90`. New signature:

```r
compartment_diff(pp, keep1, keep2, sigma, eps, filter_density,
                 rescale = TRUE, min_density = NULL, hard_threshold = FALSE)
```

Order of operations — each position is load-bearing:

1. `d1`, `d2` ← `density.ppp(pp[keep1/keep2], sigma, eps)`.
2. **`min_density` mask, on raw intensities before any rescale.** `tot <- d1$v + d2$v`;
   `thr <- min_density * (sum(keep1) + sum(keep2)) / area(Window(pp))`;
   `drop <- is.na(tot) | tot < thr`; set both images `NA` there. Everything needed is already in
   scope (`pp` carries the window), so no extra parameters. Expressed as a **fraction of the
   sample's own mean intensity** so it is scale-free across samples and coordinate units.
3. **`rescale`** each image to [0,1], with the zero-range guard below.
4. `unfiltered <- d1 - d2` — drives the per-cell sign **and** the new score column.
5. **`filter_density`** per image, *after* rescale, so a filter written against [0,1] (the author's
   `m[m < 0.025] <- NA`) works as written. `validate_filtered_im()` unchanged.
6. **`hard_threshold`** on the *filtered* image only: `v[v > 0] <- 1; v[v < 0] <- -1`. Never on the
   unfiltered field — that would destroy the score's magnitude and flatten the plot raster.
7. `list(raw = unfiltered, filtered = filtered)`.

**Zero-range guard is required, not defensive.** `scales::rescale()` returns `0.5` for zero-range
input (verified). If `class2` is absent its density is all-zero, so rescaling sends it to 0.5 and
`r1 - r2` crosses zero wherever `r1 = 0.5`, manufacturing a contour where `one_class_mif()` asserts
`nrow(boundary) == 0`:

```r
rescale01 <- function(v) {
  r <- range(v, na.rm = TRUE)
  if (!all(is.finite(r)) || diff(r) == 0) return(v)   # absent class stays 0, contour stays empty
  (v - r[1]) / diff(r)
}
```

## Step 2 — `split_tissue()` signature, validation, provenance, score column

`R/split_tissue.R:99`. Add `rescale = TRUE, min_density = NULL, hard_threshold = FALSE` as **real
formals, not via `...`** — the dots whitelist is deliberately only `filter_density` (`:109`) and
`test-split-tissue.R:436` pins that any other name errors.

Validation beside the existing block (`:136-159`): `rescale`/`hard_threshold` single non-`NA`
logical; `min_density` `NULL` or single finite `>= 0`.

**Score column is free.** `v` at `:210` is already `im_value_at(d$raw, xloc, yloc)`. Assign it in the
assembly loop (`:229-237`) as `density_score`. Bit-identical across `workers` because there is no RNG
— required by the whole-list `expect_identical(one$spatial, two$spatial)` at `:473`.

The clash check (`:162-166`) and the `overwrite = FALSE` error must name `density_score`, or a second
run would silently pass on a partially populated mif.

`call_info` (`:250-255`) gains `rescale`, `min_density`, `hard_threshold`. **Note:** assigning `NULL`
to a list element deletes it, so store "off" as `NA_real_` for `min_density` and translate at both
ends rather than round-tripping a `NULL`.

## Step 3 — `split_tissue_settings()` backward compatibility

`R/utils-density-boundary.R:231-264`. The `required` vector (`:232-234`) is the gate, and both `:242`
and `:263` do `settings[required]`, which **drops any field not listed** — so the three new fields
must be registered there or `plot_tissue_split()` will never see them.

But a mif saved by 2.1.0 has none of them, and adding them to `required` unconditionally makes
`split_tissue_settings()` error on every such mif. Backfill before the missing-field check:

```r
if (is.null(call_info$rescale))        call_info$rescale        <- FALSE
if (is.null(call_info$hard_threshold)) call_info$hard_threshold <- FALSE
if (is.null(call_info$min_density))    call_info$min_density    <- NA_real_
```

`FALSE`/`NA` are not arbitrary — they are what 2.1.0 actually computed, so an old mif replots
correctly rather than merely not erroring.

## Step 4 — `plot_tissue_split()` reproduces the field

`R/plot_tissue_split.R:135` — pass `cfg$rescale`, `cfg$min_density`, `cfg$hard_threshold` into
`compartment_diff()`. The `@param settings` list at `:25-31` must gain the three fields (it already
omits `spatialTIME_version`; fix while there).

Document one cosmetic caveat: when `raster_max_pixels` caps the raster, `eps_plot != cfg$eps`
(`:134`), so with `rescale = TRUE` the [0,1] map is computed over a coarser grid and the raster's
colour mapping shifts slightly. The drawn boundary is the stored polyline and is unaffected.

## Step 5 — tests

**Re-measure, do not hand-edit:** the five targets at `test-split-tissue.R:20-26` under the new
default. Run and paste.

**Confirm these survive `rescale = TRUE`** (all are properties rather than values, and all are
expected to hold — verify, don't assume):

- swap invariance (`:182-199`) — each image is rescaled by its own range independently of argument
  order, so `r1-r2 → -(r1-r2)` under swap; zero set identical, level order flips with it
- `filter_density = identity` bit-identical (`:390-400`)
- unfiltered drives the sign, so `mask_half` leaves `density_compartment` unchanged (`:402-418`)
- `one_class_mif()` empty boundary (`:103-113`) — this is the zero-range guard's regression test
- `test-window-invariant.R:171,178` call `compartment_diff()` **positionally** with 6 args; pass the
  new ones explicitly so the hand-rolled reference cannot drift from `split_tissue()`'s call

**New tests:**

- `density_score` on every frame, no NAs, sign agrees with `density_compartment`, bit-identical
  `workers = 1` vs `2`
- `rescale = FALSE` reproduces the 2.1.0 numbers exactly — the old targets become its regression
- `min_density`: `NULL` identical to `0`; a small value cuts pieces and length while changing few
  labels; a large value changes many. Pins both regimes, so the safe plateau is a test, not a note
- `hard_threshold = TRUE` leaves the piece count unchanged, increases length, leaves `density_score`
  untouched
- `call_info` round-trip for the three new fields (extend `:262-281`)
- a hand-built 2.1.0-style `call_info` lacking all three still plots (Step 3 regression)

## Step 6 — docs

- `R/utils-density-boundary.R:8-11` calls per-image rescaling "prototype bug #2 … moves and destroys
  contour components". Rewrite rather than delete: it is now an option, default on, and the finding
  still stands as the *cost* — rescaling makes the zero set depend on each sample's own density
  range, so `Boundary Length` stops being an absolute cross-sample criterion.
- `R/split_tissue.R:63-72` (`@section Boundary length`) — the convergence figures (6016.8 / 6756.5 /
  6799.1, −11.5% / −0.6%, "9 pieces at every resolution") were measured on the unrescaled field.
  **Annotate as measured with `rescale = FALSE`** rather than re-deriving: the argument they support
  is about `eps/sigma` ratios, which rescaling does not change. Mirrored in
  `man/split_tissue.Rd:107-118` (regenerated).
- `NEWS.md` — new 2.1.1 section, leading with the fact that `Boundary Length` changes for everyone
  under the new default: the three options, the score column, why `rescale` defaults on and what it
  costs, the `min_density` regimes with the A7 numbers, and the Step 0 fixture fix.
- Regenerate `man/` and pkgdown (`docs/reference/split_tissue.*`, `plot_tissue_split.*`, `news/`).
- `R/global.R` needs **no** entry — `split_tissue.R` uses only `[[ ]]`/`$` and
  `plot_tissue_split.R` uses `.data$` (already listed). Add only if `R CMD check` notes it.

## Accepted side effects

- `merge_mifs()` compares whole `call_info` with `identical()` (`R/merge_mifs.R:134-135`), so merging
  a 2.1.0 mif with a 2.1.1 one now warns "different settings". Correct — they *were* computed
  differently — but state it in NEWS.
- `Boundary Length` changes for every user under `rescale = TRUE`.

## Verification

1. **The overlay figure first, before any of this ships.** On `110361 A7` and `120189 C1` at
   sigma 500, draw the smooth and `hard_threshold` contours over the same density raster, zoomed on a
   region where both classes are near zero. This is the artifact the author asked to see, and it
   decides whether `hard_threshold` is worth keeping at all.
2. Re-measure and paste the five targets; `testthat::test_local()` green.
3. `devtools::check(document = FALSE)` — 0 errors / 0 warnings; only the pre-existing
   `_pkgdown.yml` top-level NOTE.
4. Cross-check on all five WSIs: `min_density = NULL` vs `0.05`, reporting pieces, `Boundary Length`,
   Interface count and cells-changed. Confirm the 1–10% plateau reproduces beyond A7.
5. `density_score`: sign agrees with `density_compartment` everywhere; `plot_tissue_split()`'s
   neutral fill band still sits on the drawn polyline.
6. Run the suite under **both** R installs (`/usr/local/bin/Rscript` 4.5.0 and the conda 4.6.1 once
   its deps are installed). This is the actual regression test for Step 0's portability claim.

## Effort

Step 0 ~30 min (isolated, ships first) · Steps 1–4 ~2–3 h (additive; `min_density` and
`hard_threshold` default off, so `rescale` is the only behaviour change) · Step 5 ~2–3 h (mostly
waiting on 5-WSI runs) · Step 6 ~1 h. **About a day.** The risk is not the code — it is that the
re-measured numbers must be generated rather than guessed, and that Step 0 lands first.

---

# Appendix — 2.1.0 plan (completed, shipped in 7803bef)

Retained because `tests/testthat/test-split-tissue.R:1-26` cites the measured numbers below.

## Context

`SplittingTissue/functions.R` holds an untracked prototype (`get_interface()`, `boundary_length()`,
`plot_get_interface()`) that splits a sample into tumor / stroma / interface compartments from the
difference of two kernel density estimates. It works, but it is not packageable as-is:

- **It extracts a band, not a boundary.** `get_interface()` keeps grid cells with
  `-boundary_threshold < z < boundary_threshold`, so the "boundary" is a point pattern of pixels whose
  width is a user-tuned parameter. The zero level set should be extracted exactly instead.
- **The one place that does use `contourLines()` is silently wrong.**
  `contourLines(im$xcol, im$yrow, im$v)` — an `im` stores `v` with `dim = c(length(yrow), length(xcol))`,
  and R's `contourLines()` has **no dimension check**, so it reinterprets the matrix transposed and
  returns a plausible-looking garbage contour. Verified: it returns 10 pieces where the correct call
  returns 5. This is why `boundary_length()` needed that "rotate back to match" `scales::rescale()`
  swap. Fix: `contourLines(x = im$xcol, y = im$yrow, z = t(im$v), levels = 0)`.
- **The per-image `scales::rescale(v, c(0,1))` before differencing moves the boundary.** Verified on
  `example_spatial[["TMA3_[9,K].tif"]]` at sigma 40: the raw difference gives 9 contour pieces /
  length 6756.5; the rescaled difference gives **7 pieces / 6723.8** — two whole components destroyed.
  (The zero set *is* invariant to any monotone transform of the ratio: `edge=FALSE` → 6756.6,
  `log(l1)-log(l2)` → 6756.4, relative risk `l1/(l1+l2)-0.5` → 6756.4. The per-image rescale is not
  one of those, because each image gets a different affine map.)
- **Too many knobs, and the memory-critical one is `dimyx`.** Also `dimyx` gives *non-square* pixels on
  a non-square window (verified 4× anisotropy on `owin(c(0,400), c(0,100))` at `dimyx = 64`), which
  distorts contour length.
- **It retains the density images and point patterns**, which would balloon the mif.

Outcome: two exported functions that annotate each cell with its density compartment, store only the
boundary polyline, record boundary length per sample, and recompute the density at plot time.

Answers already given: name is `split_tissue()`; `sigma` is required (no default); **all** cells get the
new columns, not just `class1`/`class2` cells; the plotting function ships in the same pass.

---

## Design decisions (settled, with evidence)

### Pixel size is derived, not exposed: `eps = sigma / 8`

`dimyx` and `boundary_threshold` disappear from the API. Resolution is set by `eps = sigma/8`, which
also guarantees square pixels. Justification, measured on `TMA3_[9,K].tif`, Tumor vs Stroma, sigma 40,
convex-hull window (1803 cells, area 1.37e6):

| eps | pixels | contour pieces | boundary length | bias |
|---|---|---|---|---|
| sigma/1 | 1,156 | 9 | 6016.8 | −11.5% |
| sigma/2 | 4,556 | 9 | 6516.5 | −4.2% |
| sigma/4 | 18,224 | 9 | 6692.1 | −1.6% |
| **sigma/8** | **72,628** | **9** | **6756.5** | **−0.6%** |
| sigma/16 | 290,512 | 9 | 6771.9 | −0.4% |
| sigma/32 | 1,162,048 | 9 | 6799.1 | ref |

**Topology is set by sigma, not resolution** — 9 pieces at every resolution. Resolution only adds a
sub-1% length bias past sigma/8, and it is a *consistent* bias at fixed sigma, so cross-sample
comparison is unaffected. There is nothing for the user to tune, which is the argument for hiding it.

`eps` is the only lever on memory and the relationship is quadratic and invisible
(`npixels = area(Frame(win)) * 64 / sigma^2`): a 40000×30000 WSI at sigma 25 needs ~1.2e9 pixels ≈
10 GB *per image*, ×3 images, ×`workers` forks. So a **hard pixel budget** is required — see
`check_density_budget()` below. Without it, a small sigma is an OOM kill rather than an error.

### The field is the raw intensity difference `lambda_class1 − lambda_class2`

No rescaling. Positive means class1. `density.ppp(edge = TRUE)` stays (the default); the zero set does
not care (6756.5 vs 6756.6) but the raster the plot draws does, so both functions must use the same call.

### Cells are labelled by the *sign of the field*, not by point-in-polygon

The contour **is** the zero level set, so `sign(diff)` at a cell's location is exactly "which side of the
boundary". This avoids all point-in-polygon logic on contours that are a mix of closed loops and open
pieces terminating at the mask edge (verified: real data returns both).

`spatstat.geom::lookup.im(diff, x, y, naok = TRUE)` returns `NA` for 9–16 cells per example core — cells
inside the hull whose *pixel centre* falls outside the mask. `spatstat.geom::nearest.valid.pixel()`
resolves all of them.

### `filter_density` affects the boundary, never the sign

Passed through `...`, a `function(im) im` applied to **each class image before differencing**. The
**filtered** difference drives the contour / length / interface distances; the **unfiltered** difference
drives the per-cell sign. So masking a tissue hole stops it inflating the boundary without ever
orphaning a cell to `NA`. Three caveats to document: the columns become mutually inconsistent inside
filtered regions (a cell can be `class1` while the polyline that separated it was filtered away); a
data-dependent filter forfeits the cross-sample comparability that `sigma`-is-required exists to
protect; and "per class image, before differencing" is not equivalent to filtering the difference for
any non-linear filter. Guard that the filter does not change the pixel grid, or `d1 - d2` misaligns.

### No `Run` column; `overwrite` is replace-or-refuse

`write_derived()` (`R/ripleys_k.R:252`) is structurally unusable here — it calls `nrow()`/`$Run`/
`bind_rows` on its argument, and our slot is a list keyed by sample. More fundamentally a cell can carry
only one compartment label and `mif$sample` one `Boundary Length` column, so appending has nowhere to
go. `overwrite = FALSE` errors listing every clashing location and says to keep two mifs to compare two
sigmas.

---

## Files

**Create:** `R/utils-density-boundary.R`, `R/split_tissue.R`, `R/plot_tissue_split.R`,
`tests/testthat/test-split-tissue.R`, `tests/testthat/test-plot-tissue-split.R`.

**Modify:** `R/utils-deprecate.R`, `R/merge_mifs.R`, `R/global.R` (only if check flags),
`tests/testthat/helper-mif.R`, `test-window-invariant.R`, `test-output-schema.R`,
`test-run-semantics.R`, `test-merge-mifs.R`, `test-subset-mif.R`, `test-deprecated-args.R`,
`NAMESPACE` + `man/` (regenerate), `DESCRIPTION` (prose only), `_pkgdown.yml`, `NEWS.md`,
`.Rbuildignore`, `.gitignore`.

**No new dependency.** `grDevices`, `spatstat.geom`, `spatstat.explore`, `ggplot2`, `scales`,
`RColorBrewer`, `dplyr`, `parallel`, `stats`, `utils` are already in Imports. The prototype's `ggpubr`
and `arrow` must not be added — the figure is buildable in base ggplot2 (see below).

---

## Step 1 — `R/utils-density-boundary.R` (internals, `@keywords internal` + `@noRd`)

```r
DENSITY_EPS_DIVISOR <- 8L    # see the convergence table in ?split_tissue
DENSITY_MAX_PIXELS  <- 4e7   # ~320 MB per im; three ims live at once, per worker

density_pixel_size(sigma)                        # sigma / DENSITY_EPS_DIVISOR
check_density_budget(win, eps, sigma, label)      # errors naming the minimum viable sigma
class_mask(values, level)                         # !is.na(values) & values == level
compartment_diff(pp, keep1, keep2, sigma, eps, filter_density)   # list(raw = im, filtered = im)
validate_filtered_im(before, after)               # must be an im on an identical grid
zero_contour(im, sample_id, label)                # data.frame(<sample_id>, piece, x, y)
empty_boundary_df(sample_id)                      # 0 rows, correct column types
boundary_psp(boundary_df, win)                    # psp over Frame(win); 0-segment safe
boundary_length(S)                                # sum(spatstat.geom::lengths_psp(S))
im_value_at(im, x, y)                             # lookup.im + nearest.valid.pixel fallback
split_tissue_settings(mif, settings = NULL)       # call_info accessor with actionable error
```

Each of `zero_contour`, `boundary_psp`, `boundary_length`, `im_value_at` encodes one prototype bug, so
each earns a named function with a comment block.

```r
zero_contour <- function(im, sample_id, label) {
  # t(v) is load-bearing. An `im` stores v as dim = c(length(yrow), length(xcol)) --
  # rows are y -- and grDevices::contourLines() has NO dimension check, so passing
  # im$v silently transposes the field and returns a plausible wrong answer.
  cs <- withCallingHandlers(
    grDevices::contourLines(x = im$xcol, y = im$yrow, z = t(im$v), levels = 0),
    warning = function(w) if (grepl("all z values are NA", conditionMessage(w)))
      invokeRestart("muffleWarning"))
  cs <- cs[vapply(cs, function(c) length(c$x) >= 2L, logical(1))]   # 1-vertex -> NA coords
  if (!length(cs)) return(empty_boundary_df(sample_id))
  out <- do.call(rbind, lapply(seq_along(cs), function(i)
    data.frame(label, piece = i, x = cs[[i]]$x, y = cs[[i]]$y, stringsAsFactors = FALSE)))
  names(out)[1] <- sample_id
  out$piece <- as.integer(out$piece)
  out
}

im_value_at <- function(im, x, y) {
  v  <- spatstat.geom::lookup.im(im, x, y, naok = TRUE)
  na <- is.na(v)
  if (any(na)) {
    # 9-16 cells per example core land in a pixel whose CENTRE is outside the mask,
    # even though the cell is inside the hull (the hull is built FROM the cells).
    nv <- spatstat.geom::nearest.valid.pixel(x[na], y[na], im)
    v[na] <- im$v[cbind(nv$row, nv$col)]
  }
  v
}

check_density_budget <- function(win, eps, sigma, label) {
  a <- spatstat.geom::area(spatstat.geom::Frame(win))
  n <- a / eps^2
  if (n > DENSITY_MAX_PIXELS) {
    sigma_min <- DENSITY_EPS_DIVISOR * sqrt(a / DENSITY_MAX_PIXELS)
    stop(sprintf(
      "Sample \"%s\" needs a %s-pixel density grid at sigma = %g.\n  The internal pixel size is sigma/%d, so halving sigma quadruples memory.\n  Use sigma >= %.0f for this sample, or split the image first.",
      label, format(round(n), big.mark = ","), sigma, DENSITY_EPS_DIVISOR,
      ceiling(sigma_min)), call. = FALSE)
  }
  invisible(n)
}
```

`contourLines()` returns vertices already ordered along each polyline, so `geom_path(aes(group = piece))`
needs no re-sorting. `sum(lengths_psp(S))` equals the naive `sum(sqrt(diff(x)^2 + diff(y)^2))` to the
last digit (verified), so either is fine; use `lengths_psp`.

## Step 2 — `R/split_tissue.R`

```r
split_tissue <- function(mif, classifier, class1, class2, sigma, interface_width,
                         workers = 1, overwrite = FALSE, xloc = NULL, yloc = NULL, ...)
```

**Dots.** `deprecated_arg_map()`'s `switch` falls through to `common`, so add explicit
`split_tissue = character(0)` and `plot_tissue_split = character(0)` branches in
`R/utils-deprecate.R:23` — otherwise `keep_perm_dis` is silently accepted. Then follow the
`R/pair_correlation.R:53-56` split-the-dots pattern with a whitelist of exactly `"filter_density"`;
anything else is an error naming the offender, plus a check for unnamed dots.

**Validation.** `inherits(mif, "mif")` (house-style message); reject `list(NA)` spatial from
`create_mif(spatial_list = NULL)` (`R/create_mif.R:80-82`); `missing(sigma)` / `missing(interface_width)`
error naming the argument; `sigma > 0`, `interface_width >= 0`, both length-1 finite; `classifier` /
`class1` / `class2` length-1 character; **`class1 == class2` errors** (identically-zero field); **either
class named `"Interface"` errors** (reserved third level, would silently merge categories);
`filter_density` must be a function.

**Clash check** before doing any work — collect which of the three destinations already exist
(`density_compartment`/`refined_density_compartment` in any spatial frame, `derived$density_boundary`,
`Boundary Length` in `mif$sample`) and error if `!overwrite`. On `overwrite = TRUE`, **drop
`mif$sample[["Boundary Length"]]` before the join** or `left_join` yields `.x`/`.y`.

**Per-sample loop** — follows `R/ripleys_k.R:100-178` **minus the seeding**: there is no RNG anywhere in
this function, so omit `seeds`/`set.seed` with a one-line comment saying so (otherwise the omission
looks like the pre-2.0.0 reproducibility bug).

```r
res <- parallel::mclapply(seq_along(mif$spatial), function(sample_i) {
  spat  <- add_cell_centres(mif$spatial[[sample_i]], xloc, yloc)   # R/ripleys_k.R:217
  label <- as.character(spat[[mif$sample_id]][1])
  cls   <- spat[[classifier]]
  if (is.null(cls)) stop("Column \"", classifier, "\" not found in sample \"", label, "\".",
                         call. = FALSE)

  # Window from EVERY cell in the sample -- pinned by test-window-invariant.R. Class
  # selection is a logical MASK over the full ppp; pp[keep] preserves the window.
  win <- spatstat.geom::convexhull.xy(spat$xloc, spat$yloc)
  pp  <- spatstat.geom::ppp(spat$xloc, spat$yloc, window = win, check = FALSE)
  check_density_budget(win, eps, sigma, label)

  # class_mask() is !is.na(v) & v == level, NOT v == level: a logical subscript with
  # NA makes pp[keep] error ("Index out of bounds in [.ppp"), and unclassified cells
  # are realistic. The two masks need not partition the sample.
  keep1 <- class_mask(cls, class1); keep2 <- class_mask(cls, class2)
  if (!any(keep1) || !any(keep2)) warning(...)   # per-sample, names the counts

  d  <- compartment_diff(pp, keep1, keep2, sigma, eps, filter_density)
  bd <- zero_contour(d$filtered, mif$sample_id, label)
  S  <- boundary_psp(bd, win)
  v  <- im_value_at(d$raw, spat$xloc, spat$yloc)          # UNFILTERED drives the sign
  list(label = label,
       code      = ifelse(v > 0, 1L, ifelse(v < 0, 2L, NA_integer_)),
       interface = spatstat.geom::nncross(pp, S)$dist <= interface_width / 2,
       boundary  = bd,
       length    = boundary_length(S))
}, mc.cores = workers, mc.preschedule = FALSE)
```

Return integer codes plus a logical, **not** factors and **not** the mutated spatial frame: building the
factor in the parent guarantees identical levels across samples, and shipping a 3800×51 frame back
through the fork is pure waste.

**Add the `mclapply` error check the metric functions lack** — with `workers > 1` `mclapply` returns
per-element `try-error` objects instead of raising, which here would leave the mif with the new columns
on only *some* frames. Detect with `vapply(res, inherits, logical(1), "try-error")` and re-raise naming
the failing sample.

**Assemble in the parent:**

```r
lev <- c(class1, class2); levr <- c(class1, "Interface", class2)
for (i in seq_along(mif$spatial)) {
  lab <- lev[res[[i]]$code]
  mif$spatial[[i]][["density_compartment"]] <- factor(lab, levels = lev)
  lab[res[[i]]$interface] <- "Interface"
  mif$spatial[[i]][["refined_density_compartment"]] <- factor(lab, levels = levr)
}
```

Building `refined` by *overwriting* `lab` means the two columns cannot disagree by construction. Do not
write it as a nested `ifelse` — verified that `ifelse(dist <= w/2, "Interface", ifelse(v > 0, ...))`
leaves 8 cells `NA` on `TMA3_[9,K]`, because the outer `ifelse` cannot rescue an inner `NA` for cells
that are not Interface.

**Exactly zero → `NA`.** Verified 0 of 1803 cells and 0 of 72,628 pixels are exactly 0 on real data, so
any tie-break rule would be an untested branch that silently biases one compartment. Exact zero is
essentially only reachable when *both* classes are absent — i.e. the user's `class1`/`class2` don't match
their data — and an all-`NA` column says that loudly. A `diff == 0` cell also lies *on* the zero set, so
it becomes `"Interface"` for any `interface_width > 0`; `NA` survives into the refined column only when
`interface_width == 0` or no contour exists. Document that pairing. Every cell still gets both columns —
`NA` is a value, not an omission.

**`mif$derived$density_boundary`** — a named list, one element per spatial sample, named by
`names(mif$spatial)` (matching `plot_immunoflo`'s tested contract), each element a data frame
`<sample_id>, piece (int), x, y`, 0 rows with correct types when no boundary exists. Warn when
`names(mif$spatial)[i]` disagrees with `spat[[sample_id]][1]` — that mismatch is exactly what makes the
summary join drop rows.

**Provenance** — `attr(mif$derived$density_boundary, "call_info")`:

```r
list(classifier =, class1 =, class2 =, sigma =, eps =, interface_width =,
     xloc =, yloc =,               # REQUIRED: without these the plot cannot rebuild the ppp
     filter_density =,             # NULL when absent
     sample_id = mif$sample_id,
     spatialTIME_version = as.character(utils::packageVersion("spatialTIME")))
```

`xloc`/`yloc` matter: the `(XMax+XMin)/2` fallback applies only when *both* are `NULL`, so a user who
passed explicit columns would otherwise get a different point pattern at plot time. Attributes survive
`[[` and `saveRDS` but are **dropped by `[`** (verified) and by `bind_rows`, hence
`split_tissue_settings(mif, settings = NULL)` — returns `settings` if supplied, else the attribute, else
errors telling the user to re-run `split_tissue()` or pass `settings` explicitly, and validates that a
supplied partial `settings` has every required name.

Rejected alternatives: `list(boundary =, call_info =)` violates the "one element per sample" shape and
makes every access `$boundary[[i]]`; a separate settings slot can't be a data frame (`filter_density` is
a closure) so it hits the same `merge_mifs` exposure plus can diverge; per-element attributes are two
sources of truth for no gain.

**`Boundary Length` into `mif$sample`** — build with `check.names = FALSE`, rename column 1 to
`mif$sample_id`, `dplyr::left_join`. Four guards, all reachable:

- duplicate sample ids across spatial frames → error (a left join would fan out `mif$sample` rows);
- id class mismatch → error naming both classes rather than coercing. `create_mif()` was *deliberately*
  relaxed to accept mismatched id types (`R/create_mif.R:63-69`), so this will happen;
- spatial frames unmatched in `mif$sample` → warn naming them (a left join drops them silently);
- `0` vs `NA`: a sample *with* spatial data and no contour gets `0` (a real measurement); a sample in
  `mif$sample` with no spatial frame gets `NA`. On the shipped data that is 5 values and 224 `NA`s.

No piece count column (`length(unique(b$piece))` is one line, and topology is resolution-invariant). No
area-normalised column, but document in `@details` that raw length is not scale-free — compare
`Boundary Length / sqrt(area)` across cores of different size.

**Roxygen.** `@param class1,class2` as a pair stating the sign convention explicitly ("the difference is
`class1` minus `class2`, so positive means `class1`"); `@param xloc,yloc` as a pair per house style;
`@param sigma` explains why there is no default; `@param interface_width` states it is a *diameter* (band
is ±half) in the same coordinate units as `sigma`; `@param ...` documents `filter_density` briefly.
Sections: **What is and is not kept** (densities and point patterns are discarded — that is the point;
use `plot_tissue_split()`), **Choosing sigma** (units differ per platform; memory is quadratic in
1/sigma), **Boundary length** (paste the convergence table verbatim as the justification for the hidden
`eps`, and note topology is identical at every resolution so no knob is needed).

## Step 3 — `R/plot_tissue_split.R`

```r
plot_tissue_split <- function(mif, which = NULL,
                              compartment = c("refined_density_compartment",
                                              "density_compartment"),
                              panels = c("both", "density", "compartment"),
                              colors = NULL, point_size = 0.4,
                              raster_max_pixels = 5e5, workers = 1,
                              filename = NULL, path = NULL, settings = NULL, ...)
```

**One faceted ggplot per sample beats arranging two panels**, and it removes the need for `ggpubr` /
`patchwork` / `gridExtra` (none are dependencies). ggplot2 keeps `fill` and `colour` on separate scales,
so give each layer its own data with a `panel` column and facet on it. Verified end to end on
ggplot2 4.0.3: 3 layers, 54,752 / 1,803 / 3,468 rows, **no warnings**.

```r
ras <- as.data.frame(diff_im); ras$panel <- "Density difference"     # x, y, value; NA pixels dropped
pts <- data.frame(x =, y =, compartment =, panel = "Compartment")
bd2 <- rbind(transform(bd, panel = "Density difference"),
             transform(bd, panel = "Compartment"))                   # boundary on both panels

ggplot2::ggplot() +
  ggplot2::geom_raster(data = ras, ggplot2::aes(.data$x, .data$y, fill = .data$value)) +
  ggplot2::scale_fill_gradient2(name = sprintf("%s - %s\ndensity", class1, class2),
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-m, m), oob = scales::squish) +     # m = max(abs(value))
  ggplot2::geom_point(data = pts,
    ggplot2::aes(.data$x, .data$y, colour = .data$compartment), size = point_size) +
  ggplot2::scale_colour_manual(NULL, values = colors, drop = FALSE, na.value = "grey70") +
  ggplot2::geom_path(data = bd2,
    ggplot2::aes(.data$x, .data$y, group = .data$piece), colour = "black", linewidth = 0.3) +
  ggplot2::facet_wrap(~ panel) + ggplot2::coord_equal() +
  ggplot2::scale_y_reverse(breaks = scales::pretty_breaks(5)) +
  ggplot2::ggtitle(paste0("ID: ", label)) +
  ggplot2::theme_bw(base_size = 18) +
  ggplot2::theme(axis.title = ggplot2::element_blank())
```

Two things that are easy to get wrong:

- **The fill scale's zero must coincide with the drawn polyline** — `midpoint = 0` *and* symmetric
  `limits = c(-m, m)` with `oob = scales::squish`. Without symmetric limits the neutral colour lands
  somewhere other than 0 and the raster implies a boundary in a different place from the line, which is
  precisely what this figure exists to check. (The prototype's `gradient2(mid = "red", midpoint = 0.5)`
  over a per-image `rescale(0,1)` has no defined relationship to the zero set at all.)
- **Draw the STORED boundary, never a recomputed contour.** Only the *raster* is recomputed, so its
  resolution is purely cosmetic — hence `raster_max_pixels`, recomputing at
  `eps_plot <- max(eps, sqrt(area(Frame(win)) / raster_max_pixels))`. Above the cap the raster is a
  display approximation while the polyline stays exact. Without it a WSI is unrenderable.

`.data$` everywhere (as `R/plot_immunoflo.R:132` already does) so `R/global.R` need not grow. `colors`
default: `class1`/`class2` from `RColorBrewer::brewer.pal(3, "Set1")[1:2]`, `"Interface"` `"grey20"`,
`na.value = "grey70"` so tie/NA cells are visible rather than silently dropped. `drop = FALSE` keeps
panels comparable across samples. Keep the `filename`/`path` PDF block from
`R/plot_immunoflo.R:157-176` verbatim, including the single `on.exit(grDevices::dev.off(), add = TRUE)`
and the trailing-`.pdf` strip. Use `parallel::mclapply` with `workers` (`pbmcapply` appears only in
`plot_immunoflo`). Same dots-typo guard (whitelist empty → any `...` is an error).

**Returns a named list of ggplots, not the mif** — a deliberate deviation from `plot_immunoflo`, called
out in the roxygen and `NEWS.md`. Reasons: a ggplot captures its input data, so each plot carries a
~54k-row / 1.3 MB raster frame and putting five of those in the mif triples its size (untenable for a
WSI cohort); `merge_mifs()` is already broken by exactly this shape for `spatial_plots` and a second
list slot doubles the exposure; `create_mif`'s own `@return` says `derived` is "data derived using the
MIF object" and a plot is a view; and `plots[["TMA3_[9,K].tif"]]` composes better than
`mif$derived$tissue_split_plots[[4]]`.

## Step 4 — `merge_mifs()` list-slot guard

`R/merge_mifs.R:108-118` unconditionally `bind_rows()`es every derived slot. Verified this "succeeds"
on `density_boundary`, silently collapsing per-sample frames into one nameless 2-row table — the worst
outcome. It already mangles `spatial_plots`. Replace the body with a type check: `bind_rows` when every
part is a data frame, else `do.call(c, parts)` with a duplicate-name error under `check.names`, and
carry `call_info` forward from the first mif with a **warning when settings disagree** (merging cohorts
run at different sigma yields a `Boundary Length` column that is not one measurement, and nothing else
records that).

## Step 5 — Metadata

- **`R/utils-deprecate.R`** — `split_tissue = character(0)`, `plot_tissue_split = character(0)` branches
  *before* the fall-through default, with a comment explaining why they must not inherit `common`.
- **`R/global.R`** — aim for **zero additions** via `.data$` and base `[[<-`. If check still reports "no
  visible binding", add only the names it reports.
- **`_pkgdown.yml`** — new group after "Other Spatial Statistics":
  `- title: Tissue Architecture / desc: Density-based segmentation of a sample into tissue compartments
  / contents: - split_tissue`; append `plot_tissue_split` to the existing `Visualization` group.
- **`.Rbuildignore`** — add `^SplittingTissue$`. It is currently absent, so the prototype would ship in
  a CRAN tarball. **`.gitignore`** — add `SplittingTissue/` so `git add -A` can't commit it. Keep the
  directory locally as provenance until both functions land, then delete it (`1.0.splitting.R` needs
  `arrow` and absolute `data/parquet/` paths, so nobody else can run it anyway).
- **`DESCRIPTION`** — one clause in `Description:`. No Imports/Suggests change.
- **`NEWS.md`** — a new `# spatialTIME 2.1.0` section in the existing register (bolded claim, then the
  concrete numbers and the dataset it was reproduced on): the two functions; `sigma` required and why;
  `eps = sigma/8` with the convergence table; the two spatial columns and the reserved `"Interface"`
  level; the new `mif$sample` column; `overwrite` has no `Run` and why; `plot_tissue_split()` returns a
  list rather than the mif, unlike `plot_immunoflo()`, and why; the `merge_mifs()` list-slot fix (which
  also fixes `spatial_plots`, a pre-existing bug); the `subset_mif()` interaction.

## Step 6 — Tests

**Two new fixtures in `tests/testthat/helper-mif.R`:**

`halfplane_mif()` — deterministic lattice `seq(10, 990, by = 20)` squared (2500 cells),
`Classifier.Label = ifelse(x < 500, "Tumor", "Stroma")`. The difference field is exactly antisymmetric,
so at **sigma = 40 / eps = 5** I verified: **1 contour piece**, `range(x)` a single value `502.5`
(= 500 + eps/2, the half-pixel discretisation offset), **length exactly 975.00 = 980 − eps** where 980
is the hull height, 0 lookup-NAs, 0 exact zeros, and labels matching ground truth for **all 2500 cells**.
Assert `all(abs(bd$x - 500) <= eps)`, `length(unique(bd$x)) == 1`, and
`expect_equal(len, height - eps, tolerance = 0.01)` — derived, not fudge-factored. Use sigma 40, **not
60**: at sigma 60 numerical noise splits it into 2 pieces (970.86). Comment that in the fixture so
nobody "improves" the bandwidth.

`one_class_mif()` — all `Classifier.Label == "Tumor"`; exercises the no-boundary path.
(`nncross(pp, <0-segment psp>)` returns `Inf`, verified — so no special-casing is needed and nothing
becomes Interface.)

**`test-split-tissue.R`**

- *Structure*: returns a `mif`; both columns on every frame; `levels()` exactly `c(class1, class2)` and
  `c(class1, "Interface", class2)` **in that order**; `nrow` unchanged, no other column touched;
  `derived$density_boundary` is a **list, not a data frame**, named by `names(mif$spatial)`, elements
  with columns `c(sample_id, "piece", "x", "y")` and classes character/integer/numeric/numeric.
- *`halfplane_mif()`*: the assertions above, plus every Interface cell has
  `abs(x - 500) <= interface_width/2 + step` and every non-Interface cell `>= interface_width/2 - step`.
- *Real data*: **per-`piece` polyline length summed equals `Boundary Length` to 1e-9** — one assertion
  covering both the vertex ordering and the reported value, and exactly what the transpose bug broke;
  the interface reproduced independently by rebuilding the psp from the stored frame and re-running
  `nncross` (exact equality); `refined == density` wherever `refined != "Interface"`;
  `interface_width = 0` gives at most a handful of Interface cells and the count is **monotone
  non-decreasing** over `c(0, 25, 50, 100, 200)`; **swap invariance** — swapping `class1`/`class2` gives
  identical `Boundary Length` and geometry with labels exchanged (the prototype computed `str - tum`
  while its signature read `get_interface(tumor_ppp, stroma_ppp, ...)`, so its labels were inverted
  relative to argument order); `sum(is.na(density_compartment)) == 0` on all five cores.
- *`Boundary Length`*: exactly 5 non-`NA` of 229; a sample with data and no contour gets `0` not `NA`;
  two `overwrite = TRUE` runs create no `.x`/`.y`; id class mismatch errors naming both classes;
  duplicate ids error rather than fanning out `mif$sample`.
- *Provenance*: `call_info` round-trips every field including `xloc`/`yloc` and the closure.
- *`overwrite`*: second call with `FALSE` errors and the message names all three clashing locations;
  `TRUE` replaces with column count and `nrow` unchanged; clash detection fires on a hand-built partial
  state where only the spatial columns exist.
- *Errors and edges*: missing `sigma`; missing `interface_width`; `class1 == class2`; a class named
  `"Interface"`; `classifier` absent (error names the sample); **`NA` in the classifier column** (the
  `pp[keep]` "Index out of bounds in [.ppp" trap, verified); both classes absent → warning + length 0;
  `filter_density = function(im) im` **bit-identical** to omitting it; a filter NA-ing a half-plane →
  shorter boundary with `density_compartment` **unchanged** (proves unfiltered drives the sign); a filter
  returning a coarser `im` → grid error; returning a matrix → `im` error; `workerss = 4` → "Unknown
  argument"; `sigma = 1` → pixel-budget error naming a minimum sigma; explicit `xloc`/`yloc` equals the
  `XMin/XMax` fallback; `workers = 2` **exactly** equals `workers = 1` (no RNG);
  `create_mif(spatial_list = NULL)` mif → clear error not a subscript failure.

**`test-plot-tissue-split.R`** — named list of ggplots, one per sample, correctly named, **not** a mif,
with `mif$derived` unchanged; `which = 2` returns one plot; built layer row counts (raster > 0, points
== `nrow(spatial)`, path == `2 * nrow(boundary)` for `panels = "both"`); **the path is the stored
boundary** — shift `mif$derived$density_boundary[[1]]$x` by 50, rebuild, assert the drawn path moved;
fill `limits[1] == -limits[2]`; `panels = "compartment"` drops the raster and `"density"` the points;
`compartment = "density_compartment"` gives 2 colour levels vs 3; y reversed and exactly one y scale
(mirroring the `plot_immunoflo` assertions at `test-plot-immunoflo.R:109-110`); missing
`density_boundary` → error telling the user to run `split_tissue()`; stripped `call_info` → error
mentioning `settings`, then supplying `settings` works; partial `settings` → error naming the missing
field; `filename`/`path` writes a file, no doubled `.pdf`, `dev.list()` unchanged;
`expect_no_warning(ggplot_build(g))`; `raster_max_pixels = 100` → far fewer raster rows, **identical**
path rows.

**Cross-cutting edits**

- **`test-window-invariant.R`** — the important one. Add a fixture with a three-level classifier where
  `class1 ∪ class2` sits in one corner and a third level fills the rest, so
  `area(hull(all)) / area(hull(class1 ∪ class2)) > 4` (assert that guard first, as the file already
  does). Then assert the boundary and length match the full-window reference and **differ** from the
  `class1 ∪ class2`-window computation beyond tolerance. "Build a ppp from just the two classes" is the
  shorter, wrong refactor and its failure would be silent.
- **`test-output-schema.R`** — do **not** add to `metric_table()` (no `r`, no `iter`, no `Observed *`).
  Add a test asserting the slot is a list and not a data frame, so the exclusion reads as intentional.
- **`test-run-semantics.R`** — likewise not in `run_writers()`. Add a test that no `Run` appears
  anywhere and that `overwrite = FALSE` on a populated mif **errors** instead of appending.
- **`test-merge-mifs.R`** — two mifs each with `density_boundary`: merged slot is a list of 4 named
  elements with per-sample frames intact; duplicate sample under `check.names = TRUE` errors; differing
  `call_info` warns. Add a `spatial_plots` case, since the same guard fixes that pre-existing bug.
- **`test-subset-mif.R`** — pin current behaviour: after `split_tissue()`, `subset_mif()` **keeps** the
  two spatial columns (they ride along with the row filter) but **drops** the boundary slot and
  `Boundary Length`, because it calls `create_mif()` afresh (`R/subset_mif.R:102-104`). Don't change it;
  document "subset first, then split" in both man pages — a `density_compartment` column surviving into
  a mif whose boundary slot is gone is a genuinely confusing state.
- **`test-deprecated-args.R`** — `split_tissue(..., keep_perm_dis = TRUE)` and the same for
  `plot_tissue_split` must **error**, not warn. Guards the `switch` fall-through.

---

## Verification

```r
devtools::document(); devtools::load_all()
mif <- create_mif(example_clinical, example_summary, example_spatial,
                  "deidentified_id", "deidentified_sample")
out <- split_tissue(mif, "Classifier.Label", "Tumor", "Stroma",
                    sigma = 40, interface_width = 100, overwrite = TRUE)
```

Must reproduce exactly (measured twice, ~3 s total at `workers = 1`) — cells / pieces / boundary length:
`TMA1_[3,B]` 3803 / 10 / 3238.2 · `TMA2_[3,B]` 3008 / 4 / 5189.2 · `TMA3_[7,B]` 1850 / 6 / 2346.8 ·
`TMA3_[9,K]` 1803 / 9 / 6756.5 · `TMA3_[8,U]` 2318 / 8 / 4469.4.

```r
# 1. self-consistency: stored polyline length == the reported number
naive <- vapply(out$derived$density_boundary, function(b)
  if (!nrow(b)) 0 else sum(vapply(split(b, b$piece), function(p)
    sum(sqrt(diff(p$x)^2 + diff(p$y)^2)), numeric(1))), numeric(1))
stopifnot(all(abs(naive - stats::na.omit(out$sample$`Boundary Length`)) < 1e-9))

# 2. no orphaned cells (the 9/16/13/10/14 lookup NAs must all be resolved)
stopifnot(!any(vapply(out$spatial, function(s) anyNA(s$density_compartment), logical(1))))

# 3. agreement with Classifier.Label -- expect 83.7-96.3%, NOT 100%
vapply(out$spatial, function(s)
  mean(as.character(s$density_compartment) == s$Classifier.Label), numeric(1))

# 4. sign convention: swapping classes must not move the boundary
sw <- split_tissue(mif, "Classifier.Label", "Stroma", "Tumor", sigma = 40,
                   interface_width = 100, overwrite = TRUE)
stopifnot(isTRUE(all.equal(sw$sample$`Boundary Length`, out$sample$`Boundary Length`)))

# 5. join shape; 6. determinism under forking (no RNG, so exact)
stopifnot(nrow(out$sample) == 229, sum(!is.na(out$sample$`Boundary Length`)) == 5)
stopifnot(identical(out$spatial, split_tissue(mif, "Classifier.Label", "Tumor", "Stroma",
          sigma = 40, interface_width = 100, workers = 2, overwrite = TRUE)$spatial))

# 7. plots
p <- plot_tissue_split(out)
stopifnot(identical(names(p), names(out$spatial)))
invisible(lapply(p, ggplot2::ggplot_build))          # must be warning-free
ggplot2::ggsave(file.path(tempdir(), "check.png"), p[[4]], width = 14, height = 7, dpi = 150)

# 8. the memory guard errors instead of hanging
try(split_tissue(mif, "Classifier.Label", "Tumor", "Stroma", sigma = 1,
                 interface_width = 10, overwrite = TRUE))

# 9. suite + check
testthat::test_local(); devtools::check(document = FALSE)
```

```bash
# 10. the prototype must not ship
R CMD build . && tar -tzf spatialTIME_*.tar.gz | grep -i splittingtissue    # expect nothing
```

On the saved PNG, check three things against the target figure: the neutral band of the diverging fill
sits **on** the drawn polyline; Interface cells form a band of the requested width; the two panels are
aligned and equal-aspect. Watch `R CMD check` for "no visible binding" notes (the `global.R` decision
returning) and a `\usage` line-width note on the long `plot_tissue_split` signature.
