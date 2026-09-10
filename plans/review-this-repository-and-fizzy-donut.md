# spatialTIME 2.0.0 — Cleanup & Refactor

## Context

`spatialTIME` has accumulated three generations of spatial-statistic implementations. Each
generation was written to be faster or more memory-frugal than the last, but none of the
earlier ones were removed. The result:

- **~82% of `R/utils-helpers.R` (815 of 998 lines) is unreachable from the exported API.** The
  entire gen-1 permutation engine hangs off `compute_metrics()`, which has a man page and a
  pkgdown entry but **is not in `NAMESPACE`** — it is not exported and never called.
- **Ripley's K contains ~250 lines of hand-rolled distance/edge-correction/tiling code** whose
  only purpose is to avoid materialising an `n × n` matrix. It works by chunking the full
  distance matrix, and it pays for that with real accuracy costs: `nrow(spat) > big` silently
  rewrites `edge_correction` to `"none"`, the tiled branch hard-codes translation regardless of
  the requested correction, `d < r` is used where spatstat uses `d <= r`, `NA` is not removed
  from the correction sums, and spatstat's `rmax` truncation is not applied.
- **Bivariate K is duplicated** across `bi_ripleys_k()` and `bi_ripleys_k_WSI()`.
- **Bivariate NN G hand-rolls `rs`/`han`** on a full `as.matrix(dist())`, which is `O(n²)` and
  returns `NaN` where spatstat returns `0`.
- Only 1 of 16 exports has any test (4 assertions on `create_mif`).

The whole memory problem is self-inflicted. `spatstat.geom::closepairs(X, rmax)` returns only
the pairs within `rmax`, and `edge.Trans(..., paired = TRUE)` / `edge.Ripley(X, r = d)` compute
edge weights **per pair** instead of as an `n × n` matrix. Measured on this machine:

| | current approach | `closepairs` + paired weights |
|---|---|---|
| n = 200,000, small `rmax` | `n²` matrix = **319 GB** | 3.1M pairs = **75 MB**, 0.23 s |
| n = 150,000 full K | tiled, correction forced to `none` | **0.67 s**, translation *or* isotropic |

And it is **exact**, not approximate — verified against `spatstat`:

| check | max abs difference |
|---|---|
| `closepairs` + paired `edge.Trans` vs `Kest(correction="translation")` | 1.4e-20 |
| `closepairs` + `edge.Ripley` vs `Kest(correction="isotropic")` | 1.4e-20 |
| `crosspairs` + paired `edge.Trans` vs `Kcross(correction="translation")` | 5.2e-16 |
| pairs/weights computed **once** and reused across permutations vs per-permutation `Kcross` | 1.0e-17 (and **26× faster**) |
| current `bi_NN_G` `rs`/`han` vs `Gcross(correction=...)` | 0 (identical) |

That last row is the key permission slip for the refactor: **the hand-rolled bivariate G already
agrees exactly with `Gcross`**, so replacing it is numerically free. And because translation /
isotropic edge weights depend only on a pair's displacement (and the window), not on which other
points are in the pattern, the weights can be computed once over all cells and reused for every
marker subset and every permutation. That single fact is what lets the tiling machinery go away.

**Outcome:** one exact, memory-bounded engine per statistic; ~1,100 lines of dead code deleted;
consistent parameters and output columns; real test coverage; version 2.0.0.

## Ground rules

- All work on a **local-only** `feature-upgrade` branch off `feature-alex`. **Never push.**
- Commit at each numbered phase boundary (commit messages listed per phase).
- **Do not touch** `pkgdown/`, `docs/`, or `.github/` — they serve GitHub Pages. `docs/` will go
  stale (it lists `compute_metrics`); note that in NEWS and leave regeneration to the user/CI.
- User API stays put except the parameter/column renames below. `mif` object structure unchanged.
- No new hard dependencies. `spatstat.geom`, `spatstat.explore`, `spatstat.univar` are already in
  `Imports`; everything needed is there.
- **Window invariant (non-negotiable).** The observation window for a sample is always the convex
  hull of **every cell in that sample's spatial data frame** — never of a marker-positive subset,
  never of an anchor/counted subset, never of a dual-positive-filtered subset. The window and its
  area then stay fixed across all markers, all marker pairs, and all permutations within that
  sample. Audited: all 8 exported metric functions do this correctly today
  (`ripleys_k.R:104`, `bi_ripleys_k.R:97`, `bi_ripleys_k_WSI.R:96`, `nn_g.R:65`,
  `bi_nn_g2.R:74`, `pair_correlation.R:65`, `bi_pair_correlation.R:79`, `interaction.R:60`) —
  the refactor must not regress it. It is easy to break, because the natural-looking
  "build a `ppp` per marker" refactor gets it wrong. Measured cost of getting it wrong: with
  marker-positive cells confined to one corner of a sample, area goes 0.989 → 0.115 and K is
  inflated by the area ratio, **~8.5× at every radius**. Enforced by
  `test-window-invariant.R` (Phase 6).

## Environment

Already verified present in the system library (`/opt/homebrew/lib/R/4.6/site-library`, R 4.6.1):
all `DESCRIPTION` Imports/Suggests plus `covr` 3.6.5, `devtools` 2.5.2, `roxygen2` 8.1.0,
`pkgdown` 2.2.1, `SpatialExperiment` 1.22.0, `SummarizedExperiment` 1.42.0,
`VectraPolarisData` 1.16.0. **No conda env or project-local library is needed.** Defensively add
`^\.Rlib$`, `^\.conda.*$`, `^plans$` to `.Rbuildignore` and `.Rlib/`, `.conda*/` to `.gitignore`
so any later local tooling cannot leak into a commit or tarball.

---

## Phase 0 — Branch and baseline

1. `git checkout -b feature-upgrade` (from `feature-alex`).
2. Stage the 13 already-deleted working-tree files (`.Rprofile`, `R_dev/*`, `doc/*`) and commit.
   **Flag first:** `R_dev/ripleys_k_auc.R` is the only deleted file with functionality that
   exists nowhere in `R/` — it computed AUC of the K curve into `mif$derived$univariate_AUC` and
   was parked because it needed `flux::auc`. Record it as a future item (a base-R trapezoid rule
   removes the dependency); do not silently lose it.
3. Capture a numeric **baseline** from the current code before changing anything: run
   `ripleys_k`, `bi_ripleys_k`, `NN_G`, `bi_NN_G`, `pair_correlation`, `bi_pair_correlation`,
   `interaction_variable`, `dixons_s` on `example_spatial[["TMA3_[9,K].tif"]]` (the only sample
   with substantial marker counts: 536 CD3+, 83 CD8+) and save the `derived` slots to
   `tests/testthat/fixtures/baseline-v1.4.0.rds`. Every later phase is checked against this.

> Commit: `chore: remove R_dev/ and doc/ build artifacts; add v1.4.0 numeric baseline fixture`

## Phase 1 — Delete dead code

All verified by repo-wide grep across `R/`, `tests/`, `vignettes/`, `man/`, `NAMESPACE`.

| Delete | Where | Lines | Evidence |
|---|---|---|---|
| `compute_metrics()` + `man/compute_metrics.Rd` | `R/compute_metrics.R` | 259 | Not in `NAMESPACE`, no `@export`, zero call sites, `@examples` has no runnable call |
| `uni_Rip_K`, `bi_Rip_K`, `uni_NN_G`, `bi_NN_G_sample` | `R/utils-helpers.R:65,223,406,526` | 449 | Sole caller is `compute_metrics` |
| `uni_K`, `K_out`, `bi_K`, `uni_G`, `G_out`, `bi_G`, `perm_data` | `R/utils-helpers.R:44,1,187,386,371,486,34` | 153 | Called only by the four above |
| `get_exactK` | `R/utils-helpers.R:640` | 31 | Zero call sites; cannot run — line 665 uses `areapp`, defined nowhere |
| `dix_s_z`, `dix_s_c` | `R/utils-helpers.R:815,917` | 182 | Zero call sites; `dixons_s()` calls `dixon::dixon()` directly at `R/dixons_s.R:112`. `dix_s_z` also hard-codes `nsim = 1` and discards its own result |
| duplicate `get_bi_rows` | `R/bi_pair_correlation.R:3-11` | 9 | Byte-identical to `R/utils-helpers.R:24`; masked by collation order |

Survivors in `R/utils-helpers.R`: `get_bi_rows`, `list.append`, `getTile`, `calculateK` (151
lines). `getTile` and `calculateK` are then removed in Phase 2 when the tiling goes away.

Then prune `R/global.R`'s `globalVariables()` vector — `areapp`, `Naa`/`Nab`/`Nba`/`Nbb`/`Na`/`Nb`,
`'    Obs.Count'` (note the leading spaces), `un`, `W` were only masking dead code or genuine
undefined-variable bugs. Also drop the three unused `importFrom(spatstat.univar, ...)` in
`R/bi_pair_correlation.R:24` — `dkernel`/`match.kernel`/`unnormdensity` appear only in
commented-out lines.

> Commit: `refactor!: remove unexported compute_metrics and the gen-1 helper engine (~1,100 lines)`

## Phase 2 — Ripley's K: one exact engine

### New internal engine — `R/utils-k-engine.R`

```
k_pairs(pp, r_range, edge_correction, block = NULL)
  -> list(i, j, d, w)        # pair list + per-pair edge weight, computed ONCE per sample
k_from_pairs(pairs, keep, n_i, n_j, area, r_range)
  -> numeric(length(r_range))  # cumsum(whist(d[keep], breaks=r, weights=w[keep])) * area / denom
```

- `pp` is built **once per sample from all cells**, and `k_pairs` runs on that full pattern. Marker
  and marker-pair selection happens later as a logical `keep` mask over the pair list, so the
  window, its area, and the edge weights are structurally incapable of being derived from a
  subset. This is the main reason to mask rather than re-subset — it makes the window invariant
  hold by construction instead of by discipline. `area` is likewise computed once per sample and
  threaded through, never recomputed from a subset.
- `k_pairs` = `closepairs(pp, rmax = max(r_range), what = "all")`, then
  `edge.Trans(ppp(xi,yi,W), ppp(xj,yj,W), paired = TRUE)` for translation,
  `edge.Ripley(ppp(xi,yi,W), r = d)` for isotropic, `1` for none, `Kount`-based for border.
- `k_from_pairs` denominator: `n(n-1)` univariate, `n_i * n_j` bivariate — matching `Kest`/`Kcross`.
- **Must replicate spatstat's `rmax` truncation** (`Ktrans[r >= rmax.Trans(W)] <- NA`,
  `Kiso[r >= boundingradius(W)] <- NA`). The current manual paths omit this and emit numbers where
  spatstat deliberately returns `NA`.
- Use `d <= r` (via `whist` bin semantics), not the current `d < r`. Matters for integer
  HALO/Vectra coordinates.
- Do **not** drop `d == 0` pairs. The current code's `dists > 0` / `dists[dists==0] <- NA` silently
  discards genuinely co-located distinct cells that `Kest` counts.
- `block`: the large-sample flag. When `n > big`, iterate `crosspairs(pp[block], pp, rmax)` over
  index blocks and accumulate `whist` counts (additive, therefore exact). ~15 lines replacing the
  ~120 lines of nested `mclapply` tiling. **Edge correction is never downgraded.**

### `R/ripleys_k.R`

Replace all four branches (`:106`, `:193`, `:228`, plus the `:267-301` tiling) with:
compute `k_pairs` once per sample → observed K per marker via a logical `keep` mask → permuted K
by re-masking the same pair list (26× faster, exact) → `Exact CSR` = `k_from_pairs` over all cells.

Delete: `getTile` (`R/utils-helpers.R:673`), the `nrow(spat) > big → edge_correction = 'none'`
downgrade at `:97-99`, the `as.matrix(dist(...))` at `:107`, the positional column reorder
`final[,c(4,8,7,1,2,3,5,6)]` at `:187`, the two `spat[1,1] #hard coded for example` lines,
and the `gc(full=T)` at `:260`.

### `R/bi_ripleys_k.R`

Same engine via `crosspairs`. Absorbs the WSI case, so `big`/`nlarge`/`force` collapse into one
`big` argument. `Exact CSR` stays `Kest`-of-all-cells — verified valid: mean of 500 permuted
`Kcross` values agrees with `Kest(all)` to within Monte Carlo noise (<0.6%), so it is a correct
closed form for `E[Kcross]` under random labelling.

Also remove: `cat(spatial_name, "\t", combo, ...)` at `:131` and `cat(perm_n)` at `:183` (raw
console spam), the commented-out block at `:200-206`, and the fragile
`dplyr::rename('xloc' := xloc)` at `:93` (works only by data-masking fallback — breaks if a
spatial file already has an `xloc` column; use `!!xloc` as `ripleys_k.R:94` does).

### `R/bi_ripleys_k_WSI.R`

Reduce to an exported stub that `stop()`s with a pointer to `bi_ripleys_k()`. Keep
`man/bi_ripleys_k_WSI.Rd` documenting the removal. Delete `calculateK` (`R/utils-helpers.R:694`),
which also removes a live bug: its chunked branch returns `NULL` for `edge_correction =
"isotropic"` — a value its own docs advertise — yielding garbage rather than an error.

> Commit: `refactor!: single exact memory-bounded engine for Ripley's K; fold WSI into bi_ripleys_k`

## Phase 3 — Nearest Neighbour G

`R/nn_g.R` — `NN_G()` already delegates to spatstat. Two fixes: delete the `exact_G` computation
at `:67-71` (assigned, never used — pure wasted compute on every sample), and replace the
positional reorder `res[,c(7,6,4,1,2,5,3)]` at `:110` with named selection.

`R/bi_nn_g2.R` → rename file to `R/bi_nn_g.R`. Replace the entire ~150-line hand-rolled body
(`as.matrix(dist(...))` at `:111`, `bdist.points`, `km.rs`, `eroded.areas`, and three
copy-pasted permutation blocks) with `spatstat.explore::Gcross()`. Verified identical (`max|diff|
= 0` for both `rs` and `han` on `TMA3_[9,K]`), fixes the `NaN`-where-spatstat-gives-`0` edge case
on sparse markers, and drops memory from `O(n²)` to `O(n)`. Also fixes the undefined `W` at
`:170`/`:185` (currently harmless only because R never forces it — `handle.r.b.args` ignores its
window argument when `r` is supplied) and the invalid `"hans"` correction spelling.

Per requirement #10, NN G gains an **all-`NA` `Exact CSR`** column. There is deliberately no
value to put there: unlike K, `G` is intensity-dependent, so `Gest(all cells)` is *not* an
estimate of the permuted G — measured 0.52 vs 0.13 at the same radius. Document this in the Rd.

> Commit: `refactor!: delegate bivariate NN G to spatstat::Gcross; add Exact CSR column`

## Phase 4 — Consistent parameters and output columns

### Parameters

Canonical set, applied across `ripleys_k`, `bi_ripleys_k`, `NN_G`, `bi_NN_G`,
`pair_correlation`, `bi_pair_correlation`, `interaction_variable`:

`mif`, `mnames`, `r_range`, `num_permutations`, `permute`, `keep_permutation_distribution`,
`edge_correction`, `workers`, `overwrite`, `xloc`, `yloc`, `big`

| Old | Canonical | In |
|---|---|---|
| `keep_perm_dis` | `keep_permutation_distribution` | `NN_G`, `bi_NN_G` |
| `nlarge`, `force` | `big` | `bi_ripleys_k_WSI`, `bi_ripleys_k` |
| `method` | *(dropped)* | `ripleys_k` — documented "not used currently", never referenced in the body |

Per your decision, each function gains `...` and maps old names to new with a one-time
deprecation warning naming the replacement. Implement once as an internal
`resolve_deprecated_args()` in `R/utils-deprecate.R` and call it from every entry point — do not
repeat the mapping per function.

Also unify inconsistent **defaults**: `workers = 1` (was `6` in `bi_ripleys_k`/`bi_ripleys_k_WSI`
— a CRAN policy problem in examples) and `overwrite = FALSE` (was `TRUE` in `bi_ripleys_k`).
**Leave `num_permutations` defaults as they are** (50 for K/G, 100 for pcf/interaction, 1000 for
`dixons_s`) — changing them would silently alter results. Note the divergence in NEWS instead.

### Output columns

One schema for every metric:

```
<sample_id> | Marker            (univariate)
            | Anchor, Counted   (bivariate)
iter | r
Theoretical CSR | Permuted CSR | Exact CSR | Observed <K|G|g|Interaction>
Permutations Larger than Observed
Degree of Clustering Theoretical | Degree of Clustering Permutation | Degree of Clustering Exact
Run
```

Renames (`K`'s current schema is the target, so `ripleys_k`/`bi_ripleys_k` barely move):

| Function | Change |
|---|---|
| `ripleys_k` | gains `Permutations Larger than Observed` (only `bi_ripleys_k` has it today) |
| `NN_G` | `Theoretical G`→`Theoretical CSR`, `Permuted G`→`Permuted CSR`, `+Exact CSR`, `+iter`, `+Permutations Larger than Observed`, `+Degree of Clustering Exact` |
| `bi_NN_G` | `Permuted G`→`Permuted CSR`, `+Exact CSR`, `+iter`, `+Permutations Larger than Observed`, `+Degree of Clustering Exact` |
| `pair_correlation` | `Theoretical g`→`Theoretical CSR`, `Permuted g`→`Permuted CSR`, `Permuted_larger_than_Observed`→`Permutations Larger than Observed`, `Degree of Correlation *`→`Degree of Clustering *`, `+Exact CSR` |
| `bi_pair_correlation` | as above, plus `From`,`To`→`Anchor`,`Counted` |
| `interaction_variable` | `From`,`To`→`Anchor`,`Counted`, `Permuted Interaction`→`Permuted CSR`, `Degree of Interaction Permuted`→`Degree of Clustering Permutation`, `+Exact CSR` |

`Observed K` / `Observed G` / `Observed g` / `Observed Interaction` keep their statistic-specific
names — that is the one place a per-metric name is genuinely informative. Column **order** is
identical across all seven functions so `bind_rows` across metrics behaves predictably.

Implement the ordering once as an internal `standard_metric_cols()` helper rather than repeating
`relocate()` chains.

> Commit: `refactor!: unify parameter names and output column schema across all metric functions`

## Phase 5 — Real bugs found while reading

Each is independently confirmed; none is stylistic.

| File:line | Bug |
|---|---|
| `R/pair_correlation.R:125-127` | Non-overwrite path writes to `univaraite_pair_correlation` (misspelled) — appended runs land in a slot nothing reads |
| `R/bi_pair_correlation.R:166-168` | Same, `bivaraite_pair_correlation` |
| `R/interaction.R:~150` | Non-overwrite path reads `mif$derived_intraction_variable` — wrong slot *and* misspelled |
| `R/bi_pair_correlation.R:62` | `if(FALSE %in% unique(unlist(mnames)) %in% colnames(spat))` — `%in%` is left-associative, so this evaluates `(FALSE %in% mnames) %in% colnames(spat)`. The marker-existence check never fires |
| `R/bi_pair_correlation.R:95-103` | Insufficient-cell early return omits `Observed g`/`Theoretical g`, so the downstream `out$'Observed g' - out$'Theoretical g'` errors when every pair is sparse |
| `R/spatial_exp_to_mif.R:76-80` | Calls `colData()`, `spatialCoords()`, `metadata()` unqualified. These live in `SummarizedExperiment`, `SpatialExperiment`, `S4Vectors` — **none declared in `DESCRIPTION`, none imported in `NAMESPACE`**. The function fails unless the user has separately attached them. Fix: add `SpatialExperiment`/`SummarizedExperiment` to `Suggests`, guard with `requireNamespace()`, and namespace-qualify |
| `R/spatial_exp_to_mif.R:37-60` | `@examples` calls `VectraPolarisData::HumanOvarianCancerVP()` **unwrapped** — `R CMD check` will try to run it. Wrap in `\dontrun{}` |
| `R/spatial_exp_to_mif.R:72` | `class(spatial_exp) == "SpatialExperiment"` → use `inherits()` |
| `R/bi_pair_correlation.R:61`, `R/bi_ripleys_k.R:131`, `R/bi_ripleys_k_WSI.R:128` | Bare `cat()` progress spam. Route through `message()` behind a `verbose` argument, or delete |
| `R/ripleys_k.R:303` | Unsuppressed `dplyr` join message ("Joining with `by = join_by(r)`") — supply `by=` explicitly |
| **24 nested `mclapply` call sites** across all 9 metric files | **`workers` is not honored, and permuted results are irreproducible.** Only the outer per-sample loop passes `mc.cores = workers`; every inner loop over markers / marker-pairs / permutations / tiles omits `mc.cores`, so each silently falls back to `getOption("mc.cores", 2L)`. Two consequences: (a) actual process count is `workers × 2 × 2` in the nested case, so `workers = 1` is *not* serial — a CRAN 2-core policy violation; (b) forked children get non-deterministic RNG streams, so `set.seed()` does **not** make permutation results reproducible. Verified: two `ripleys_k(permute = TRUE, workers = 1)` runs under the same seed differ by up to 5909 in `Permuted CSR`, while `permute = FALSE` is bit-identical. Isolated the cause — `mclapply(mc.cores = 1)` alone *is* reproducible, so it is specifically the unbounded inner calls. **Fix:** parallelise at exactly one level (per sample, honoring `workers`), run everything inside it serially, and draw all permutation indices in the parent process under the user's seed. The Phase 2 engine makes this free — re-masking a precomputed pair list is ~26× cheaper than the per-permutation `Kcross` calls the inner parallelism existed to hide |
| `DESCRIPTION` | `parallel` and `stats` used via `::` (10× in `ripleys_k.R` alone) but absent from `Imports` |
| `inst/CITATION` | `pages = 4584-4586` is unquoted arithmetic evaluating to `-2`; `textVersion = paste()` yields `""`. `R CMD check` parses this file |
| `R/data-example-spatial.R` | Documented names (`TMA_[3,B].tiff`, `TMA_[6,F].tiff`, …) do not match reality (`TMA1_[3,B].tif`, `TMA2_[3,B].tif`, `TMA3_[7,B].tif`, `TMA3_[9,K].tif`, `TMA3_[8,U].tif`). `TMA_[6,F]` does not exist |

> Commit: `fix: correct misspelled derived slots, undeclared imports, and broken guards`

## Phase 6 — Test coverage

Create `tests/testthat/helper-mif.R` owning the shared fixture (currently duplicated at the top
of `test_creatMIF.R`): a small `mif` built from `example_spatial[["TMA3_[9,K].tif"]]` thinned to
~300 cells to keep the suite fast, plus `mnames_good`.

Rename `test_creatMIF.R` → `test-create-mif.R` (testthat 3 convention) and broaden it: the
`patient_id`/`sample_id`/`derived` slots, the 10 `stopifnot` messages, `print.mif`.

New files, one per exported function:

- `test-ripleys-k.R`, `test-bi-ripleys-k.R` — **the important ones.** Assert the engine equals
  `spatstat` directly: for each of `translation`/`isotropic`/`none`/`border`, and for both a
  convex-hull window and a rectangle, `expect_equal(observed, Kest(...)$<corr>, tolerance = 1e-12)`.
  Assert weight-reuse permutations equal per-permutation `Kest`/`Kcross`. Assert `big` blocking
  gives bit-identical results to unblocked. Assert `edge_correction` is *never* silently changed.
- `test-nn-g.R`, `test-bi-nn-g.R` — equality with `Gest`/`Gcross` for `rs`/`km`/`han`; the sparse
  marker case that currently returns `NaN` returns `0`.
- `test-pair-correlation.R`, `test-interaction-variable.R`, `test-dixons-s.R`,
  `test-subset-mif.R`, `test-merge-mifs.R`, `test-marker-freq-diff.R`, `test-plot-immunoflo.R`
  (returns a ggplot; check layers not pixels), `test-spatial-exp-to-mif.R`
  (`skip_if_not_installed("SpatialExperiment")`, build a tiny SPE by hand — do **not** pull
  `VectraPolarisData` into the suite).
- `test-window-invariant.R` — enforces the Ground-rules window invariant for **all 8** metric
  functions. Fixture: a sample whose cells fill a large region but whose marker-positive cells are
  confined to one corner. For each function, assert the result matches a `spatstat` reference
  computed with the **full-sample** window, and assert it does *not* match the subset-window
  reference (which differs by the ~8.5× area ratio, far outside any tolerance). Also assert that
  the window/area used for `Observed`, `Permuted CSR`, and `Exact CSR` within one sample are
  identical, and that adding an extra marker column to the input does not change any existing
  marker's result.
- `test-output-schema.R` — table-driven: for every metric function, assert exact column names
  **and order** against the Phase 4 schema. This is what keeps the schema from drifting again.
- `test-deprecated-args.R` — each old name warns once and still produces the new behaviour.
- `test-regression-vs-baseline.R` — compare against `baseline-v1.4.0.rds` from Phase 0.
  Deliberate changes get an explicit documented tolerance or an `expect_failure` with a comment
  naming the fix; everything else must match.

Guard every test with `workers = 1` (CRAN limits to 2 cores) and small `r_range`/`num_permutations`.

Target: `covr::package_coverage()` above 70%, all 16 exports exercised.

> Commit: `test: add helper fixture and coverage for all 16 exported functions`

## Phase 7 — Metadata, docs, NEWS

- `DESCRIPTION`: `Version: 2.0.0`; add `parallel`, `stats` to `Imports`; add `SpatialExperiment`,
  `SummarizedExperiment`, `covr` to `Suggests`; add `Collate:` so the collation order that makes
  `get_bi_rows` resolve is explicit rather than accidental.
- `.Rbuildignore`: add `^vignettes/spatialexperiment\.Rmd$` (requirement #6 — pkgdown builds it
  for GitHub Pages from `vignettes/`, `.Rbuildignore` only affects the CRAN tarball, so this is
  exactly the right lever). Add `^plans$`, `^\.Rlib$`, `^\.conda.*$`. Drop the stale entries
  whose targets no longer exist (`^renv$`, `^R/merge_for_jordan\.R$`, `^R/simulate_mif\.R$`,
  `^data/example_tma\.rda$`, `^bugs`, `^NAMESPACE_bkup`, `^\.ccache$`, `^_notes$`).
- `.gitignore`: add `.Rlib/`, `.conda*/`. Note that `dev/`, `/doc/`, `/Meta/` are listed but
  already tracked, so those rules do nothing today.
- Regenerate `NAMESPACE` and `man/` with `roxygen2::roxygenise()`. Delete
  `man/compute_metrics.Rd`.
- Update `vignettes/intro.Rmd` and `vignettes/deriving_functions.Rmd` for the renamed
  arguments/columns (`intro.Rmd:284,305` use `keep_perm_dis`; `deriving_functions.Rmd:138,158`
  read pcf columns). `vignettes/spatialexperiment.Rmd:180,213` call `ripleys_k`/`NN_G`.
- **Create `NEWS.md`** (none exists anywhere in the repo). Under `# spatialTIME 2.0.0`, sections:
  Breaking changes (removed `compute_metrics`; `bi_ripleys_k_WSI` now errors; column renames;
  `method` dropped) / Deprecations (old argument names, with the replacement for each) /
  Refactor (single exact K engine; `Gcross` for bivariate G; edge corrections no longer silently
  downgraded — this **changes results** for any sample over `big`, for the better) / Bug fixes
  (the Phase 5 table) / Internal (~1,100 lines removed; test coverage). Mention that `docs/`
  needs a `pkgdown::build_site()` refresh.

> Commit: `docs: bump to 2.0.0, add NEWS.md, update vignettes and build config`

---

## Verification

Run in order; each must pass before the next phase's commit.

1. **Engine correctness (the whole point).** Standalone script, not just the suite:
   ```r
   devtools::load_all(".")
   # for each correction x window x marker density:
   #   engine value == Kest/Kcross/Gest/Gcross value, tolerance 1e-12
   ```
   Floating-point-only deviation is the acceptance bar, per your requirement. Measured today the
   engine lands at 1e-20–1e-16.
2. **Regression vs baseline.** `test-regression-vs-baseline.R` against
   `tests/testthat/fixtures/baseline-v1.4.0.rds` (17 entries, captured in Phase 0 from the
   unmodified v1.4.0 code). Compare **exactly** only the deterministic columns — `Observed K` /
   `Observed G` / `Observed g` / `Observed Interaction`, `Theoretical CSR`, `Exact CSR` — because
   the v1.4.0 permutation paths are irreproducible (see the `mclapply` row in Phase 5), so the
   `Permuted CSR` baseline values are not a stable reference. Check permuted columns loosely
   instead (same order of magnitude, monotone in `r`), and add a *new* reproducibility assertion
   that post-refactor permutations **are** bit-identical under a fixed seed — something v1.4.0
   cannot satisfy. Every exact-column difference must be an intended fix with a NEWS entry — in
   particular, samples over `big` now keep their requested edge correction, so their K values
   legitimately change.
3. **Memory/scale.** Simulate 200,000 cells; confirm `bi_ripleys_k` completes with
   `edge_correction = "translation"` under a few hundred MB, and that `big` blocking returns
   bit-identical numbers to unblocked.
4. **Suite + coverage.** `devtools::test()` clean; `covr::package_coverage()` > 70%.
5. **Package check.** `devtools::check(document = TRUE)` — target 0 errors / 0 warnings, and
   confirm no note about undeclared imports, orphan Rd, or examples requiring
   `VectraPolarisData`.
6. **Vignettes build.** `devtools::build_vignettes()` — confirms the renames propagated. Then
   `R CMD build .` and check the tarball does **not** contain
   `vignettes/spatialexperiment.Rmd`.
7. **SpatialExperiment path.** With `VectraPolarisData` (installed, but never a dependency),
   run `spatial_exp_to_mif()` on `HumanOvarianCancerVP()` once by hand to confirm the
   `requireNamespace` guards and namespace qualification work end to end.

## Designed for what comes next

Deliberate choices that make the planned adaptability/integration work cheap:

- `k_pairs()` / `k_from_pairs()` split means a new count-based statistic (`Kinhom`, `Lest`, cross-K
  variants, marked K) is a new denominator over the same pair list — no new distance or
  edge-correction code.
- `standard_metric_cols()` plus `test-output-schema.R` means adding a metric means opting into an
  existing schema, so cross-metric `bind_rows` and downstream tooling keep working.
- `resolve_deprecated_args()` is the reusable seam for the next rename round.
- Explicit `Collate:` removes the accidental-collation-order hazard that currently makes one
  `get_bi_rows` definition silently shadow another.
- `Exact CSR` exists as a column for every metric even where it is `NA`, so filling it in later
  (pcf has the same random-labelling closed form as K) is additive, not breaking.
