# spatialTIME 2.0.0 — Cleanup & Refactor

> **Status: Phases 0–4 are implemented and committed on local branch `feature-upgrade`.**
> Phases 5–7 remain. Jump to [Remaining work](#remaining-work). Sections above that are
> kept as the record of what was done and what the measurements actually showed.
>
> | Commit | Phase | Net |
> |---|---|---|
> | `ec510b0` | 0 — branch, artifact removal, v1.4.0 numeric baseline | — |
> | `160918f` | 1 — delete dead code | **−1,231 lines** |
> | `5129551` | 2 — exact memory-bounded Ripley's K engine, WSI folded in | +1,056 / −952 |
> | `1b4a94c` | 3 — bivariate NN G delegated to `Gcross` | +630 / −347 |
> | `2f3265e` | 4 — unified parameters and output columns | +888 / −398 |
>
> Cumulative vs `feature-alex`: **39 files, +3,070 / −6,312 lines.** Test suite: 195
> assertions, all passing (was 4). Nothing pushed.

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
  the requested correction, bin edges do not match spatstat's on tied distances, `NA` is not
  removed from the correction sums, and spatstat's `rmax` truncation is not applied.
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
  sample. Audited at the start: all 8 exported metric functions already did this correctly, so the
  job was to not regress it. It is easy to break, because the natural-looking
  "build a `ppp` per marker" refactor gets it wrong. Measured cost of getting it wrong: with
  marker-positive cells confined to one corner of a sample, area goes 0.989 → 0.115 and K is
  inflated by the area ratio, **~8.4× at every radius**. Now structurally guaranteed rather than
  merely intended: `k_pairs()` builds the pair list once from the full-sample pattern and marker
  selection is a logical mask over it, so a subset can no longer reach the window at all. Enforced
  by `tests/testthat/test-window-invariant.R` (11 assertions, already written), which also checks
  that a marker's result is unchanged by which *other* markers were requested.

## Environment

Already verified present in the system library (`/opt/homebrew/lib/R/4.6/site-library`, R 4.6.1):
all `DESCRIPTION` Imports/Suggests plus `covr` 3.6.5, `devtools` 2.5.2, `roxygen2` 8.1.0,
`pkgdown` 2.2.1, `SpatialExperiment` 1.22.0, `SummarizedExperiment` 1.42.0,
`VectraPolarisData` 1.16.0. **No conda env or project-local library is needed.** Defensively add
`^\.Rlib$`, `^\.conda.*$`, `^plans$` to `.Rbuildignore` and `.Rlib/`, `.conda*/` to `.gitignore`
so any later local tooling cannot leak into a commit or tarball.

---

## Phase 0 — Branch and baseline  ✅ DONE (`ec510b0`)

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

**Outcome:** baseline captured for all 8 metric functions plus `subset_mif` (17 entries). It
carries a `deterministic_cols` attribute, because the v1.4.0 permutation paths turned out to be
irreproducible (see the `mclapply` row in Phase 5) — so only `Observed *`, `Theoretical CSR` and
`Exact CSR` are a stable reference. `plans/future-ripleys-k-auc.md` records the deferred
`ripleys_k_auc()` capability.

## Phase 1 — Delete dead code  ✅ DONE (`160918f`)

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

**Outcome:** −1,231 lines. `utils-helpers.R` went 998 → 26 lines. Verified behaviour-neutral:
all 15 exports still resolve, all 15 deleted names are gone, and `ripleys_k(permute = FALSE)`
output was identical to the baseline across every deterministic column.

## Phase 2 — Ripley's K: one exact engine  ✅ DONE (`5129551`)

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
- Bin edges must come from `handle.r.b.args()` + `whist()`, not a hand-rolled comparison.
  **Correction to this plan's original claim:** I had written that spatstat uses `d <= r` where
  the old code used `d < r`. That is not right, and the truth matters. `Kest` has *two* internal
  paths. Ask only for `"none"` or `"border"` with evenly spaced `r` and it takes a fast C path
  binning right-closed (`d <= r`); ask for `"none"` alongside any other correction, or ask for
  translation/isotropic, and it routes through `whist()`, binning left-closed. On an integer grid
  those disagree — 390.48 vs 0 at the same radius. The engine therefore always uses the `whist`
  path, which makes it bit-identical to `Kest(correction = c("none", "translation"))$un`
  (verified `max|diff| = 0`) and, unlike the fast path, internally consistent with its own
  translation output. Only visible on tied distances, i.e. exactly the integer/half-integer
  coordinates HALO and Vectra emit.
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

**Outcome, measured:** 200,000 cells with translation correction in **1.1 s / 738 MB peak**,
against 320 GB for the n² matrix the tiling existed to avoid. Agreement with spatstat is 0 in
most combinations, worst 9.1e-13 absolute on K ≈ 11,000 (8e-17 relative). `big = 100` and
`big = 1e9` give identical numbers. Permutations are now reproducible under `set.seed()` **and**
independent of `workers`.

One residual, documented in the engine header and bounded by a test: bivariate **isotropic** can
differ from `Kcross` by ~1e-6 relative at a single radius. Cause is a genuine geometric
degeneracy — `Kest` uses `closepairs()`, `Kmulti` uses `crosspairs()`, they differ by 1 ulp
(1.4e-17) in a pair distance, and the isotropic weight jumps discontinuously when the circle of
radius `d` passes exactly through a window vertex. Since the window is the convex hull of the
cells, its vertices *are* cells, so such pairs occur. A 2e6-point Monte Carlo integration of the
true arc fraction confirmed spatstat's value is the accurate one and the discrepancy is confined
to that one pair. Translation — the default — is exact.

## Phase 3 — Nearest Neighbour G  ✅ DONE (`1b4a94c`)

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

**Outcome:** ~150 hand-rolled lines replaced by one `Gcross()` call, and it was numerically
free — `Observed G` and `Theoretical CSR` are **identical** to the v1.4.0 hand-rolled values
(`max|diff| = 0`), and both agree with `Gest`/`Gcross` to 0 for `rs`, `km`, `han` and `none`.
Also fixed: `NaN` on sparse markers, the undefined `W`, the invalid `"hans"` spelling, and
selecting spatstat's estimate by column *position* — which silently returns `rs` when you ask
for `km`, because `Gcross(correction = "km")` puts `rs` in column 3.

## Phase 4 — Consistent parameters and output columns  ✅ DONE (`2f3265e`)

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

**Outcome:** all seven metrics now emit identical columns in identical order, differing only in
`Observed K`/`G`/`g`/`Interaction`; `bind_rows()` across all seven produces a clean frame. Old
argument names warn and forward; unknown names in `...` are a hard error, so adding `...` did
not open a typo hole. Four Phase 5 bugs fell out here as side effects (misspelled derived slots,
the `%in%` precedence guard, the sparse-pair early return, and the `cat()` spam) — struck from
Phase 5 below.

<a name="remaining-work"></a>
# Remaining work

Everything below is still to do. Phases 0–4 above are committed.

## Phase 5 — Remaining real bugs

This grew a lot. A follow-up audit of the eight exports the refactor had *not* touched turned up
considerably more than the original pass, and the worst of it is statistical rather than
structural. Everything below marked ✔ I re-verified myself by running the code — I did not take
the audit's word for it, and in one case (below) the audit's claim was overstated.

Four items originally in this phase were fixed as side effects of Phase 4 and are struck: the three
misspelled derived slots, the `%in%`-precedence marker guard, the sparse-pair early return, and the
`cat()` spam.

### Tier 1 — produces wrong numbers or wrong plots

| File | Bug |
|---|---|
| `R/marker_freq_diff.R:96-103` | ✔ **Every p-value this function has ever produced is wrong.** The Fisher table is built from `contains(c(marker, "Total"))`, so its two rows are *positives* and *totals* — not positives and negatives. The totals row double-counts the positives in the margin. Verified on `example_spatial[[1]]` / `CD3..Opal.570..Positive`: the function returns **7.949897e-07**, which is exactly `fisher.test(rbind(pos, total))`; the correct `fisher.test(rbind(pos, neg))` is **6.001431e-07**. Fix: second row must be `Total - pos`. |
| `R/marker_freq_diff.R:98` | ✔ **`contains(marker)` is a substring match, so prefix-colliding markers are silently merged into one r×c test.** With `mnames = c("CD3..CD8.", "CD3..CD8..FOXP3.")` — both real columns of the shipped data — the selection for `"CD3..CD8."` pulls in `CD3..CD8..FOXP3.` as well and hands `fisher.test` a **3×2** table, stored as if it were that marker's 2×2 p-value. *Audit nuance corrected:* the audit claimed this changes the answer; on the shipped data both come out at p = 1 because the counts are degenerate (0/3 and 0/0), so it does **not** reproduce numerically there. The 3×2 table is real and confirmed, and it will change answers on non-degenerate counts. Fix: exact name matching, not `contains()`. |
| `R/plot_immunoflo.R:126-129` | ✔ **`cell_type` has no effect on the plot.** `ggplot2::aes(..., shape = cell_type)` maps the *string* `"Classifier.Label"`, not the column. Verified: `p$mapping$shape` is `~cell_type` and exactly **one** shape (3) is drawn for every cell despite two levels in the data — so the documented example renders Tumor and Stroma indistinguishably, with a one-entry legend labelled with the column name. Fix: `shape = .data[[cell_type]]`. Then `scale_shape_manual(values = c(3, 16))` at `:141` hard-codes two shapes, so a classifier with >2 levels will error at draw time — widen it at the same time. |

### Tier 2 — silent data loss or hard failure on the default path

| File | Bug |
|---|---|
| `R/marker_freq_diff.R:126-139` | ✔ **The `overwrite` branches are inverted.** With `overwrite = FALSE` on a mif that already has results, the previous run is **destroyed** (verified: `nrow` stays 1, `Run` stays 1). On a *fresh* mif it takes the other branch and ships `Run = -Inf` (verified) — `max()` of an empty slot, and it takes `max()` of the whole data frame rather than of `$Run`. `overwrite = FALSE` is the default, so both paths are the common case. |
| `R/dixons_s.R:158,167` | ✔ **`overwrite = FALSE` — the default — always errors.** The append path does `mutate(Run = max(Run) + 1)` where `Run` is expected to be a column of the freshly computed table, but `Run` is only ever created in the *overwrite* branch. `max(Run)` resolves lexically and fails. `'Run'` is in `global.R`'s `globalVariables()` list, which is precisely the "a stale entry hides a genuine undefined-variable bug" failure mode that file's own comment warns about. |
| `R/subset_mif.R:41-59` | **Two failure modes from one line.** `out` is assigned only inside `if(nrow(tmp) > 2)` while `rbind.data.frame(summary, t(out))` runs unconditionally. First sample with ≤2 cells at the level → `object 'out' not found` (✔ reproduced). A *later* sample → silently re-appends the **previous** sample's row, so the summary gains a row matching no retained spatial frame. The silent case is the dangerous one. |
| `R/subset_mif.R:45-60` | ✔ **The whole new `sample` table comes back character.** `c(patient, id, unlist(counts), unlist(percent))` coerces at `c()`, so counts and proportions ship as `"3611"`, `"0.00249238438105788"`. Build the row as a typed one-row data frame instead of a `c()` + `t()`. Separately, the `% ` columns hold *proportions* (`sum/nrow`) while `marker_freq_diff` puts true *percentages* in its `%` columns. **Decided: percentage everywhere** — multiply `subset_mif`'s values by 100 so `%` matches its own label and the other export. This changes `subset_mif` output by 100×, so it needs a loud NEWS entry and an explicit assertion in `test-subset-mif.R`. |
| `R/plot_immunoflo.R:155,161` | ✔ **`dev.off()` twice.** `on.exit(dev.off())` plus an explicit `grDevices::dev.off()`. With no other device open the function **errors on exit** after correctly writing the pdf; with a user device open it silently closes the caller's device. Keep one, not both. |
| `R/dixons_s.R:102-108` | **`tablaC`'s two code paths have colliding schemas.** `dixon::dixon()`'s real column names carry padding spaces (`"  df "`, `"  P.rand"`) and only `tablaZ` is de-spaced at `:114`; the early return builds clean names, so `bind_rows` yields *both* `"  P.rand"` and `"P.rand"`. Worse, `colnames(final_df_c)[1] = mif$sample_id` renames the **`df`** column, which `:130` then overwrites with the sample name — so degrees-of-freedom is destroyed rather than reported. |
| `R/create_mif.R:58-69` vs `:85-90` | ✔ **`sample_data_clean` and `clinical_data_clean` are computed and thrown away.** The returned mif uses the raw inputs, so the documented `sample_string` feature does not exist in the product, the work is wasted on every call, and — worst — the discarded `full_join` can *fail* on incompatible id types, making `create_mif` refuse to build a mif because of a join whose result it never uses. That is almost certainly why every example in the package carries `mutate(deidentified_id = as.character(...))`. **Decided: delete the block** (~30 lines) and correct `@return`, which documents a `sample_string` that was never populated. Side benefit: integer-vs-character ids stop erroring, so the `as.character()` boilerplate can come out of the examples and vignettes. |
| `R/print.R:10-12` | ✔ **A partially built mif cannot be printed.** `NULL[["id"]]` errors in R, so `clinical = NULL`, `sample = NULL` or `patient_id = NULL` each make auto-printing throw (`subscript out of bounds`) — the object cannot even be echoed at the console. Also returns `NULL` instead of `invisible(x)`, so `y <- print(mif)` yields NULL and print cannot sit in a pipe. |
| `R/merge_mifs.R:36,91` | **An empty list passes the guard and dies obscurely.** `if(is.null(mifs) | length(mifs) == 1)` misses `length == 0`, then `names(sizes) = seq(length(sizes))` hits `seq(0)`, which is `c(1, 0)` — length 2 — giving `'names' attribute [2] must be the same length as the vector [0]`. Use `length(mifs) < 2` and `seq_along`. |
| `R/merge_mifs.R:40-52,72-74` | With `check.names = FALSE`, differing `patient_id`s are stored as a length-2 vector rather than resolved, and duplicate spatial names survive — which the documented example `merge_mifs(list(x, x), check.names = FALSE)` produces. Derived slots are also bound with no `Run` renumbering, so two mifs each holding `Run = 1` yield duplicate `(sample_id, Run)` keys and break every metric's `max(Run) + 1`. |

### Tier 3 — declarations, docs, and check cleanliness

| Item | Fix |
|---|---|
| ✔ `Roxygen: list(markdown = TRUE)` is **absent** from `DESCRIPTION` | **This one affects work already done.** Without it roxygen does not process markdown, so 11 existing man pages ship literal backticks and dead `[spatialTIME::create_mif()]` links — and the new Rd blocks written in Phases 2–4 lean on markdown heavily (`[spatstat.explore::Kest()]`, `**bold**`), so they would ship broken too. Add the field *before* the Phase 7 roxygenise. |
| ✔ `utils` and `methods` used but undeclared | `utils::globalVariables`/`packageVersion`; `is()` at `dixons_s.R:58`, `marker_freq_diff.R:36`, `plot_immunoflo.R:58`. `methods` is not attached under `R --vanilla`. Replace `is()` with `inherits()` and drop `methods` entirely; declare `utils`. |
| `R/spatial_exp_to_mif.R:76-80,37-60,72` | Unqualified `colData()`/`spatialCoords()`/`metadata()` from three undeclared Bioconductor packages; `@examples` that calls `installed.packages()`, `BiocManager::install()` and downloads `VectraPolarisData` (an outright CRAN policy violation, not just a check failure); `class(x) == "SpatialExperiment"` → `inherits()`. |
| `R/dixons_s.R:99` | ✔ `` `Image Location` = spat$`Image Location`[1] `` — that column is not in the shipped data, so the expression is `NULL` and `mutate()` silently *removes* rather than adds. On real Vectra data it appears only on early-return rows, making the table ragged. |
| `R/dixons_s.R:86` | Nested `mclapply` with no `mc.cores` — the same bug fixed in the metrics. `workers` cannot control it, and when it forks, errors surface as `$ operator is invalid for atomic vectors` instead of the real message. |
| `R/dixons_s.R:96` | `nrow(df_tab) == 0 | nrow(df_tab) == 1` is dead: `Marker` is a factor with both levels, so `table()` always has 2 rows. Only the `< 3` clause fires. `type` is also never validated, so `type = "z"` runs the whole permutation loop and returns unchanged. |
| `inst/CITATION` | `pages = 4584-4586` unquoted arithmetic → `-2`; `textVersion = paste()` → `""`; title typo `"immunofluorescnece"` — the string users paste into papers; no `doi` field. |
| `R/data-example-*.R` | `example_spatial` names are all wrong (actual: `TMA1_[3,B].tif`, `TMA2_[3,B].tif`, `TMA3_[7,B].tif`, `TMA3_[9,K].tif`, `TMA3_[8,U].tif`); `example_clinical` documents `deidenitifed_sample`; `example_summary` ships a literal `...` item. Document that the clinical/summary ids are **integer** while all spatial ids are character — that mismatch is why every example needs an `as.character()` call. |
| `DESCRIPTION` | Stale `Packaged:` field checked into source; `Description:` still describes the package as Ripley's K only and omits G, pcf, Dixon and interaction. |

Given Tier 1, `marker_freq_diff` and `plot_immunoflo` are effectively broken as shipped, so this
phase is no longer optional polish — split it into two commits, Tier 1 separately, so the
statistical fix is reviewable on its own.

> Commits: `fix!: marker_freq_diff Fisher table and plot_immunoflo cell_type mapping`
> then `fix: Run bookkeeping, device handling, undeclared imports, docs`

### Explicitly *not* a regression

The audit flagged that dropping `furrr`/`future` removes the Windows parallel path. Checked against
`feature-alex`: those packages appear **zero** times in the exported metric functions and only in
the unexported `compute_metrics()`. The exported metrics have always used `parallel::mclapply`,
which simply runs serially on Windows. So this is a **pre-existing** limitation, not something the
refactor introduced, and removing the two unused Imports is safe. Whether to add a Windows-capable
backend is a separate future decision — noted, not actioned.
## Phase 6 — Test coverage

Current state: **7 of 16 exports are exercised** (195 assertions across 5 files). Untested:
`bi_ripleys_k_WSI`, `dixons_s`, `marker_freq_diff`, `merge_mifs`, `plot_immunoflo`,
`spatial_exp_to_mif`, `subset_mif`, `print.mif` — and `create_mif` has only the original four weak
assertions.

1. **`tests/testthat/helper-mif.R`** — own the shared fixture that `test_creatMIF.R` and the four
   new test files currently each rebuild. Thin `example_spatial[["TMA3_[9,K].tif"]]` to a few
   hundred cells so the suite stays fast, and expose the corner-marker fixture from
   `test-window-invariant.R` for reuse.
2. **Rename `test_creatMIF.R` → `test-create-mif.R`** (testthat 3 convention; note the missing "e"
   in the current name) and broaden it: the `patient_id`/`sample_id`/`derived` slots, each of the
   ten `stopifnot()` messages, and `print.mif`. The existing assertions are near-tautological —
   `mif$sample` is compared to the unmodified input even though `create_mif` joins and reorders.
3. **One file per remaining export**, each asserting behaviour rather than just "does not error":
   - `test-subset-mif.R` — must include the regression for the `out` bug above: a first sample with
     ≤2 matching cells, and a *middle* sample with ≤2, asserting `nrow(sample)` equals
     `length(spatial)` and that no row is duplicated.
   - `test-dixons-s.R` — a marker pair that hits the early return alongside one that does not, so
     the `Image Location` path is actually executed.
   - `test-plot-immunoflo.R` — returns a `mif` whose `derived$spatial_plots` holds ggplot objects;
     assert layers/mappings, not pixels. Assert that calling it with `filename = NULL` leaves
     `grDevices::dev.list()` unchanged (the double-`dev.off()` regression).
   - `test-merge-mifs.R` — derived slots with *differing* column sets, plus the `check.names`
     guards.
   - `test-marker-freq-diff.R`, `test-spatial-exp-to-mif.R`
     (`skip_if_not_installed("SpatialExperiment")`, hand-build a tiny SPE — do **not** pull
     `VectraPolarisData` into the suite), `test-bi-ripleys-k-wsi.R` (asserts the stub errors with a
     pointer to `bi_ripleys_k()`).
4. **Two cross-cutting contract tests** — each catches several Phase 5 bugs at once, so write these
   before the per-function files:
   - `test-run-semantics.R` — for *every* function that writes a derived slot, call it with
     `overwrite = TRUE`, then twice with `overwrite = FALSE`, asserting `Run` is `1`, then
     `c(1, 2)`, with `expect_no_warning` and row counts doubling. `dixons_s` and
     `marker_freq_diff` fail this today in two different ways; the seven refactored metrics already
     pass it via `write_derived()`, which is the implementation the other two should adopt.
   - Extend `test-output-schema.R` to `Dixon_Z`, `Dixon_C` and `frequency_difference`: assert no
     column name has leading/trailing whitespace, none is duplicated, and the column set is
     invariant to whether any sample or marker pair hit a degenerate path. That single assertion
     covers the `tablaC` schema collision, the `Image Location` raggedness, and the
     `"  P.rand"`/`"P.rand"` split.
5. **A golden-value test for `marker_freq_diff`** — the highest-value single assertion in the whole
   phase. Assert its p-value equals `fisher.test()` on the hand-built positives-vs-negatives table.
   It fails today (7.949897e-07 produced vs 6.001431e-07 correct), and it is the only thing that
   will keep the Tier 1 fix honest.
6. **`test-regression-vs-baseline.R`** against `tests/testthat/fixtures/baseline-v1.4.0.rds`.
   Compare **exactly** only the columns named in the fixture's `deterministic_cols` attribute,
   because the v1.4.0 permutation paths are irreproducible. Two differences are expected and must
   be asserted *as* differences with a comment naming the cause: `edge_correction = "none"` now
   uses the `whist` binning (the example data has half-integer centres, so ties are common), and
   any sample over `big` now keeps its requested correction. Add the converse assertion too — that
   post-refactor permutations **are** bit-identical under a fixed seed and across `workers` values,
   which v1.4.0 could not satisfy.

Guard every test with `workers = 1` and small `r_range`/`num_permutations` (CRAN allows two cores).
Target `covr::package_coverage()` > 70% with all 16 exports exercised.

> Commit: `test: helper fixture plus coverage for all 16 exported functions`

## Phase 7 — Metadata, docs, NEWS

- **`DESCRIPTION`**: `Version: 2.0.0`. **Add `Roxygen: list(markdown = TRUE)` first** — it is
  absent today, which is why 11 man pages ship literal backticks, and without it the markdown-heavy
  Rd blocks written in Phases 2–4 would ship broken too. Add `parallel`, `stats`, `utils` to
  `Imports` (not `methods` — replace the three `is()` calls with `inherits()` instead). Add
  `SpatialExperiment`, `SummarizedExperiment`, `covr`, `spatstat.random` to `Suggests`
  (`spatstat.random` is used by the new tests). Add an explicit `Collate:`.
  **Drop four now-unused Imports** — `purrr`, `furrr`, `future` and `tidyselect` have **zero**
  remaining `::` uses anywhere in `R/`; they existed only for `compute_metrics()` and the gen-1
  helpers deleted in Phase 1. Removing `furrr`/`future` in particular drops a substantial
  dependency subtree. Verify with a fresh `R CMD check` rather than by grep alone.
- **`.Rbuildignore`**: add `^vignettes/spatialexperiment\.Rmd$` (requirement #6 — pkgdown builds
  articles from `vignettes/` for GitHub Pages, and `.Rbuildignore` only affects the CRAN tarball,
  so this is exactly the right lever). Add `^plans$`, `^\.Rlib$`, `^\.conda.*$`. Drop the stale
  entries whose targets no longer exist (`^renv$`, `^R/merge_for_jordan\.R$`,
  `^R/simulate_mif\.R$`, `^data/example_tma\.rda$`, `^bugs`, `^NAMESPACE_bkup`, `^\.ccache$`,
  `^_notes$`, `^R_dev$`, `^doc$`).
- **`.gitignore`**: add `.Rlib/`, `.conda*/`.
- **Regenerate `NAMESPACE` and `man/`** with `roxygen2::roxygenise()`. Expect a large diff:
  `RoxygenNote` is 7.3.2 but the installed roxygen2 is 8.1.0, so every `.Rd` gets rewritten. Do it
  as its own commit so the mechanical churn does not obscure the hand-written changes. `man/` has
  no orphans left (`compute_metrics.Rd` was deleted in Phase 1); `bi_nn_g2.R` was renamed but its
  Rd is `bi_NN_G.Rd`, which is still correct.
- **Vignettes** — these will now error or mislead, so all three need updating:
  - `intro.Rmd:182,204` pass the defunct `method = "K"`; `:278` documents `keep_perm_dis`.
  - `deriving_functions.Rmd:90,110` pass `keep_perm_dis`; `:140,160` plot
    `Degree of Correlation Permuted`; `:161,180,199` `facet_grid(~From)`; `:196` filters `From != To`
    — all renamed. Note this vignette is `.Rbuildignore`d but **is** published by pkgdown, so it
    still has to work.
  - `spatialexperiment.Rmd:219` passes `keep_perm_dis`.
- **Create `NEWS.md`** (none exists anywhere in the repo). Sections under `# spatialTIME 2.0.0`:
  *Breaking changes* — `compute_metrics()` removed; `bi_ripleys_k_WSI()` now errors with a pointer;
  the column renames per metric; `method` and `force` defunct; `workers` and `overwrite` defaults
  changed; `subset_mif()`'s `%` columns are now true percentages (×100 from before) and its summary
  table is properly typed rather than all-character; `create_mif()` no longer computes the unused
  `sample_string`, so it no longer rejects mismatched integer/character ids. *Deprecations* — each old argument name with its replacement and the 3.0.0 removal note.
  *Behaviour changes that alter numbers* — `edge_correction` is no longer silently downgraded
  (results change, for the better, on any sample over `big`); `"none"` now uses spatstat's `whist`
  binning, which differs from `Kest(correction = "none")` alone on tied coordinates; permutations
  are now reproducible and `workers`-independent, so a fixed seed gives different numbers than
  v1.4.0 did. *Bug fixes* — the Phase 5 tiers plus the four already fixed. Lead with the two that invalidate
  prior output rather than burying them: **`marker_freq_diff()`'s Fisher p-values were all wrong**
  (the contingency table used totals where it needed negatives — 7.949897e-07 where the correct
  answer is 6.001431e-07 on the shipped example), so any analysis relying on them must be re-run;
  and **`plot_immunoflo()`'s `cell_type` argument never had any effect**, so every plot produced
  with it drew one shape for all cells. *Internal* — the engine,
  the dependency reduction, coverage. Note that `docs/` needs a `pkgdown::build_site()` refresh
  (deliberately not done here, since `docs/` is off-limits).

> Commit: `docs: bump to 2.0.0, add NEWS.md, refresh vignettes and build config`

## Open questions for a maintainer

Not blocking, but neither should be decided by a refactor:

1. **`interaction_variable()` normalisation.** The numerator counts anchor cells but the
   denominator is anchor **+** counted, so `Observed Interaction` tops out near
   `100 * n_anchor / (n_anchor + n_counted)` rather than 100, and its ceiling moves with the
   relative abundance of the two markers. Preserved exactly as-is and documented in its Rd,
   including how to rescale. Confirm against Steinhart et al. whether that is intended.
2. **`ripleys_k_auc()`.** Recorded in `plans/future-ripleys-k-auc.md`. Worth reviving as a
   schema-generic `metric_auc(mif, slot)` now that all seven metrics share one output shape — and
   without the `flux` dependency, since it is a three-line trapezoid rule.
## Verification

Run in order; each must pass before the next phase's commit.

1. ✅ **Engine correctness (the whole point).** Now lives in `tests/testthat/test-k-engine.R`
   (38 assertions) and `test-nn-g.R` (41): engine value == `Kest`/`Kcross`/`Gest`/`Gcross` at
   `tolerance = 1e-12`, for every correction, univariate and bivariate, on both continuous and
   integer coordinates, plus chunked-vs-unchunked and permutation-reuse equivalence.
   Floating-point-only deviation was the acceptance bar; achieved 0 in most combinations, worst
   9.1e-13 absolute on K ≈ 11,000. The one documented exception (bivariate isotropic at a
   circle-through-vertex degeneracy) is bounded by its own test.
2. **Regression vs baseline** (still to write). `test-regression-vs-baseline.R` against
   `tests/testthat/fixtures/baseline-v1.4.0.rds` (17 entries, captured in Phase 0 from the
   unmodified v1.4.0 code). Compare **exactly** only the deterministic columns — `Observed K` /
   `Observed G` / `Observed g` / `Observed Interaction`, `Theoretical CSR`, `Exact CSR` — because
   the v1.4.0 permutation paths are irreproducible, so the `Permuted CSR` baseline values are not a
   stable reference. Check permuted columns loosely instead (same order of magnitude, monotone in
   `r`), and add a *new* reproducibility assertion that post-refactor permutations **are**
   bit-identical under a fixed seed *and* across `workers` values — something v1.4.0 cannot
   satisfy. Two exact-column differences are already known, expected, and must be asserted as
   differences with a comment naming the cause: `edge_correction = "none"` (whist binning vs
   `Kest`'s fast path, visible because the example data has half-integer centres) and any sample
   over `big` (which now keeps its requested correction). Spot-checked so far: translation matches
   the baseline to 7.3e-12, isotropic to 0, `none` differs by 19.1 on K ≈ 11,000 — that last one
   being the intended binning change.
3. ✅ **Memory/scale.** Verified: 200,000 cells, `edge_correction = "translation"`, **1.1 s and
   738 MB peak**; `big = 100` vs `big = 1e9` bit-identical. Worth re-running once after Phase 5–7
   to confirm nothing regressed.
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
- `apply_deprecated_args()` plus `deprecated_arg_map()` (`R/utils-deprecate.R`) is the reusable
  seam for the next rename round, and it rejects unknown `...` names so the seam cannot become a
  typo sink.
- The duplicate `get_bi_rows` is gone, and an explicit `Collate:` (Phase 7) will stop
  collation order from silently deciding which definition wins if it ever recurs.
- `Exact CSR` exists as a column for every metric even where it is `NA`, so filling it in later
  (pcf has the same random-labelling closed form as K) is additive, not breaking.
