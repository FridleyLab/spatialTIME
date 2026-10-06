# Establishing output trust for spatialTIME 2.0.0 (`feature-upgrade`)

## Context

`feature-upgrade` is a large refactor of `feature-alex` (which is **v1.4.0** and a strict
*ancestor* — `git merge-base` is feature-alex's own tip `20960eb`). The refactor replaced the
Ripley's K implementation with one shared pair-list engine, delegated nearest-neighbour G to
spatstat, unified every metric's output schema, rewrote the permutation/RNG scheme, and added
two new capabilities (density-based tissue segmentation, and — still uncommitted — disk-backed
`mif` objects).

The question this plan answers is not "does it run" but **"if someone runs `feature-upgrade`,
can they trust the degree of clustering and the exact results?"** Permutations are random, but
the underlying process must produce the same *expected* results, and the disk-backed path must
be byte-for-byte reproducible.

Investigating that question surfaced **four numerical bugs in the statistics** and **ten
correctness gaps in the disk-backed store**, all verified empirically (evidence inline below).
Two of the statistics bugs silently produce wrong numbers with no warning, so this plan fixes
code first and only then pins it with tests — testing a broken estimator would just certify
the breakage.

There is already a narrow v1.4.0 baseline (`tests/testthat/fixtures/baseline-v1.4.0.rds`:
1 sample, 3 markers, 7 radii) consumed by `tests/testthat/test-regression-vs-baseline.R`.
That file's per-delta pattern is right; it is too narrow to establish trust. This plan widens
it and makes it contract-driven.

## Decisions taken

| Decision | Choice |
|---|---|
| `edge_correction = "border"` | Implement the real reduced-sample estimator **inside the engine**, keeping pair-list reuse |
| `derived` lost on disk-backed mifs | **Fix** — persist it, then test both directions |
| Parity deliverable | **Wide committed fixture + committed, auditable generator**, excluded from CRAN checks |
| `Permutations Larger than Observed` | Switch the count to `>=` (ties raise the p-value) **and** add a proper p-value column, to all 7 metrics, named in the `Degree of Clustering *` style |

## Keeping it out of CRAN's way

`.Rbuildignore` is exactly the right mechanism — a path listed there is excluded from
`R CMD build`, so `R CMD check` on the tarball never sees it. Two notes:

- `dev/` is already in **both** `.Rbuildignore` and `.gitignore`, so the generator cannot live
  there (it would not be committed). Use **`tools/parity/`** and add `^tools$` to `.Rbuildignore`.
- The wide fixture goes in `tests/testthat/fixtures/` with an `.Rbuildignore` entry, and its
  test opens with `skip_if_not(file.exists(path))` — **already the house idiom**
  (`test-regression-vs-baseline.R:22`). Net effect: runs under `devtools::test()` from a git
  checkout, skips under `R CMD check` on the tarball. A development-only suite, zero CRAN weight.

---

## Phase 0 — Statistics bugs (must land before Phase 1)

### BUG-1 · `k_from_pairs()` early return drops the NA truncation
`R/utils-k-engine.R:186` — `if (!any(sel)) return(rep(0, length(r)))` bypasses the
`k[r >= pairs$rmax_valid] <- NA_real_` on line 189.

Verified: window with `boundingradius = 67.1`, `r = c(0,20,40,60,80)`, isotropic. Selecting the
two most distant points (denominator `2*1 > 0`, but no pair within `rmax`) returns
`0 0 0 0 0`; the normal path returns `0 1368.08 5246.19 11313.72 NA`.

Consequences: the NA contract is inconsistent *within one column*; under `permute = TRUE`,
`rowMeans(permuted, na.rm = TRUE)` mixes those zeros with NAs so `Permuted CSR` beyond
`rmax_valid` collapses toward 0 and biases `Degree of Clustering Permutation`; **and it breaks
the unbiasedness identity Phase 1 relies on.**

Fix: build the vector on one path, apply the truncation once before returning.

### BUG-2 · `Permutations Larger than Observed` is 0 where the observation is NA
`rowSums(permuted > observed, na.rm = TRUE)` returns `0` when every term is NA — which reads
as "no permutation exceeded the observation", i.e. maximal clustering. Call sites:
`R/ripleys_k.R:161`, `R/bi_ripleys_k.R:164`, `R/nn_g.R:114`, `R/bi_nn_g.R:132`,
`R/pair_correlation.R:103`, `R/bi_pair_correlation.R:122`, `R/interaction.R:136`.

Fix as part of the p-value rework below.

### BUG-3 · `edge_correction = "border"` silently returns the uncorrected estimator
`k_pairs()` has branches for `"translation"` (`R/utils-k-engine.R:134`) and `"isotropic"`
(`:147`) only. `"border"` falls through with `w = rep_len(1, npair)`, so it is **bit-identical
to `"none"`** — while `match_edge_correction()` accepts it (`:230`) and the roxygen documents
it as supported (`:107`, `R/ripleys_k.R:8`).

Verified on 400 uniform points: `max|border − none| = 0`, `max|border − Kest$border| = 174`.
feature-alex passed it through to `Kest(correction = 'border')` genuinely, so this is a
**silent regression**, not merely a missing feature. `test-k-engine.R:13`'s `CORR` vector
excludes border, which is why nothing caught it.

Fix — reduced-sample estimator in the engine. Border is not a per-pair weight, so it needs a
second denominator path:
- In `k_pairs()`, compute `bd <- spatstat.geom::bdist.points(pp)` once per sample and carry it.
- In `k_from_pairs()`, for each `r` restrict to eligible points (`bd > r`), count pairs from
  eligible `i`, and divide by the eligible count rather than `n_i(n_i - 1)`.
- `rmax_valid` stays as spatstat sets it for border.

This keeps the pair list and therefore the permutation reuse. `border` becomes the only
correction whose denominator depends on `r` — new structure in `k_from_pairs()`, so it needs
its own validation (below). Mirror into `bi_ripleys_k()`'s bivariate denominator.

**DONE, with one consequence the plan did not anticipate: `Exact CSR` must be `NA` for
border.** The unbiasedness identity in Phase 1 needs a denominator fixed once `m` is fixed.
Translation, isotropic and none all have one (`m(m-1)/area`); border's counts the *selected*
cells still beyond `r` from the edge, so it varies across relabellings and
`E[numer/denom] != E[numer]/E[denom]`. Measured on 600 cells, `m = 120`, `B = 3000`: the K of
all cells sits 0.3% below the mean permuted K at `r = 20` rising to 1.0% at `r = 80`
(z = −0.9 to −7.9) — real bias, growing with `r`, not Monte Carlo noise. A column named
*Exact* CSR cannot be 1% wrong, so border returns `NA` there, matching the existing precedent
in `NN_G`/`pair_correlation`. Implemented in `R/ripleys_k.R` and `R/bi_ripleys_k.R`, documented
in a new "The CSR columns, and which to trust" roxygen section, and pinned by the final test in
`test-csr-null.R` (which also asserts the bias, so a future change that "helpfully" fills the
column in fails with an explanation).

Verified: border now agrees with `Kest(correction = c("border","translation"))$border` and the
matching `Kcross()` at `max|diff| = 0` — continuous and integer coordinates, whole samples and
marker subsets, including the `NaN` region past the inradius.

### BUG-4 · `Exact CSR` is hard-NA whenever `permute = TRUE`
`R/ripleys_k.R:131` and `R/bi_ripleys_k.R:127`:
`exact = if(permute) rep(NA_real_, length(r_range)) else k_from_pairs(pairs, rep(TRUE, n))`.

So a user cannot see both CSR references in one run — precisely the comparison that answers
"are my permutations behaving". Fill it unconditionally: it is one extra `k_from_pairs()` call
per sample against an already-built pair list (negligible). This also puts the Phase 1
convergence diagnostic in the user's own output.

### The p-value rework (all 7 metrics)

Context that shaped this: the strict-`>` raw count is **inherited verbatim from feature-alex**,
where it existed in only 4 of 7 metrics — `bi_ripleys_k` (as
`Permutations Larger than Observed`) and `pair_correlation` / `bi_pair_correlation` /
`interaction_variable` (as `Permuted_larger_than_Observed`). It was **absent from
`ripleys_k`, `NN_G` and `bi_NN_G`**. The refactor did not change the formula; the shared
`standard_metric_cols()` promoted it to all 7.

Changes:
1. Count switches to `>=`, so permutations tying the observed value count as "not more
   clustered" and raise the p-value. Conservative, and the standard choice.
2. New column `Permutation p-value` = `(1 + #{perm >= obs}) / (B + 1)`, space-separated
   title case to match the `Degree of Clustering *` family. Added to
   `standard_metric_cols()` in `R/utils-schema.R` so all 7 metrics get it uniformly.
   *(Name is a one-line change if you'd prefer another.)*
3. Both columns are `NA` where `observed` is `NA` (BUG-2), and `NA` when `permute = FALSE`.
4. Divide by the count of non-NA permutations, not blindly by `B`.

Touches: `R/utils-schema.R` (`standard_metric_cols`), the seven `*_result_frame()` helpers,
`R/global.R` (globalVariables), `tests/testthat/test-output-schema.R` (table-driven — add the
column once), and a NEWS entry.

---

## Phase 1 — Seed-free validation of the CSR null

This is the direct answer to the user's question. Since
`Degree of Clustering X = Observed − X CSR`, correctness decomposes into **(i)** `Observed` is
exactly spatstat — already pinned by `test-k-engine.R` — and **(ii)** the CSR reference is an
unbiased estimate of the null expectation. Phase 1 establishes (ii).

**The key fact, proven and verified: `E[Permuted CSR] = Exact CSR` is an *exact* identity, not
asymptotic.** `sample.int(n, m)` is SRSWOR so `P(i and j both marked) = m(m−1)/(n(n−1))`, and
the denominator is `m(m−1)/area`, giving
`E[K_perm(r)] = area·Σ w·1{d≤r}/(n(n−1)) = k_from_pairs(pairs, rep(TRUE, n))`. The bivariate
scheme gives `P(i∈anchor ∧ j∈counted) = n_i·n_j/(n(n−1))` against a denominator of
`n_i·n_j/area` — same result, **which also proves `bi_ripleys_k`'s `Exact CSR` is right.**
Holds for every `r`, every `m ≥ 2`, every edge correction.

Verified two ways: by complete enumeration (max relative error 5.4e-16 univariate,
3.2e-15 bivariate) and by Monte Carlo at B=4000 (z ≈ 1.0–1.7 across radii, both schemes).

New file **`tests/testthat/test-csr-null.R`**, in value-per-runtime order. Total < 2 s.

1. **Enumeration unbiasedness** (~0.04 s). `n = 10`, all `C(10,3) = 120` labellings; mean of
   `k_from_pairs()` over every labelling must equal `k_from_pairs(pairs, rep(TRUE, n))` to
   `1e-12`, for `translation`/`isotropic`/`none`/**`border`**. Bivariate twin: `n = 8`,
   `n_i = 2`, `n_j = 3`, all 560 configurations. Exact proof, no RNG, no guessed tolerance.
   **Include a radius beyond `rmax_valid`** — that makes this the cheapest BUG-1 detector, and
   it is also the primary validation of the new border denominator.
2. **Same-estimator checks** (~ms). A marker positive for *every* cell must give
   `Observed K == Exact CSR` (both are `k_from_pairs(pairs, rep(TRUE, n))`); if the observed
   path ever forks from the CSR path these stop being equal. This kills the v1.4.0
   separate-hand-rolled-estimator class of bug.
3. **Non-degenerate null** (~ms). With `keep_permutation_distribution = TRUE`,
   `length(unique(Permuted CSR)) == num_permutations` at a non-degenerate radius. Only
   detector of "all permutations identical", which has a *correct mean* and would pass 1, 4
   and 5. Must assert the precondition: at small `r` permuted values tie heavily (measured: 20
   of 300 distinct at r=5, exactly 1 at r=0), so pick `r >= 20` and fail loudly if the fixture
   drifts.
4. **Two-sided power** (< 1 s). Clustered marker → `Degree of Clustering Exact > 0` and
   `Permutation p-value` small; inhibited marker → `< 0` and p large. Measured at r ≥ 20 the
   separation is unambiguous (clustered `+5111…+48689`; inhibited `−1200…−3106` with
   `larger/B ≥ 0.987`). Assert inequalities at `r >= 20`, never `== 0`/`== B`.
   **Needs a new `inhibited_marker_mif()` in `tests/testthat/helper-mif.R`** (greedy hard-core
   thinning, ~10 lines) — *no under-dispersion fixture exists anywhere in the suite today.*
5. **Self-studentised convergence** (~0.03 s). `z(r) = sqrt(B)·(mean_B − Exact)/sd_B`, which is
   `N(0,1)` under the CLT, so the tolerance is *derived* rather than guessed. Assert
   `max|z| < 5` at B=200 (per-radius flake 5.7e-7; radii are strongly positively correlated so
   family-wise flake stays < 1e-6). Two mandatory guards: drop radii where `sd_B == 0` (r=0 is
   always `NaN`) and where `Exact CSR` is `NA`.
6. **p-value calibration** (~2 s, `skip_on_cran()`). Build the marker *by* random labelling so
   the null is true by construction, then the p-values must be ~Uniform. Measured design:
   n=500 continuous coords, m=80, r₀=50, R=200 replicates, B=100 — `mean(p) = 0.497`,
   5-bin counts `48 37 31 42 42`, χ² = 4.05. Assert `chisq_5bin < qchisq(0.9999, 4)` (= 23.51)
   and `|mean(p) − 0.5| < 0.085` (4σ). Use the binned χ², not a rejection rate at α=0.05
   (almost no power, one tail only) and not KS (`p` is discrete on `{0, 1/B, …, 1}`). Add a
   small R=20 replicate driven through `ripleys_k()` itself so the wrapper's bookkeeping is
   covered, not just the engine's. Same tie precondition as test 3.

Also add to **`tests/testthat/test-k-engine.R`**:
- `border` into the `CORR` vector (`:13`) so it is checked against `Kest$border` like the others.
- A `big` invariance test for `ripleys_k`: `expect_equal(run(big = 500), run(big = 1e9))`.
  `?ripleys_k:49-51` claims `big` is memory-only, and only `bi_ripleys_k`/`_WSI` have such a test.

**Rejected:** the `1/sqrt(B)` rate-ratio assertion. `|mean_B − exact|` is half-normal, so the
ratio of single-realisation errors at B=100 vs B=400 is `2·|Cauchy(0,1)|` and
`P(error does not shrink at all) = 29.5%`. Rescuable by averaging over ~100 replicates, but
that measures the same `σ/sqrt(B)` as test 5 in 5 s instead of 0.03 s. Phase 4 only.

---

## Phase 2 — Disk-backed store: fixes then tests

All ten items below were verified empirically. Each is a **code fix first**, then a test.

| Item | Verified behaviour | Fix |
|---|---|---|
| Derived slot lost on reopen | `ripleys_k(disk_mif)` populates `derived`; `open_mif(root)` returns it **empty**. `sync_manifest()` has exactly one call site, `R/split_tissue.R:384` | Call `sync_manifest()` at the end of `write_derived()` when `is_disk_mif(mif)`, with a `sync = TRUE` opt-out |
| Overlay validation | Replaced a 2500-row overlay with a 3-row wrong-schema file; **`open_mif()` accepted it**, failing later with an opaque error | Extract `open_mif`'s per-sample checks into `validate_store_file()` and call it for overlays too, validating against the `overlay_nrow`/`overlay_columns` already in format 1 — **no `MIF_STORE_FORMAT` bump needed** |
| List columns | `mif_to_disk()` succeeds; `collect_mif()` returns an `AsIs` column **not identical** to the input. Silent corruption | Reject pre-write in `mif_to_disk()` naming the sample and columns |
| `raw` / `complex` | `raw` **silently becomes `integer`**; `complex` fails with arrow's own `"Cannot infer type from vector"` | Same pre-write rejection, so the user gets our message |
| `collect_mif(in_memory, samples = "nope")` | Returns a 1-element list named `NA`; the disk path errors correctly | Route the in-memory branch through a list twin of `store_subscript()`; assert **both paths give the same message** in one `test_that` so they cannot drift |
| `c.mif_store()` format version | Merged stores with `format` 1 and 99 without complaint; result claims `format = 1` | Error on mixed formats |
| Reference-mode staleness | Replaced the source parquet with a smaller one — `attr(store, "nrow")` still stale, nothing checks, and `mif_spatial_set()`'s row guard uses the stale number | Record `bytes`/`nrow` in `create_mif()`'s character branch (it already probes) and add an exported `verify_mif(mif)` re-probe verb, usable in both modes |
| Partial store | `manifest.json` missing → clean error. **`sample.rds` missing → raw `gzfile` warning + bare `readRDS` error** (same for `clinical.rds`, `derived/*.rds`) | `read_store_rds()` that checks existence and reports `"Store is incomplete"` |
| Missing `mif_store` methods | `[[<-` and `$` error with poor messages, but `d$spatial[1] <- list(df)` **silently yields a plain `list`** and `length(d$spatial) <- 0` **silently yields a plain `character`** | Add `$.mif_store` (must *work* — vignettes and mxfda use `mif$spatial$Name`) plus `[[<-`/`[<-`/`length<-` that error with a shared read-only message pointing at `collect_mif()` |
| `open_mif()` id validation | A manifest with a wrong `sample_id` yields a mif where every metric fails deep inside `add_cell_centres()` | Assert `m$sample_id %in% columns` per sample in the existing loop |

Then extend the tests. **`tests/testthat/test-disk-backed-parity.R`** gets the derived-reopen
round-trip (including two metrics in sequence, asserting manifest slot order is preserved).
**`tests/testthat/test-mif-store.R`** gets the rest, highest-value first:

- **Heterogeneous cohort** — the single most valuable missing store test. `collapse_schema()`
  makes `attr(spatial, "columns")` a *list* when samples differ, and
  `store_columns`/`store_types`/`[.mif_store`/`c.mif_store` all branch on `is.list()`;
  `test-mif-store.R:494` only asserts the **homogeneous** case. Build a two-sample store where
  sample 2 has an extra column and exercise `[`, `c`, `collect_mif`, and a metric.
- Projection **column order**: `mif_spatial(mif, i, columns)` vs `mif_spatial(mif, i, NULL)`.
  `plot_tissue_split()` depends on it; no test covers it.
- Overlay **shadowing** a base column — promised at `R/utils-mif-store.R:298-304` and relied on
  by `split_tissue(overwrite = TRUE)`, never tested.
- Column types, table-driven: `logical`/`integer`/`Date`/`POSIXct` (incl. `tzone`) round-trip
  exactly — measured; document that `integer64` returns as `integer` when every value fits in
  32 bits; `raw`/`complex` refused with our message.
- `mif_spatial_set()` **directly** — row-count guard, reference-mode refusal, and the
  "re-run replaces its own columns" merge logic, independent of the 21 KB `split_tissue.R`.
- Colliding sample names round-tripped through `open_mif()` on the *original* store
  (`:416` only checks `subset_mif`'s store).
- Zero-row and one-row samples: `pq_probe` on empty parquet, the all-overlay
  `data.frame(row.names = seq_len(0))` branch, `print.mif_store`'s row-count format.
- `c.mif_store()` with duplicate sample names (reachable directly as an exported S3 method).

---

## Phase 3 — Cross-branch parity harness and wide fixture

Run **after Phase 0 and 2**, so the pinned "new" side is the version you actually want.

### Mechanics

Two processes, because both trees declare `Package: spatialTIME` and one R session cannot hold
both namespaces regardless of `lib.loc`. **Do not install either side** — a CRAN spatialTIME
1.5.0 is already in the system library and a third copy makes "which did I measure?" a real
hazard. All of feature-alex's Imports are installed (purrr 1.0.4, furrr 0.3.1, future 1.70.0,
tidyselect 1.2.1) and neither tree has `src/`, so `load_all()` compiles nothing.

New, all `.Rbuildignore`d via `^tools$`:
- `tools/parity/capture.R` — load one side, run the call table, write one RDS
- `tools/parity/compare.R` — diff two RDS, print a per-`(entry, column)` verdict table, write the fixture
- `tools/parity/run.sh`:

```bash
WT=$(mktemp -d /tmp/st-alex.XXXX); OUT=$(mktemp -d /tmp/st-parity.XXXX)
git worktree add --detach "$WT" 20960eb
Rscript tools/parity/capture.R --pkg="$WT"  --side=old --out="$OUT/old.rds"
Rscript tools/parity/capture.R --pkg="$PWD" --side=new --out="$OUT/new.rds"
Rscript tools/parity/compare.R "$OUT/old.rds" "$OUT/new.rds" \
  --fixture=tests/testthat/fixtures/parity-v1.4.0.rds
git worktree remove "$WT"
```

`--pkg="$PWD"` is the **live working tree** — the disk-backed feature is untracked, and
`git stash create` has no `--include-untracked`, so git plumbing cannot snapshot the new side.
If an immutable snapshot is ever wanted, `rsync -a --exclude='.git'`.

`capture.R` requirements:
- `pkgload::load_all(pkg, export_all = TRUE, helpers = FALSE, attach_testthat = FALSE)` —
  `export_all` to reach `k_pairs`/`k_from_pairs`; `helpers = FALSE` because feature-alex has
  its own `helper-mif.R` that would shadow.
- **Version interlock, mandatory:** assert `packageVersion("spatialTIME")` is `1.4.0` / `2.0.0`
  and that `pkgload::pkg_path()` is the intended directory. **Never write `spatialTIME::`
  anywhere** — bare names only, or you silently measure the installed 1.5.0.
- One call table, one per-side argument translator (`keep_permutation_distribution` →
  `keep_perm_dis`, inject `method = "K"`). `tryCatch` each call and store the condition message
  — several old entries are *expected* to error, and that is part of the record.
- `workers = 1` everywhere (feature-alex's nested `mclapply` forks uncontrollably).

### Fixture design — `tests/testthat/fixtures/parity-v1.4.0.rds`

Store **full tables**, not checksums: ~60 entries at 30–60 rows each lands well under 200 KB
(the existing narrow fixture is 11.8 KB), and the full table is the only thing that can pin
column names, column order, `iter` values, output class and row counts — which are most of the
delta list. Exceptions stored as structural summaries: ggplots (`nrow(data)`, labels, mappings,
layer classes — a serialised ggplot carries its whole environment) and `density_boundary`
(`call_info` + classes, not the `psp` objects).

- **Samples (2, unthinned):** `TMA3_[9,K].tif` (1803 cells; 536/83/34 positives) as the primary,
  and `TMA1_[3,B].tif` (3803 cells; 17/7/3) as the degenerate case — it is where the row-count
  deltas actually live (`pair_correlation` `return(NULL)`, `interaction_variable` drops,
  `bi_NN_G` sparse stubs).
- **Markers:** CD3/CD8/FOXP3 positives univariate. Bivariate pairs CD8/FOXP3 (the only workable
  real pair) **and** CD3/CD8 (near-total nesting — pins the one-counted-cell degeneracy).
- **`r_range = c(0,10,20,30,40,50,60,650,700,1300)`.** The first seven preserve comparability
  with the existing fixture; `650`/`700` straddle `TMA3`'s isotropic `rmax_valid = 689.8` and
  `1300` crosses its translation `rmax_valid = 1276.0`. **Those three radii are the only thing
  that pins the new NA-truncation behaviour.**
- **Corrections:** `translation`, `isotropic`, `none`, **`border`**. Capturing old `border` at
  `permute = FALSE` (where feature-alex genuinely called `Kest(correction='border')`) gives the
  Phase 0 border fix an independently-verifiable target — do not skip it.
- **One `big = 500` entry** on `TMA3`: feature-alex's `if(nrow(spat) > big) edge_correction = 'none'`
  (`R/ripleys_k.R:97`) makes the old side downgrade *and* take the tiled `getTile()` branch, so
  one capture exercises both legacy behaviours. No 10k-cell synthetic fixture needed.
- **`permute = FALSE`** for every numeric entry, plus one `permute = TRUE, num_permutations = 5,
  keep_permutation_distribution = TRUE` entry per metric recorded as `NOT_COMPARABLE`
  numerically — kept to pin `iter`, column names/order, class and row counts.
- **Functions:** `ripleys_k`, `bi_ripleys_k`, `NN_G` (rs/han/none/km), `bi_NN_G` (rs/han/km),
  `pair_correlation`, `bi_pair_correlation`, `interaction_variable`, `dixons_s`,
  `marker_freq_diff`, `subset_mif`, `plot_immunoflo`.
- **Attributes:** per-side `spatialTIME_version`, `git_sha`, `git_dirty`, `R.version`, platform,
  **`spatstat.{explore,geom,univar}` versions** (agreement with spatstat *is* the claim, so a
  spatstat bump must be attributable), `dixon`, `arrow`, `md5sum` of every `R/*.R`, timestamp;
  plus the `design` list and the old→new `rename_map`.

### The contract drives the test

`attr(fx, "contract")` is a data frame of
`entry, column, verdict ∈ {MUST_MATCH, MUST_DIFFER, NOT_COMPARABLE}, tolerance, reason,
independent_check` — and **`compare.R` exits non-zero if any measured result contradicts its
declared verdict**, so the generator is itself a check.

`tests/testthat/test-parity-v1.4.0.R` is then a loop over that contract with three
`expect_*` branches, plus a guard test asserting every captured `(entry, column)` appears in
the contract. For each `MUST_DIFFER` row, `independent_check` pins the **new** value to
spatstat rather than merely asserting difference — generalising what
`test-regression-vs-baseline.R:79-87` does by hand. This is what keeps the delta list
maintainable instead of 40 hand-written blocks that drift.

Contract rows the deltas require (each with its reason recorded):

| Verdict | Entries |
|---|---|
| `MUST_MATCH` | `Observed K` (translation/isotropic, `permute=FALSE`), `Exact CSR`, `Theoretical CSR` (= old `Theoretical G`/`Theoretical g`), `Observed G` (rs/han/none), `Observed g`, `Observed Interaction`, dixon `Obs.Count`/`Exp.Count`/`S`/`Z`/`p-val.Z`/`P.asymp`, `marker_freq_diff` counts and `%`, `subset_mif` counts, **`Observed K` border** (old `Kest$border` is the Phase 0 target) |
| `MUST_DIFFER` | `none` binning (`whist` left-closed vs `Kest` fast-path `d<=r`; new == `Kest(correction=c("none","translation"))$un`), `big=500` (old downgraded + tiled), `NN_G(km)` structure (old malformed: 5-column positional reorder leaked `theohaz.x`, dropped `Marker`), `marker_freq_diff` p-values (old double-counted the margin *and* substring-matched `CD3..CD8.` into `CD3..CD8..FOXP3.`), `subset_mif` `%` (now ×100), `plot_immunoflo` `cell_type`, `Observed K` at `r >= rmax_valid` (new NA-truncates), `Permutations Larger than Observed` for the 4 metrics that had it (now `>=`), row counts for sparse pcf/interaction/`bi_NN_G`, `iter` values, output class |
| `NOT_COMPARABLE` | every `Permuted CSR` / `Degree of Clustering Permutation` (v1.4.0's 24 nested `mclapply` calls omit `mc.cores`, so `set.seed()` had no effect), `Observed K` on the `permute=TRUE` path (feature-alex used a separate hand-rolled `d < r` estimator there), dixon `P.rand`/`p-val.Nobs` |
| *new column* | `Permutation p-value`, and `Permutations Larger than Observed` for `ripleys_k`/`NN_G`/`bi_NN_G` (absent on feature-alex) |

Retire the old fixture **last**: get `test-parity-v1.4.0.R` green, then delete
`baseline-v1.4.0.rds`, `test-regression-vs-baseline.R` and
`test-disk-backed-parity.R:216-245` in one commit.

### What Phase 3 actually found

All 30 entries agree with the contract, but four results differed from what this plan
predicted. Recorded because each narrows a claim in the delta table above.

**A new correctness bug in 1.4.0's `interaction_variable()`, not previously known.**
`get_bi_rows()` returns rows marker-major, so `cells$cell` is not globally ascending
— but `subset(sample_ppp, cells$cell)` returns points in *sorted* order, and
`marks(ps) <- cells$Marker` then attached labels in marker-major order to points in
index order. 1.4.0 measured distances between **scrambled marker sets**. Verified by
brute force on `TMA3_[9,K].tif` with FOXP3/CD8 at `r = 20`: exactly 4 of 109 anchors
lie within 20 units (fifth-nearest is 20.55), so `3.669725` is right; 1.4.0 reported
`4.587156` (= 5/109), and re-running its exact code path reproduces that. Both sides
divide by 109, so the denominator was never the issue. `Observed Interaction` moves
from MUST_MATCH to MUST_DIFFER.

**The `"none"` binning delta is far narrower than documented.** It requires *evenly
spaced* `r`, because that is the condition for `Kest` to take its fast C path — with
the truncation radii in `r_range` the two versions agree exactly. It is also confined
to quantities `Kest` computes: in `bi_ripleys_k` on an even grid `Observed K` agrees
(from `Kcross`, which has no fast path) while `Exact CSR` differs (from `Kest` over
all cells, which does), in the same call. Needed two extra capture entries on an even
grid to exercise at all.

**`border`'s `Observed K` is unchanged from 1.4.0**, which this plan did not predict.
1.4.0 passed `border` through to `Kest` on its `permute = FALSE` path, so the new
reduced-sample engine reproduces those values exactly. The regression was that the
correction was silently replaced by `"none"` on the `permute = TRUE` path and above
`big` — not that the estimator itself was wrong there.

**`Exact CSR` for an unestimable marker.** 1.4.0 reported the K of all cells for a
marker with fewer than 3 positives (that quantity does not depend on the marker);
2.0.0 NAs the whole stub row. Defensible either way, left as-is, recorded as
MUST_DIFFER so the change is not silent.

Two deltas this plan predicted needed deliberate fixtures to exercise at all, because
nothing in the shipped markers triggers them: the sparse-marker row drops needed
`CD3..PD.L1.` (1 positive) for `pair_correlation`, and the non-nested `CD8`/`PD1`
pair for `interaction_variable` — `FOXP3` at 3 sits just inside the `< 3` guard, and
every other sparse marker is entirely nested inside CD3, which empties the counted
set on *both* cores rather than one.

---

## Phase 4 — Slow diagnostics (manual, never in `R CMD check`)

`tools/parity/slow-null-diagnostics.R`: the `1/sqrt(B)` rate check at R=100 replicates per B,
p-value calibration at R=2000, and a full 5-core × all-corrections × `permute=TRUE` sweep.

---

## Verification

Existing suite runtimes are small (`test-disk-backed-parity.R` 4.0 s / 51 pass,
`test-mif-store.R` 2.4 s / 112 pass, `test-k-engine.R` 0.9 s / 38 pass), so the budget is open.

1. **Phase 0 lands:** `devtools::test(filter = "k-engine|csr-null|output-schema")`.
   Border must agree with `Kest(correction = "border")` to `1e-12` on both the continuous and
   integer-grid fixtures. BUG-1 is confirmed fixed by the enumeration test including a radius
   past `rmax_valid`. Re-run the border check I used:
   `max|engine_border − Kest$border|` must go from 174 to ~0.
2. **Phase 1 lands:** `devtools::test(filter = "csr-null")` green and under ~3 s. Then run it
   20× with different `set.seed()` values to confirm no flake before committing.
3. **Phase 2 lands:** `devtools::test(filter = "mif-store|disk-backed")`. Specifically confirm
   `ripleys_k(disk_mif)` → `open_mif(root)$derived$univariate_Count` is now identical, and that
   the four store-corruption cases I verified as silently accepted now error with messages
   naming the sample.
4. **Phase 3 lands:** `bash tools/parity/run.sh` exits 0 and its verdict table is reviewed **by
   hand** against the delta table above. Then `devtools::test(filter = "parity")`.
5. **Full suite, both ways:**
   - `devtools::test()` from the checkout — the wide-fixture tests **run**.
   - `devtools::check(document = FALSE)` — the wide-fixture tests **skip** (fixture is
     `.Rbuildignore`d), and check is clean. Confirm with
     `R CMD build . && tar -tzf spatialTIME_2.0.0.tar.gz | grep -E 'tools/|parity-v1'`
     → **expect no output**.
6. **Both R installs.** Everything above was measured under R 4.5.0 at `/usr/local/bin/R`
   (arrow 23.0.1.2). `plans/create-an-implementation-plan-squishy-turing.md:218` and
   `plans/this-is-my-spatialtime-mighty-lighthouse.md:379-384` both require re-running under
   `/opt/anaconda3/envs/spatialTIME/bin/Rscript` (R 4.6.1) — the arrow thread-pool deadlock was
   only ever observed there.
7. **NEWS.md** gains entries for: the border fix (and that it was silently `"none"`), the
   `>=` tie convention, the new `Permutation p-value` column, `Exact CSR` now populated under
   `permute = TRUE`, the NA-truncation and `larger`-with-NA fixes, and derived-slot persistence
   on disk-backed mifs.

## Open items

- **Column name** `Permutation p-value` — chosen to match the `Degree of Clustering *`
  space-separated style; trivial to rename before Phase 0 commits.
- **`sync = TRUE` opt-out** on `write_derived()`: a metric run now mutates the store. Cheap
  (one manifest + RDS write per call) but it is I/O the caller did not request.
- **`interaction_variable()` normalisation** remains unconfirmed against Steinhart et al. —
  flagged at `NEWS.md:371-376` and carried forward unchanged. Out of scope here, but it is the
  one remaining place where a *documented* oddity could be a genuine error.
- **CI is stale and arrow-naive.** `.github/workflows/tic.yml` is a 2020-11-14 template
  (`actions/checkout@v2.3.4`, `setup-r@master`, `::set-output`) and pins no arrow system
  dependency, though `arrow` is now a hard `Imports` needing C++20. None of this plan's tests
  will run in CI until that is modernised.
