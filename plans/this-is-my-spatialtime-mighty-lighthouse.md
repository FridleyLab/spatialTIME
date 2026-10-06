# Disk-backed `mif` objects

## Context

**Nothing like this exists yet.** Grepping `R/` for `arrow|parquet|feather|duckdb|fst|HDF5|DelayedArray|
lazy|on_disk|mmap|readRDS|saveRDS` finds no storage layer, and no plan file proposes one —
`plans/create-an-implementation-plan-squishy-turing.md:343` mentions `arrow` only to *rule it out* as a
dependency for a `split_tissue()` figure. This is new work.

`mif$spatial` is a named `list` of in-memory `data.frame`s (`R/create_mif.R:84-91`). All seven metrics
then do the identical thing:

```r
out = parallel::mclapply(seq_along(mif$spatial), function(sample_i){
  spat = mif$spatial[[sample_i]]          # closure captures the WHOLE mif
  ...
}, mc.cores = workers, mc.preschedule = FALSE)
```

The parent must hold every sample resident for the whole run. Measured on the author's own cohort in
`SplittingTissue/whole_slide/` (283 ovarian WSIs, **342,267,952 cells**, 21 columns):

| | |
|---|---|
| parquet on disk | **2.19 GB** |
| the same data as R `data.frame`s (100.0 bytes/row, measured on 3 samples, spread < 0.1%) | **34.2 GB** |
| largest sample (3,247,618 cells), all 21 columns | 325 MB |
| largest sample, the 4 columns a metric actually reads | 91 MB |
| this machine | 38.7 GB RAM, 12 cores |

Every metric reads only `xloc`/`yloc` (or `XMin/XMax/YMin/YMax`), the `sample_id` column, and the 1–2
marker columns of the moment. 17 of 21 columns are never touched on a given pass.

This is the author's own working pattern, formalised: `SplittingTissue/1.0.splitting.R:5-15` already builds
a named character vector of parquet paths and calls `arrow::read_parquet(f)` *inside*
`mclapply(mc.cores = 16, mc.preschedule = FALSE)`, and `SplittingTissue/functions.R:107-144` already
implements `IndexedList` — a directory plus `index.rds`, with `[` reading one element on demand.

### What this does and does not fix — measured, because two analyses disagreed by 10x

Per worker, Ripley's K retains the pair list from `k_pairs()` (`R/utils-k-engine.R:125-161`: `i`, `j` int
+ `d`, `w` double = 24 bytes/pair) alongside the frame. Measured on a real 1,000,977-cell slide
(`example_whole-slides/Peres_P3_110361 A7`, λ = 9.33e-4):

| `max(r_range)` | close pairs | pair list | frame (21 col) |
|---|---|---|---|
| 10 | 20,368 | 0.0005 GB | 0.100 GB |
| 30 | 4,323,730 | 0.10 GB | 0.100 GB |
| 50 | 12,716,078 | 0.31 GB | 0.100 GB |
| **100** (the default) | 51,248,348 | **1.23 GB** | 0.100 GB |

So **the pair list overtakes the frame at `max(r_range) ≈ 30`**, and at the default `r_range = seq(0, 100, 1)`
it is ~12x the frame. Disk-backing does not touch it. Anyone claiming this change alone makes whole-slide K
cheap is wrong.

The benefit is a different term. Two memory pools add:

- **O(cohort), parent-resident, for the entire run** — 34.2 GB. Independent of `r_range`, `workers`, or
  what you compute. Disk-backing takes this to a ~115 KB index (measured on all 283 samples: the paths
  and sample names dominate it, and one shared schema is stored rather than 283 copies).
- **O(n·λ·r²), per worker, transient** — the pair list. Scales with sample size and radius, *not* with
  cohort size. Disk-backing leaves it alone; bounding it is what `big`/`chunk_index()` already exist for.

At `workers = 12`, `r_range = seq(0, 100, 1)` on this cohort:

```
today       34.2 GB parent  +  12 x (0.10 frame + 1.23 pairs)  =  50.2 GB   -> exceeds 38.7 GB, cannot start
disk-backed  ~0    parent  +  12 x (0.03 frame + 1.23 pairs)  =  15.1 GB   -> fits
```

This was then measured end-to-end rather than left as arithmetic — on the 5 files in
`example_whole-slides/` (3,353,797 cells, 335.4 MB as a spatial slot, 100.0 bytes/row), comparing the
current code against a paths-only `mclapply` closure doing `read_parquet(col_select = c("x","y","CD3+"))`
in each child, `workers = 2`:

| | `r_range = seq(0, 30, 1)` | `seq(0, 100, 1)` (the default) |
|---|---|---|
| parent peak RSS | 974 -> **347 MB (-64%)** | 952 -> **338 MB (-64%)** |
| per-worker peak | 1278 -> 1294 MB (**rises**) | 6476 -> 7206 MB (**rises**) |
| tree-total peak | -648 MB (**-22.3%**) | -38 MB (-0.4%, noise) |

So disk-backing is **necessary but not sufficient**. The win is *entirely and only* the parent's resident
copy — reliably −64% here, and that is the term that scales with cohort size: 335 MB at 5 samples, 34.2 GB
at 283. The per-worker peak does not improve and in fact rises slightly (arrow's own read buffers), so at
the default `r_range` the tree-total barely moves on a 5-sample cohort. That is expected and is not the
point: it removes a fixed admission price you cannot trade against `workers`. The honest headline is
*memory stops scaling with cohort size*, which is the whole reason a 283-sample run becomes possible at
all. Bounding the per-worker pair list is separate work, and `big`/`chunk_index()` already exist for it.

### Four constraints that shape the design

1. **`arrow` is a hard dependency — `Imports`.** Author's decision: one reliable code path, no
   `requireNamespace()` guards, no silently-degraded half-support. Two consequences that must be handled
   *before* any code lands, because they are currently blocking:
   - `arrow` is **not installed** in the development environment (`/opt/anaconda3/envs/spatialTIME`,
     R 4.6.1 — nor are `fst`/`duckdb`), and `SplittingTissue/1.0.splitting.R:12` records that it
     *"doesn't seem to want to install"* on the author's HPC. With `arrow` in `Imports`,
     `library(spatialTIME)` fails outright in both places, not just the disk path. Installing it in both
     is **Step 0** and gates everything else.
   - CRAN `arrow` requires `R (>= 4.2)` + C++20, so `DESCRIPTION:73-74` must rise
     `Depends: R (>= 4.1)` → **`R (>= 4.2)`**. That is user-visible and belongs in `NEWS.md`. The
     reverse-suggests `mxfda` and `scSpatialSIM` inherit `arrow` transitively.
2. **`read_parquet()` only — never `open_dataset()`.** See §2; this is a silent-corruption issue, verified.
3. **The six-slot list is a public contract.** spatialTIME **1.3.4-5 is on CRAN** (2024-06-04), over nine
   archived tarballs back to 2021-05-14, and `man/create_mif.Rd:32-42` documents all six slots by name in
   `\value`. The CRAN-shipped vignettes *teach* direct slot indexing — `mif$spatial[[1]]` at
   `vignettes/deriving_functions.Rmd:31`, `vignettes/intro.Rmd:361-362`,
   `vignettes/spatialexperiment.Rmd:162` — and the reverse-suggests `mxfda` reads
   `$derived$univariate_Count` by name. No slot may be added, renamed or reordered
   (`tests/testthat/test-create-mif.R:8-9` pins all six), and `mif$spatial[[i]]` must keep working.
4. **`split_tissue()` writes to the spatial slot.** `R/split_tissue.R:300,305,309` assign
   `density_compartment`, `refined_density_compartment` and `density_score` into `mif$spatial[[i]]`.
   Read-only is not sufficient.

---

## Design

### 1. A `mif_store` in the `spatial` slot — six slots preserved

Rather than add a slot, a disk-backed mif puts an S3 `mif_store` object in `mif$spatial`: a small
data frame of sample name, relative path, row count and column schema, with the store root and format
version in attributes. It gets `length()`, `names()`, `[[`, `c()` and `print()` methods.

The payoff: `length(mif$spatial)` (`R/print.R:21`) and `mif$spatial[[i]]` keep working unchanged, so
**existing third-party code — including `mxfda` and `scSpatialSIM` — works against a disk-backed mif**,
with `[[` transparently reading that sample from disk. Compatibility is maximal by construction rather
than by audit.

Internally the metrics use a *projecting* accessor instead, in a new `R/utils-mif-store.R`:

```r
mif_spatial(mif, i, columns = NULL)   # one sample, only these columns
mif_samples(mif); n_mif_samples(mif)  # metadata, loads nothing
mif_spatial_set(mif, i, new_columns)  # write path, §6
spatial_columns(mif, mnames = NULL, xloc = NULL, yloc = NULL, extra = NULL)
```

`spatial_columns()` derives the projection from what each metric already knows — `mif$sample_id`, plus
`c(xloc, yloc)` or `c("XMin","XMax","YMin","YMax")`, plus `mnames`, plus `extra` (`classifier` for
`split_tissue()`/`subset_mif()`, `cell_type` for `plot_immunoflo()`) — so the list is computed once, not
hand-maintained in nine places.

For an **in-memory** mif, `mif_spatial()` returns `mif$spatial[[i]]` and **ignores `columns` entirely**.
That is deliberate: the in-memory path stays byte-identical, so no released result can move. The cost is
that a projection missing a column a metric actually uses fails *only* on the disk path — which is what
the parity suite in §Verification exists to catch.

### 2. The reader rule: `read_parquet(col_select=)`, never `open_dataset()`

`arrow::open_dataset() |> select() |> collect()` **does not preserve row order.** Verified on
`example_whole-slides/Peres_P3_110361 A7` (1,000,977 rows, 1 row group, arrow 23.0.1.2, R 4.5.0):
**935,427 of 1,000,977 positions differ** from the full read, first divergence at index **32,769** — the
2^15 record-batch boundary — with an identical sorted multiset.

That is disqualifying here, because this codebase is positional throughout: `R/utils-k-engine.R:185`
(`sel <- keep_i[pairs$i] & keep_j[pairs$j]`) indexes masks into a pair list built in file order, and
`R/split_tissue.R:298-310` assigns three columns back by position. A reordered marker column would attach
markers to the wrong cells and return plausible, wrong numbers — the worst failure mode available.

`read_parquet()` is order-stable. Verified across **all 283 files** and the full row-group distribution
present in the cohort (116 files with 1 row group, 140 with 2, 25 with 3, 2 with 4; row groups are
2^20 rows): order-stable **283/283** — equal to row-group concatenation order, across repeated reads,
across differing column subsets, and in **191/191** `mclapply(mc.preschedule = FALSE)` children checked
against digests taken in the parent before the fork. Over the same 283 files
`open_dataset() |> collect()` matched file order in only **72/283**. `Object Id` is unique per file, so it
is a valid join key should a sidecar ever need one.

**Rule for the implementation:** the parquet backend calls `arrow::read_parquet(file, col_select = ...)`
and nothing else. No `open_dataset`, no lazy dplyr, no predicate pushdown. A comment at the call site
must say why, with the 32,769 number, or someone will "modernise" it later.

### 3. Fork safety — measured

| test | result |
|---|---|
| `read_parquet` in child, arrow cold in parent | 8/8 ok, 0.19 s |
| `read_parquet` in child, arrow warmed in parent | 8/8 ok, 0.10 s |
| parent-created `open_dataset()` pointer used in child | 4/4 ok |
| 24 files / 8 cores, arrow threads on | 24/24 ok, **0.26 s, 24.8 M cells** |
| same, `set_cpu_count(1)` in child | 24/24 ok, 0.31 s |

I/O is free relative to the statistic: the largest sample reads in **69 ms** full, **60 ms** projected.
Each backend still re-opens in the child rather than sharing the parent's pointer — that happening to
work is an arrow internal, not a contract.

**Do not pin arrow's thread count inside the child.** This was tried — `set_cpu_count(1)` per worker, as
insurance for a large-core HPC — and it **deadlocks on arrow 25.0.1 / R 4.6.1**: the child inherits a
copy of the thread pool's mutex state from the fork and then waits forever on a lock nothing in the
child holds. The run hangs at 0% CPU with no error. It also measured *slower* than leaving arrow alone
(0.31 s against 0.26 s), so it was pure cost. Reading in a child is safe on both arrow 23.0.1.2 and
25.0.1; only the resize hangs. To limit threads, set the count in the parent and let children inherit.

This is exactly the failure the "re-verify in the conda env" caveat below was written to catch: arrow
23.0.1.2 tolerated the call and 25.0.1 does not.

### 4. Three entry points

```r
create_mif(..., spatial_list = <named character vector of paths>)   # reference mode
mif_to_disk(mif, path, overwrite = FALSE)                           # `path` REQUIRED
open_mif(path);  collect_mif(mif, samples = NULL)
```

Reference mode is a **pure addition**: `spatial_list` is `stopifnot()`-ed to be a list of data frames
(`R/create_mif.R:44-45`), so a character vector errors today and no existing code changes. It points at
parquet you already have with zero copying — exactly `1.0.splitting.R:5-6`. Each file's footer is read
once for schema and row count (~1 ms, no data read).

**`mif_to_disk()`'s `path` has no default and is required** — you name exactly where the store lands, and
nothing is ever written to a directory you did not choose:

```r
mif_to_disk(m, "~/projects/ova/cohort.mif")
mif_to_disk(m, "/scratch/$USER/cohort.mif")     # HPC scratch
mif_to_disk(m)                                  # Error: `path` must name the directory
                                                #   to write the store to.
```

This follows the precedent set by `sigma` in `split_tissue()`, which was deliberately left without a
default (`NEWS.md:19-22`) so that a silent choice could never make two runs incomparable. Writing tens of
GB somewhere the user did not name is the same class of mistake. `path` is expanded (`path.expand()`) and
normalised, and `overwrite = FALSE` errors if a store already exists there rather than merging into it.

`collect_mif()` is the inverse, named after `dplyr::collect()`. `print.mif` gains a backend line and a
manifest-derived cell count, so printing a 342 M-cell mif is instant.

**Paths must be relative to the manifest and no live handle may enter the slot.** `saveRDS(mif)` has to
keep working — `R/utils-density-boundary.R`'s `split_tissue_settings()` already backfills settings for
"a mif saved before they existed", so cross-session persistence is an existing promise. That rules out
arrow `Dataset`/`externalptr` objects in the slot by construction.

### 5. Storage layout — the verifiable folder

```
cohort.mif/
  manifest.json          <- the verification surface
  clinical.parquet  sample.parquet
  spatial/  <sample>.parquet     one per sample, never rewritten
  overlay/  <sample>.parquet     columns added later (§6); absent until needed
  derived/  <slot>.parquet       data-frame derived slots
            <slot>.rds           list-valued slots (density_boundary, spatial_plots)
```

`manifest.json` holds format version, backend, `patient_id`/`sample_id`, and per sample: file name, row
count, column names and types, and a cheap integrity stamp (size + mtime + first/last-row digest, not a
checksum of 2.19 GB). `open_mif()` validates all of it *before* any compute, so a moved or truncated store
fails with `sample X: manifest says 711,676 rows, file has 32` rather than a `spatstat` error 40 minutes in.
This is `IndexedList`'s `index.rds` (`functions.R:112`) plus types, counts and integrity.

List-valued derived slots stay RDS, reusing the type dispatch at `R/merge_mifs.R:108-131` that already
exists because `bind_rows()` silently destroys them.

### 6. The write path

`split_tissue()` appends three columns per sample. An **overlay** file holds only those columns;
`mif_spatial()` reads base + overlay and `cbind`s, with overlay shadowing base on a name clash (which is
what `split_tissue(overwrite = TRUE)` needs — `R/split_tissue.R:229-232` already errors on a clash
otherwise). Rewriting each spatial file would rewrite 2.19 GB to add ~12 MB and would **mutate the
author's proprietary source data** in reference mode. Source stays immutable; the store is append-only.

`subset_mif()` (`R/subset_mif.R:49`) filters rows and builds a new mif via `create_mif()`. On a
disk-backed mif it takes a `path` and writes a new store; omitting `path` errors with the fix rather than
silently materialising the 34 GB the user adopted disk-backing to avoid.

### 7. Backend — parquet only

Because `arrow` is now `Imports` (constraint 1), **parquet is the only spatial backend.** There is no
`rds`/`fst` fallback and no backend negotiation: one format, always present, always order-stable. This is
a direct simplification bought by the dependency decision — it halves the parity matrix in §Verification
and removes the "which backend am I on?" branch from every accessor.

Two things survive that are easy to confuse with a fallback:

- **RDS is still used for list-valued `derived` slots** (`density_boundary`, `spatial_plots`) — see §5.
  That is required regardless of the spatial backend, because those slots are not tabular.
- The read/write/probe interface stays a **3-function shim** rather than inlined `arrow` calls. Not for a
  second backend today, but so the order-stability rule in §2 lives in exactly one place instead of nine.
  If `arrow` ever does prove unavailable somewhere, an `rds` backend behind that shim is ~20 lines — it
  would lose column projection (the 100/28 = 3.6x width factor) but keep the O(cohort) win, since
  per-sample granularity is what matters and projection is a bonus. Document that, do not build it.

---

## Files

**Create** — `R/utils-mif-store.R` (accessors, manifest, backend dispatch, `mif_store` methods),
`R/mif_to_disk.R` (`mif_to_disk`, `open_mif`, `collect_mif`), `tests/testthat/test-mif-store.R`,
`tests/testthat/test-disk-backed-parity.R`.

**Modify**
- `R/create_mif.R` — accept a named character vector; every existing `stopifnot()` for the
  list-of-data-frames case unchanged.
- The nine read sites — `ripleys_k.R:108`, `bi_ripleys_k.R:110`, `nn_g.R:84`, `bi_nn_g.R:94`,
  `pair_correlation.R:73`, `bi_pair_correlation.R:74`, `interaction.R:91`, `split_tissue.R:251`,
  `plot_tissue_split.R:129` — become `mif_spatial(mif, i, spatial_columns(...))`, and
  `length(mif$spatial)` becomes `n_mif_samples(mif)`. Mechanical and identical in each;
  `R/ripleys_k.R:104-113` is the reference (and already the only metric that projects columns).
  `dixons_s.R:91` and `plot_immunoflo.R:62` additionally move from iterating the list to iterating indices.
- `R/split_tissue.R:300,305,309` — `mif_spatial_set()`.
- `R/subset_mif.R` — `path` argument. `R/print.R` — backend + cell count.
- `R/merge_mifs.R` — its spatial handling is `do.call(c, .)` over the lists (`:74-77`), which needs a
  `c.mif_store` method for two disk-backed mifs (concatenate manifests, re-check duplicate names at
  `:78-80`); mixing a disk-backed and an in-memory mif errors and names `collect_mif()` as the fix.
- `DESCRIPTION` — add `arrow` and `jsonlite` to **`Imports`**; bump `Depends: R (>= 4.1)` →
  **`R (>= 4.2)`** (forced by `arrow`). No `Suggests` change; `fst`/`duckdb` are not used.
- `NEWS.md` — the new capability, and the `R (>= 4.2)` floor called out as user-visible.
- `_pkgdown.yml`, `man/` + `NAMESPACE` (regenerate), `R/global.R` only if check flags it.

---

## Verification

**1. Parity is the contract.** A helper in `tests/testthat/helper-mif.R` that, for each of the seven
metrics, runs it on `example_mif()` and on `mif_to_disk(example_mif(), withr::local_tempdir())` and
asserts `expect_equal()` on the derived table — at `workers = 1` and `workers = 2` (CRAN's limit,
`NEWS.md:145`). Extend `tests/testthat/test-regression-vs-baseline.R` so the disk path is also pinned
against `fixtures/baseline-v1.4.0.rds`. Nothing else matters if this fails.

**2. Row order, explicitly.** A test that reads a fixture's coordinate columns and marker columns in
separate calls and asserts `identical()` against the full read. This is the §2 corruption mode; assert it
rather than trusting it, and it will also catch an arrow upgrade that changes the guarantee.

**3. Round-trip.** `collect_mif(mif_to_disk(m, tmp))$spatial` `expect_identical` to `m$spatial` including
column order, types and **factor levels** — `split_tissue()` stores factors (`R/split_tissue.R:300`) and
parquet factor round-tripping is exactly where this breaks.

**4. Store validation fails loudly.** Delete a spatial file; truncate one; edit a row count in
`manifest.json`; open a directory with no manifest; bump the format version; move the store after writing
it (relative paths, §4). Each must error naming the sample and the discrepancy. `withr::local_tempdir()`,
matching existing `withr` usage.

**5. Compatibility.** Assert `names(mif)` is still the six slots for both modes (mirroring
`test-create-mif.R:8-9`), and that `length(mif$spatial)`, `names(mif$spatial)` and `mif$spatial[[i]]` all
work on a disk-backed mif — the `mxfda`/`scSpatialSIM` surface. Also `saveRDS()`/`readRDS()` a
disk-backed mif and re-run a metric.

**6. Memory, measured — not extrapolated.** The two analyses behind §Context's table disagreed by 10x
until measured, so measure this too. Generate a synthetic cohort into a temp dir (never the proprietary
data), then record four numbers per run — parent resident, per-worker peak, total peak, and the share
attributable to the frame versus the pair list — via `/usr/bin/time -l` plus a background `ps` sampler.
Do **not** measure a forked child's RSS with `ps` alone: macOS reports shared pages separately and made an
earlier attempt at this uninterpretable. Assert `object.size(mif$spatial)` is O(KB) and flat in cohort size.

**7. Against the real cohort, locally, never committed.** `SplittingTissue/` is already in both
`.gitignore` and `.Rbuildignore`.
```r
paths <- list.files("SplittingTissue/whole_slide", pattern = "parquet$", full.names = TRUE)
names(paths) <- sub("\\.gz\\.parquet$", "", basename(paths))
m <- create_mif(clinical, sample_tbl, spatial_list = paths, sample_id = "Image Tag")
m <- ripleys_k(m, mnames = "CD3+", r_range = seq(0, 30, 1), workers = 8, permute = FALSE)
```
Success: it completes with peak RSS near `workers x (0.03 + pair list)` rather than approaching 34.2 GB.
Use `r_range = seq(0, 30, 1)` first — at `seq(0, 100, 1)` the pair list is 1.23 GB/worker and *that*, not
the frame, sets the ceiling. Today neither runs on a 38.7 GB machine.

**8. Step 0 — install `arrow` in both environments; nothing else can start until this passes.**
`arrow` in `Imports` means the package does not load without it, and it is currently **missing from the
development environment**. So, in order:

```bash
conda install -n spatialTIME -c conda-forge r-arrow     # R 4.6.1 env, currently has no arrow
/opt/anaconda3/envs/spatialTIME/bin/Rscript -e \
  'library(arrow); packageVersion("arrow"); nrow(read_parquet(<one WSI file>))'
```
Then the same on the HPC, which is the known-hard case (`1.0.splitting.R:12`). If it will not install
there, **stop and revisit constraint 1** — that is the decision point where a `Suggests` + `rds` fallback
comes back on the table, and it is much cheaper to find out now than after nine call sites are converted.
Note `conda search` could not be reached from this sandbox (SSL cert failure against
`conda.anaconda.org`), so availability of `r-arrow` for `osx-arm64` is assumed, not verified.

**Environment.** Development and verification happen in the **`spatialTIME` conda env:
`/opt/anaconda3/envs/spatialTIME`, R 4.6.1**, library `/opt/anaconda3/envs/spatialTIME/lib/R/library`
(84 packages; `jsonlite` 2.0.0, `spatstat.*` 3.8.3, `testthat` 3.3.2, `withr` 3.0.3 present —
`arrow`, `fst`, `duckdb`, `SpatialExperiment` absent). This is also the R 4.6.1 that
`plans/review-this-repository-and-fizzy-donut.md:91` was reaching for; that plan's
`/opt/homebrew/lib/R/4.6/site-library` path is an orphan (homebrew R is gone) — the real one is conda.

**Caveat on every number in this plan:** all measurements above were taken in **R 4.5.0 at
`/usr/local/bin/R`** with arrow 23.0.1.2, because that is the only R here that currently has arrow.
Re-run the row-order assertions (§2) and the memory measurements (Verification 6) in the conda env once
arrow is installed there, and treat any disagreement as a finding rather than noise — the row-order
behaviour in particular is an arrow-version property, which is exactly why Verification 2 asserts it in
the suite rather than trusting this plan's numbers.

**Fixture and benchmark home: `dev/`.** It is `.Rbuildignore`d (`^dev$`) so it never ships to CRAN, but
despite the `dev/` line in `.gitignore` it *is* tracked — `git ls-files dev/` returns
`example_data_creation.Rmd` and `.md`, added before the ignore rule — so files there are versioned.
Nothing else qualifies: `R_dev/` was deleted in `ec510b0`, there is no `data-raw/`, and
`data/example_spatial.rda` has **never had a generating script** (committed as a binary in `48a4fa5`, no
`@source` tag anywhere; `dev/example_data_creation.Rmd` is dead documentation for a different, deleted
dataset). Tests themselves generate parquet fixtures into `tempdir()` and commit nothing.
