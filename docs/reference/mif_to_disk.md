# Write a mif's spatial data to disk

Materialises an in-memory `mif` into a self-describing directory and
returns a disk-backed `mif` that reads from it. Metric functions then
load only the columns they need, for the one sample a worker is
handling, instead of the parent process holding every sample resident
for the whole run.

## Usage

``` r
mif_to_disk(mif, path, overwrite = FALSE)
```

## Arguments

- mif:

  object of class `mif` created with
  [`create_mif()`](https://fridleylab.github.io/spatialTIME/reference/create_mif.md).

- path:

  directory to write the store to. **Required, with no default** – a
  store can be tens of gigabytes, so it is always written where you say
  and never somewhere chosen for you. `~` is expanded.

- overwrite:

  whether to replace an existing store at `path`.

## Value

object of class `mif` whose `spatial` slot reads from `path`.

## Details

Use this when the spatial data no longer comfortably fits in memory,
which in practice means whole slide images. On a 283-slide cohort of 342
million cells the in-memory spatial slot is 34.2 GB against 2.19 GB of
parquet, and since the parent process must hold all of it while
`mclapply()` forks, a multi-core run cannot start at all. Disk-backing
makes peak memory scale with `workers` rather than with the number of
samples.

It is not a general speed-up. For
[`ripleys_k()`](https://fridleylab.github.io/spatialTIME/reference/ripleys_k.md)
the per-worker pair list overtakes the spatial frame at roughly
`max(r_range) = 30`, so at the default `r_range` the frame is a minor
part of a worker's footprint; what disk-backing removes is the parent's
copy of the whole cohort.

The store is a plain directory you can inspect, copy and archive:


    cohort.mif/
      manifest.json      format, ids, and per-sample rows/columns/size
      clinical.rds  sample.rds
      spatial/  <sample>.parquet    one file per sample, never rewritten
      overlay/  <sample>.parquet    columns added later by split_tissue()
      derived/  <slot>.rds

Paths inside `manifest.json` are relative, so the whole directory can be
moved or copied and still open.

## What is and is not preserved

Column names, column order, types and factor levels – including unused
levels and `NA`s – round-trip exactly, so
`collect_mif(mif_to_disk(x, p))` returns the spatial data `x` had.

**Row names are not preserved**; they come back as `1:nrow`. A columnar
file has nowhere to put them, and storing them would mean an extra
column per sample for something no function in this package reads –
every metric indexes cells by position. Non-default row names only arise
from having subset a data frame (`spat[keep, ]`), and `dplyr` resets
them anyway, so this is visible only if you compare with
[`identical()`](https://rdrr.io/r/base/identical.html) after such a
subset. Use `all.equal(check.attributes = FALSE)`, or reset them
yourself, if you need that comparison to hold.

## See also

[`open_mif()`](https://fridleylab.github.io/spatialTIME/reference/open_mif.md)
to reopen one,
[`collect_mif()`](https://fridleylab.github.io/spatialTIME/reference/collect_mif.md)
to load it back into memory.

## Examples

``` r
x <- create_mif(clinical_data = spatialTIME::example_clinical,
  sample_data = spatialTIME::example_summary,
  spatial_list = spatialTIME::example_spatial[1],
  patient_id = "deidentified_id", sample_id = "deidentified_sample")
store <- file.path(tempdir(), "example.mif")
xd <- mif_to_disk(x, store, overwrite = TRUE)
xd
#> 229 patients spanning 229 samples and 1 spatial data frames were found 
unlink(store, recursive = TRUE)
```
