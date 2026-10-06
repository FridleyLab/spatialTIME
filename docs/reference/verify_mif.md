# Re-probe every spatial file a mif points at

Checks that each file still has the row count and columns the mif
expects, and reports the first that does not.

[`open_mif()`](https://fridleylab.github.io/spatialTIME/reference/open_mif.md)
does this automatically, but two cases are not covered by it. A mif
built with `create_mif(spatial_list = <parquet paths>)` references your
files directly and has no manifest, so there is nothing to validate at
construction and nothing stops those files changing afterwards. And a
long-running session holds a store open across hours of computation,
during which the files can be replaced underneath it. In both cases the
recorded row count goes stale, and because the metrics index cells by
position a stale count means markers attached to the wrong cells rather
than an error.

## Usage

``` r
verify_mif(mif)
```

## Arguments

- mif:

  object of class `mif`. In-memory mifs have nothing to verify and pass
  trivially, so this is safe to call unconditionally.

## Value

invisibly `TRUE`; errors naming the first sample that fails.

## See also

[`mif_to_disk()`](https://fridleylab.github.io/spatialTIME/reference/mif_to_disk.md),
[`open_mif()`](https://fridleylab.github.io/spatialTIME/reference/open_mif.md)

## Examples

``` r
x <- create_mif(clinical_data = spatialTIME::example_clinical,
  sample_data = spatialTIME::example_summary,
  spatial_list = spatialTIME::example_spatial[1],
  patient_id = "deidentified_id", sample_id = "deidentified_sample")
store <- file.path(tempdir(), "verify.mif")
xd <- mif_to_disk(x, store, overwrite = TRUE)
verify_mif(xd)
unlink(store, recursive = TRUE)
```
