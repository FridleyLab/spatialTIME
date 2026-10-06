# Reopen a mif store written by `mif_to_disk()`

Reads and validates `manifest.json`, re-probes every spatial file, and
returns the `mif`. Validation happens up front so that a moved,
truncated or half-written store fails immediately with a message naming
the sample, rather than as a confusing error partway through a long
computation.

## Usage

``` r
open_mif(path)
```

## Arguments

- path:

  directory containing the store.

## Value

object of class `mif`.

## See also

[`mif_to_disk()`](https://fridleylab.github.io/spatialTIME/reference/mif_to_disk.md),
[`collect_mif()`](https://fridleylab.github.io/spatialTIME/reference/collect_mif.md)
