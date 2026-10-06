# Load a disk-backed mif's spatial data back into memory

The inverse of
[`mif_to_disk()`](https://fridleylab.github.io/spatialTIME/reference/mif_to_disk.md),
named after
[`dplyr::collect()`](https://dplyr.tidyverse.org/reference/compute.html)
for the same reason: it is the point where lazily-referenced data is
materialised.

## Usage

``` r
collect_mif(mif, samples = NULL)
```

## Arguments

- mif:

  object of class `mif`.

- samples:

  optionally a subset of sample names or positions to load, so that a
  few samples can be inspected interactively without materialising a
  cohort that does not fit in memory.

## Value

object of class `mif` with an ordinary named list of data frames in its
`spatial` slot.

## See also

[`mif_to_disk()`](https://fridleylab.github.io/spatialTIME/reference/mif_to_disk.md),
[`open_mif()`](https://fridleylab.github.io/spatialTIME/reference/open_mif.md)
