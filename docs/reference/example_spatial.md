# Example list of 5 spatial TMA data frames

A named list of 5 spatial data frames, one per TMA core, in the form
produced by HALO cell phenotyping. Each has 51 columns: the sample
identifier, cell bounding-box coordinates (`XMin`, `XMax`, `YMin`,
`YMax`), 12 binary marker/phenotype indicator columns, per-marker
intensity and area columns, and a `Classifier.Label` column with levels
`Tumor` and `Stroma`.

## Usage

``` r
example_spatial
```

## Format

A named list of 5 data frames, each with 51 columns:

- TMA1\_\[3,B\].tif:

  3803 cells

- TMA2\_\[3,B\].tif:

  3008 cells

- TMA3\_\[7,B\].tif:

  1850 cells

- TMA3\_\[9,K\].tif:

  1803 cells

- TMA3\_\[8,U\].tif:

  2318 cells

## Details

Marker positives are sparse and vary widely between cores, which matters
when choosing a core for examples or tests: `TMA3_[9,K].tif` has 536
CD3+ and 83 CD8+ cells, whereas `TMA1_[3,B].tif` has 17 and 7. Several
markers are entirely absent from some cores.

Note that `deidentified_sample` here is **character**, while
`deidentified_id` in
[example_clinical](https://fridleylab.github.io/spatialTIME/reference/example_clinical.md)
and
[example_summary](https://fridleylab.github.io/spatialTIME/reference/example_summary.md)
is **integer**.

## See also

[example_clinical](https://fridleylab.github.io/spatialTIME/reference/example_clinical.md),
[example_summary](https://fridleylab.github.io/spatialTIME/reference/example_summary.md)
