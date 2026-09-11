# Marker summaries of 229 samples

Per-sample counts and percentages for the marker phenotypes measured on
each sample, one row per sample.

## Usage

``` r
example_summary
```

## Format

A tibble with 229 rows and 29 variables:

- deidentified_id:

  patient-level id (integer)

- deidentified_sample:

  sample-level id

- Total Cells:

  number of cells analysed in the sample

- ...:

  12 per-phenotype positive-cell counts and their 12 matching percentage
  columns, then `Area Analyzed (m)` and `Area Analyzed (mm)`

## Details

Beyond the two identifiers the columns come in matched pairs: a count
column such as `CD3 (Opal 570) Positive Cells` and its percentage
counterpart `% CD3 (Opal 570) Positive Cells`, for 12 phenotypes, plus
`Total Cells` and two analysed-area columns. Note that these names
contain spaces and parentheses, unlike the syntactic names used in
[example_spatial](https://fridleylab.github.io/spatialTIME/reference/example_spatial.md).

## See also

[example_clinical](https://fridleylab.github.io/spatialTIME/reference/example_clinical.md),
[example_spatial](https://fridleylab.github.io/spatialTIME/reference/example_spatial.md)
