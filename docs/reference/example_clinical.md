# Clinical variables of 229 patients

A tibble of clinical characteristics for 229 patients, one row per
patient.

## Usage

``` r
example_clinical
```

## Format

A tibble with 229 rows and 6 variables:

- age:

  age at diagnosis

- race:

  self-identified race

- sex:

  patient biological sex

- status:

  disease status

- deidentified_sample:

  sample identifier

- deidentified_id:

  patient identifier

## Details

`deidentified_id` is **integer** here and in
[example_summary](https://fridleylab.github.io/spatialTIME/reference/example_summary.md),
while the `deidentified_sample` column of
[example_spatial](https://fridleylab.github.io/spatialTIME/reference/example_spatial.md)
is character.

## See also

[example_spatial](https://fridleylab.github.io/spatialTIME/reference/example_spatial.md),
[example_summary](https://fridleylab.github.io/spatialTIME/reference/example_summary.md)
