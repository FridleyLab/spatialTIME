# Marker Frequency Difference

marker_freq_diff allows users to easily calculate the difference in
marker abundance (frequency) between different tissue compartments. For
us this is typically looking at tumor vs stroma marker abundance.

## Usage

``` r
marker_freq_diff(
  mif,
  classifier,
  ref_level,
  diff_level,
  mnames,
  overwrite = FALSE
)
```

## Arguments

- mif:

  object of class `mif` created with
  [`create_mif()`](https://fridleylab.github.io/spatialTIME/reference/create_mif.md)

- classifier:

  character that specified which column in the spatial data indicates
  the different tissue compartments profiled

- ref_level:

  character found in the `classifier` column of spatial data indicating
  the first value when calculating frequency difference

- diff_level:

  character found in the `classifier` column of spatial data indicating
  the second value when calculating frequency difference

- mnames:

  character vector of marker names that are found in the spatial data
  for which to calculate abundance differences

- overwrite:

  boolean/logical for whether to overwrite previously calculated marker
  frequency differences

## Value

an object of class `mif` with the `derived$frequency_difference` slot
filled

## Details

To calculate a p-value associated with the difference in proportion of
positive markers identified in different classifier compartments, we
implemented the use of the Fisher's Exact test.

## Examples

``` r
mif <- create_mif(clinical_data = example_clinical %>% 
                                 dplyr::mutate(deidentified_id = as.character(deidentified_id)),
                               sample_data = example_summary %>% 
                                 dplyr::mutate(deidentified_id = as.character(deidentified_id)),
                               spatial_list = example_spatial,
                               patient_id = "deidentified_id", 
                               sample_id = "deidentified_sample")
mif = marker_freq_diff(mif = mif,
                       classifier = "Classifier.Label",
                       ref_level = "Tumor", 
                       diff_level = "Stroma",
                       mnames = c("CD3..FOXP3.", "CD3..CD8.", "CD3..PD1.", "CD3..PD.L1.",
                                  "CD3..Opal.570..Positive", "PD1..Opal.650..Positive"),
                       overwrite = TRUE)
```
