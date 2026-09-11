# Create Multiplex Immunoflourescent object

Creates an MIF object for use in spatialIF functions

## Usage

``` r
create_mif(
  clinical_data,
  sample_data,
  spatial_list = NULL,
  patient_id = "patient_id",
  sample_id = "image_tag"
)
```

## Arguments

- clinical_data:

  A data frame containing patient level data with one row per
  participant.

- sample_data:

  A data frame containing sample level data with one row per sample.
  Should at a minimum contain a 2 columns: one for sample names and one
  for the corresponding patient name.

- spatial_list:

  A named list of data frames with the spatial data from each sample
  making up each individual data frame

- patient_id:

  A character string indicating the column name for patient id in sample
  and clinical data frames.

- sample_id:

  A character string indicating the column name for sample id in the
  sample data frame

## Value

Returns a custom MIF

- clinical:

  Data frame of clinical data

- sample:

  Data frame of sample data

- spatial:

  Named list of spatial data

- derived:

  List of data derived using the MIF object

- patient_id:

  The column name for sample id in the sample data frame with the
  clinical data

- sample_id:

  The column name for sample id in the sample data frame to merge with
  the spatial data

## Examples

``` r
#Create mif object
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
x <- create_mif(clinical_data = example_clinical %>% 
mutate(deidentified_id = as.character(deidentified_id)),
sample_data = example_summary %>% 
mutate(deidentified_id = as.character(deidentified_id)),
spatial_list = example_spatial,
patient_id = "deidentified_id", 
sample_id = "deidentified_sample")
```
