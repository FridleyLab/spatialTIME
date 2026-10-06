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

  Either a named list of data frames with the spatial data from each
  sample making up each individual data frame, or a **named character
  vector of paths to parquet files**, one sample per file, which builds
  a disk-backed mif that reads each sample only when a computation needs
  it. The files are referenced where they are and never copied or
  modified. Use this for whole slide images, where holding every sample
  in memory at once is what stops a multi-core run from starting – see
  [`mif_to_disk()`](https://fridleylab.github.io/spatialTIME/reference/mif_to_disk.md)
  for the sizes involved and for building a store from data already in
  memory. If the vector is unnamed, names are taken from each file's
  `sample_id` column.

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
x <- create_mif(clinical_data = example_clinical %>% 
mutate(deidentified_id = as.character(deidentified_id)),
sample_data = example_summary %>% 
mutate(deidentified_id = as.character(deidentified_id)),
spatial_list = example_spatial,
patient_id = "deidentified_id", 
sample_id = "deidentified_sample")
```
