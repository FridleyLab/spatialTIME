# Merge several MIF objects together

This function merges MIF objects that were run separately so they can be
used as a single MIF. MIF objects don't *need* but *should* have the
same column names in the summary file and clinical data file. The MIF
objects **DO** need to have the same patient_id and sample_id.

## Usage

``` r
merge_mifs(mifs = NULL, check.names = T)
```

## Arguments

- mifs:

  A list of MIF objects to merge together

- check.names:

  whether to check names of spatial files and summary enttries

## Value

Returns a new MIF object list

- clinical_data:

  clinical information from all

- sample:

  cell level summary data from all

- spatial:

  contains all spatial files from all MIFs

- derived:

  appended derived variables

- patient_id:

  patient_id from the first MIF - this is why it is important to have
  the same patient_id for all MIFs

- sample_id:

  sample_id from the first MIF - also important for all MIFs to have the
  same sample_id

## Examples

``` r
#merge several MIF objects
library(dplyr)
x <- create_mif(clinical_data = example_clinical %>% 
mutate(deidentified_id = as.character(deidentified_id)),
sample_data = example_summary %>% 
mutate(deidentified_id = as.character(deidentified_id)),
spatial_list = example_spatial,
patient_id = "deidentified_id", 
sample_id = "deidentified_sample")
x <- merge_mifs(mifs = list(x, x), check.names = FALSE)
#> No variables have been derived yet
```
