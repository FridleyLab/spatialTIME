# Subset mif object on cellular level

This function allows to subset the mif object into compartments. For
instance a mif object includes all cells and the desired analysis is
based on only the tumor or stroma compartment then this function will
subset the spatial list to just the cells in the desired compartment

## Usage

``` r
subset_mif(mif, classifier, level, markers)
```

## Arguments

- mif:

  An MIF object

- classifier:

  Column name for spatial dataframe to subset

- level:

  Determines which level of the classifier to keep.

- markers:

  vector of

## Value

mif object where the spatial list only as the cell that are the
specified level.

## Examples

``` r
#' #Create mif object
library(dplyr)
x <- create_mif(clinical_data = example_clinical %>% 
mutate(deidentified_id = as.character(deidentified_id)),
sample_data = example_summary %>% 
mutate(deidentified_id = as.character(deidentified_id)),
spatial_list = example_spatial,
patient_id = "deidentified_id", 
sample_id = "deidentified_sample")

markers = c("CD3..Opal.570..Positive","CD8..Opal.520..Positive",
"FOXP3..Opal.620..Positive","PDL1..Opal.540..Positive",
"PD1..Opal.650..Positive","CD3..CD8.","CD3..FOXP3.")

mif_tumor = subset_mif(mif = x, classifier = 'Classifier.Label', 
level = 'Tumor', markers = markers)
```
