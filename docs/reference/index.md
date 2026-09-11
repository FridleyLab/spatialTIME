# Package index

## Nearest Neighbor G

Univariate and Bivariate Nearest Neighbor G(r) functions

- [`NN_G()`](https://fridleylab.github.io/spatialTIME/reference/NN_G.md)
  : Univariate Nearest Neighbor G(r)
- [`bi_NN_G()`](https://fridleylab.github.io/spatialTIME/reference/bi_NN_G.md)
  : Bivariate Nearest Neighbor G(r)

## Ripley’s K

Univariate and Bivariate Ripley’s K functions

- [`ripleys_k()`](https://fridleylab.github.io/spatialTIME/reference/ripleys_k.md)
  : Calculate Ripley's K
- [`bi_ripleys_k()`](https://fridleylab.github.io/spatialTIME/reference/bi_ripleys_k.md)
  : Bivariate Ripley's K
- [`bi_ripleys_k_WSI()`](https://fridleylab.github.io/spatialTIME/reference/bi_ripleys_k_WSI.md)
  : Bivariate Ripley's K for Whole Slide Images (removed)

## Pair Correlation

Univariate and Bivariate Pair Correlation functions

- [`pair_correlation()`](https://fridleylab.github.io/spatialTIME/reference/pair_correlation.md)
  : Univariate Pair Correlation Function
- [`bi_pair_correlation()`](https://fridleylab.github.io/spatialTIME/reference/bi_pair_correlation.md)
  : Bivariate Pair Correlation Function

## MIF Object Management

Functions for creating, merging, and subsetting Multiplex
Immunoflourescent objects

- [`create_mif()`](https://fridleylab.github.io/spatialTIME/reference/create_mif.md)
  : Create Multiplex Immunoflourescent object
- [`spatial_exp_to_mif()`](https://fridleylab.github.io/spatialTIME/reference/spatial_exp_to_mif.md)
  : Create Multiplex Immunoflourescent object from a SpatialExperiment
  object
- [`merge_mifs()`](https://fridleylab.github.io/spatialTIME/reference/merge_mifs.md)
  : Merge several MIF objects together
- [`subset_mif()`](https://fridleylab.github.io/spatialTIME/reference/subset_mif.md)
  : Subset mif object on cellular level

## Other Spatial Statistics

- [`dixons_s()`](https://fridleylab.github.io/spatialTIME/reference/dixons_s.md)
  : Dixon's S Segregation Statistic
- [`interaction_variable()`](https://fridleylab.github.io/spatialTIME/reference/interaction_variable.md)
  : Bivariate Interaction Variable
- [`marker_freq_diff()`](https://fridleylab.github.io/spatialTIME/reference/marker_freq_diff.md)
  : Marker Frequency Difference

## Tissue Architecture

Density-based segmentation of a sample into tissue compartments

- [`split_tissue()`](https://fridleylab.github.io/spatialTIME/reference/split_tissue.md)
  : Split a sample into tissue compartments by density difference

## Visualization

- [`plot_immunoflo()`](https://fridleylab.github.io/spatialTIME/reference/plot_immunoflo.md)
  : Generate plot of TMA point process
- [`plot_tissue_split()`](https://fridleylab.github.io/spatialTIME/reference/plot_tissue_split.md)
  : Plot density-based tissue compartments alongside the density
  difference

## Data

Example datasets included with the package

- [`example_clinical`](https://fridleylab.github.io/spatialTIME/reference/example_clinical.md)
  : Clinical variables of 229 patients
- [`example_spatial`](https://fridleylab.github.io/spatialTIME/reference/example_spatial.md)
  : Example list of 5 spatial TMA data frames
- [`example_summary`](https://fridleylab.github.io/spatialTIME/reference/example_summary.md)
  : Marker summaries of 229 samples

## Internal

- [`` `%>%` ``](https://fridleylab.github.io/spatialTIME/reference/pipe.md)
  : Pipe operator
