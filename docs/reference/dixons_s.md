# Dixon's S Segregation Statistic

This function processes the spatial files in the mif object, requiring a
column that distinguishes between different groups i.e. tumor and stroma

## Usage

``` r
dixons_s(
  mif,
  mnames,
  num_permutations = 1000,
  type = c("Z", "C"),
  workers = 1,
  overwrite = FALSE,
  xloc = NULL,
  yloc = NULL
)
```

## Arguments

- mif:

  An MIF object

- mnames:

  vector of markers corresponding to spatial columns to check Dixon's S
  between

- num_permutations:

  Numeric value indicating the number of permutations used. Default is
  1000.

- type:

  a character string for the type that is wanted in the output which can
  be "Z" for z-statistic results or "C" for Chi-squared statistic
  results

- workers:

  Integer value for the number of workers to spawn

- overwrite:

  Logical value determining if you want the results to replace the
  current output (TRUE) or be to be appended (FALSE).

- xloc:

  a string corresponding to the x coordinates. If null the average of
  XMin and XMax will be used

- yloc:

  a string corresponding to the y coordinates. If null the average of
  YMin and YMax will be used

## Value

An object of class `mif`. Depending on `type`, one or both of these
tables is written to the `derived` slot.

`Dixon_Z` (`type = "Z"`), one row per ordered marker pair per sample:

- \<sample_id\>:

  sample identifier

- From, To:

  the two cell types being compared, as
  [`dixon::dixon()`](https://rdrr.io/pkg/dixon/man/dixon.html) labels
  them

- Anchor, Counted:

  which requested marker pair the row belongs to. The per-marker count
  columns are pair-specific, because dual positives are dropped per
  pair, so without this the counts cannot be interpreted

- Obs.Count, Exp.Count, S, Z, p-val.Z, p-val.Nobs:

  Dixon's statistics

- columns:

  cells of each marker retained for that pair

- Simulations, Run:

  permutations requested, and which call produced the row

`Dixon_C` (`type = "C"`), three rows per marker pair per sample:

- \<sample_id\>:

  sample identifier

- Direction:

  "Overall segregation", or "From "

- Anchor, Counted:

  as above

- df, Chi-sq, P.asymp, P.rand:

  chi-squared segregation test

- Simulations, Run:

  as above

Marker pairs with fewer than 3 cells of either type are returned with
`NA` statistics rather than dropped, so the table shape does not depend
on cell counts.

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
```
