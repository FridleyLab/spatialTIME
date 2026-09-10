#' Marker summaries of 229 samples
#'
#' Per-sample counts and percentages for the marker phenotypes measured on each
#' sample, one row per sample.
#'
#' Beyond the two identifiers the columns come in matched pairs: a count column
#' such as `CD3 (Opal 570) Positive Cells` and its percentage counterpart
#' `% CD3 (Opal 570) Positive Cells`, for 12 phenotypes, plus `Total Cells` and two
#' analysed-area columns. Note that these names contain spaces and parentheses,
#' unlike the syntactic names used in [example_spatial].
#'
#' @format A tibble with 229 rows and 29 variables:
#' \describe{
#'   \item{deidentified_id}{patient-level id (integer)}
#'   \item{deidentified_sample}{sample-level id}
#'   \item{Total Cells}{number of cells analysed in the sample}
#'   \item{...}{12 per-phenotype positive-cell counts and their 12 matching
#'     percentage columns, then `Area Analyzed (m)` and `Area Analyzed (mm)`}
#' }
#' @seealso [example_clinical], [example_spatial]
"example_summary"
