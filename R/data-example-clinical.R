#' Clinical variables of 229 patients
#'
#' A tibble of clinical characteristics for 229 patients, one row per patient.
#'
#' `deidentified_id` is **integer** here and in [example_summary], while the
#' `deidentified_sample` column of [example_spatial] is character.
#'
#' @format A tibble with 229 rows and 6 variables:
#' \describe{
#'   \item{age}{age at diagnosis}
#'   \item{race}{self-identified race}
#'   \item{sex}{patient biological sex}
#'   \item{status}{disease status}
#'   \item{deidentified_sample}{sample identifier}
#'   \item{deidentified_id}{patient identifier}
#' }
#' @seealso [example_spatial], [example_summary]
"example_clinical"
