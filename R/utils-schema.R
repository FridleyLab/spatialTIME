# One output schema for every spatial metric in the package.
#
# Before 2.0.0 each metric invented its own column names for the same concepts:
# Theoretical CSR / Theoretical G / Theoretical g; Permuted CSR / Permuted G /
# Permuted g / Permuted Interaction; Anchor+Counted / From+To; Degree of
# Clustering / Degree of Correlation / Degree of Interaction; and
# Permuted_larger_than_Observed vs "Permutations Larger than Observed". Anything
# wanting to treat two metrics uniformly -- plotting, joining to clinical data,
# feeding a model -- had to special-case each one.
#
# 2.0.0 settles on the schema below. Only the "Observed" column keeps a
# statistic-specific name, because that is the one place the distinction is
# informative: Observed K, Observed G, Observed g, Observed Interaction.
#
# Columns a given metric cannot fill are present and NA rather than absent, so
# that dplyr::bind_rows() across metrics lines up instead of producing a ragged
# frame. In particular `Exact CSR` is NA for everything except Ripley's K -- see
# the "Why there is no exact CSR for G" section of ?NN_G.


#' Canonical column order for a metric results table
#'
#' @param sample_id the mif's sample id column name.
#' @param observed name of the statistic-specific observed column.
#' @param bivariate `TRUE` for Anchor/Counted, `FALSE` for Marker.
#' @return character vector of column names in canonical order.
#' @keywords internal
#' @noRd
standard_metric_cols <- function(sample_id, observed, bivariate) {
  c(sample_id,
    if (bivariate) c("Anchor", "Counted") else "Marker",
    "iter",
    "r",
    "Theoretical CSR",
    "Permuted CSR",
    "Exact CSR",
    observed,
    "Permutations Larger than Observed",
    "Degree of Clustering Theoretical",
    "Degree of Clustering Permutation",
    "Degree of Clustering Exact")
}


#' Put a results table into canonical shape
#'
#' Adds any missing schema columns as NA, drops nothing, and orders columns
#' canonically with any extra columns kept at the end rather than silently
#' discarded.
#'
#' @param out results data frame.
#' @param sample_id,observed,bivariate as for [standard_metric_cols()].
#' @keywords internal
#' @noRd
as_standard_metric <- function(out, sample_id, observed, bivariate) {
  want <- standard_metric_cols(sample_id, observed, bivariate)
  for (nm in setdiff(want, names(out))) out[[nm]] <- NA
  extra <- setdiff(names(out), want)
  out[, c(want, extra), drop = FALSE]
}


#' Fill in the three Degree of Clustering columns
#'
#' Every metric derives these the same way -- observed minus each reference -- so
#' the arithmetic lives here rather than being repeated per function.
#'
#' @keywords internal
#' @noRd
add_degrees_of_clustering <- function(out, observed) {
  obs <- out[[observed]]
  out[["Degree of Clustering Theoretical"]] <- obs - out[["Theoretical CSR"]]
  out[["Degree of Clustering Permutation"]] <- obs - out[["Permuted CSR"]]
  out[["Degree of Clustering Exact"]]       <- obs - out[["Exact CSR"]]
  out
}
