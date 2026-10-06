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
    "Permutation p-value",
    "Degree of Clustering Theoretical",
    "Degree of Clustering Permutation",
    "Degree of Clustering Exact")
}


#' Summarise a permutation distribution into a count and a p-value
#'
#' Both quantities come from the same comparison, so they are derived together in
#' one place rather than per metric.
#'
#' `Permutations Larger than Observed` counts permutations at least as extreme as
#' the observation (`>=`). Ties count toward the count, and therefore raise the
#' p-value: with a heavily tied permutation distribution -- which happens at small
#' radii, where many relabellings produce the identical statistic -- a strict `>`
#' silently treats every tie as evidence of clustering. Before 2.0.0 this column
#' existed in only four of the seven metrics and used strict `>`.
#'
#' `Permutation p-value` is the standard Monte Carlo p-value,
#' `(1 + #{perm >= obs}) / (B + 1)`. The `+1`s are what make it a valid p-value for
#' any number of permutations: it can never be exactly 0, which a raw `count / B`
#' can, and which reads as infinite significance rather than "no permutation in
#' this sample was more extreme".
#'
#' `B` is the number of permutations that actually produced a value at that radius,
#' not `num_permutations`. Those differ wherever the estimator returns NA -- past
#' `rmax_valid` for Ripley's K, for instance -- and dividing by the requested count
#' there would understate the p-value.
#'
#' Where the observation itself is NA, both columns are NA. Previously the count
#' was computed with `rowSums(..., na.rm = TRUE)`, which returned `0` when every
#' term was NA and so reported "no permutation exceeded the observation" -- i.e.
#' maximal clustering -- for radii where nothing had been estimated at all.
#'
#' @param permuted numeric matrix of permuted statistics, radii down the rows and
#'   permutations across the columns, as the `vapply(..., numeric(n_r))` in every
#'   metric produces. Reshaped defensively, because `vapply` returns a bare vector
#'   when there is only one radius.
#' @param observed numeric vector of observed statistics, one per radius.
#' @return list with `larger` (integer), `p_value` (numeric) and `n_perm`
#'   (integer), each of length `length(observed)`.
#' @keywords internal
#' @noRd
permutation_summary <- function(permuted, observed) {
  n_r <- length(observed)
  permuted <- matrix(as.numeric(permuted), nrow = n_r)

  usable <- !is.na(permuted)
  n_perm <- rowSums(usable)
  at_least <- rowSums(permuted >= observed & usable, na.rm = TRUE)

  larger <- as.integer(at_least)
  p_value <- (1 + at_least) / (n_perm + 1)

  unusable <- is.na(observed) | n_perm == 0L
  larger[unusable] <- NA_integer_
  p_value[unusable] <- NA_real_

  list(larger = larger, p_value = p_value, n_perm = as.integer(n_perm))
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
