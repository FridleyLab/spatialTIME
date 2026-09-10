# Shared plumbing for nearest-neighbour G statistics.
#
# Unlike Ripley's K, G has no pair-list-reuse trick available: G(r) is built from
# each cell's distance to its NEAREST neighbour, and the nearest neighbour changes
# whenever the point set changes. So permutations genuinely have to recompute.
# That is fine -- spatstat's Gest/Gcross use a kd-tree and are O(n log n) -- and it
# is why these functions delegate to spatstat rather than hand-rolling the
# estimator.
#
# Delegating is also a straight upgrade in correctness. The pre-2.0.0 bivariate
# implementation hand-rolled the rs and han estimators on a full
# as.matrix(dist(...)), i.e. O(n^2) memory, via bdist.points() + km.rs() +
# eroded.areas(). It agreed with Gcross exactly on well-populated markers (verified
# max|diff| = 0 for both rs and han) but returned NaN where Gcross returns 0 on
# sparse markers, and its "rs" branch referenced an undefined variable `W`. That
# last one never errored only because R's lazy evaluation meant the argument was
# never forced -- handle.r.b.args() ignores its window argument when `r` is
# supplied. Deleting ~150 lines in favour of one Gcross() call removes all of it.
#
# One trap worth recording: Gcross(correction = "rs") returns BOTH an `rs` and a
# `km` column, so selecting the estimate by column POSITION silently returns the
# wrong estimator for correction = "km". Always select by name -- that is what
# g_column() is for.


#' Normalise a nearest-neighbour edge correction name
#'
#' Pre-2.0.0 code accepted the misspelling "hans" (in `bi_NN_G_sample()` and the
#' `compute_metrics()` docs) and passed it straight to spatstat, which only knows
#' "han". Reject unknown values with a message naming the valid ones.
#'
#' @keywords internal
#' @noRd
match_g_correction <- function(x) {
  if (length(x) != 1L || !is.character(x)) {
    stop("`edge_correction` must be a single string.", call. = FALSE)
  }
  switch(x,
    rs = , border = , reduced = "rs",
    km = , kaplan = "km",
    han = , hans = , Hanisch = , hanisch = "han",
    none = , raw = "none",
    stop("Unsupported `edge_correction` for nearest-neighbour G: \"", x, "\". ",
         "Use one of \"rs\", \"km\", \"han\" or \"none\".", call. = FALSE)
  )
}


#' spatstat's output column name for a given G correction
#' @keywords internal
#' @noRd
g_column <- function(correction) {
  switch(correction, rs = "rs", km = "km", han = "han", none = "raw")
}


#' Univariate G for a subset of cells, using the whole-sample window
#'
#' @param pp `ppp` of **all** cells in the sample, carrying the sample window.
#' @param keep logical vector over all cells selecting the marker-positive ones.
#' @param r_range radii, including 0.
#' @param correction normalised correction name.
#' @return list with `theo` and `est`, each of length `length(r_range)`.
#' @keywords internal
#' @noRd
g_univariate <- function(pp, keep, r_range, correction) {
  # Subsetting a ppp keeps its window, so the observation window remains the
  # convex hull of every cell in the sample rather than of this marker's cells.
  out <- as.data.frame(spatstat.explore::Gest(
    pp[keep], r = r_range, correction = correction
  ))
  list(theo = out$theo, est = out[[g_column(correction)]])
}


#' Bivariate (cross-type) G, using the whole-sample window
#'
#' @param pp `ppp` of **all** cells in the sample, carrying the sample window.
#' @param keep_i,keep_j disjoint logical vectors over all cells giving the anchor
#'   and counted sets.
#' @param r_range radii, including 0.
#' @param correction normalised correction name.
#' @return list with `theo` and `est`.
#' @keywords internal
#' @noRd
g_bivariate <- function(pp, keep_i, keep_j, r_range, correction) {
  keep <- keep_i | keep_j
  Y <- pp[keep]
  spatstat.geom::marks(Y) <- factor(ifelse(keep_i[keep], "i", "j"), levels = c("i", "j"))
  out <- as.data.frame(spatstat.explore::Gcross(
    Y, i = "i", j = "j", r = r_range, correction = correction
  ))
  list(theo = out$theo, est = out[[g_column(correction)]])
}
