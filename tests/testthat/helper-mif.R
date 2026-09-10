# Shared fixtures.
#
# Four test files each rebuilt the same mif at the top of the file before this
# existed. Everything here is a function rather than a top-level object so that a
# test which mutates a mif cannot leak that mutation into the next test -- which
# matters a lot here, because most functions in this package take a mif and return
# a modified copy.

#' A small mif built from the shipped example data
#'
#' @param which which of the five example cores to include (index or name).
#' @param n_cells thin each spatial frame to at most this many cells. The full
#'   cores are 1800-3800 cells; a few hundred keeps permutation tests fast without
#'   changing what is being tested.
#' @param seed thinning is deterministic so failures are reproducible.
example_mif <- function(which = "TMA3_[9,K].tif", n_cells = 400, seed = 1) {
  sp <- example_spatial[which]
  if (!is.null(n_cells)) {
    set.seed(seed)
    sp <- lapply(sp, function(x) if (nrow(x) > n_cells) x[sort(sample(nrow(x), n_cells)), ] else x)
  }
  ids <- names(sp)
  create_mif(
    clinical_data = data.frame(deidentified_id = as.character(seq_along(ids)),
                               stringsAsFactors = FALSE),
    sample_data   = data.frame(deidentified_id = as.character(seq_along(ids)),
                               deidentified_sample = ids, stringsAsFactors = FALSE),
    spatial_list  = sp,
    patient_id = "deidentified_id",
    sample_id  = "deidentified_sample"
  )
}

#' The full-size mif, for tests that need real marker counts
#'
#' `TMA3_[9,K].tif` is the only core with substantial positives (536 CD3+, 83
#' CD8+), so anything comparing against spatstat on non-degenerate data wants this
#' rather than a thinned copy.
example_mif_full <- function(which = "TMA3_[9,K].tif") {
  example_mif(which = which, n_cells = NULL)
}

#' Markers that exist in every example core
#'
#' Ordered most- to least-abundant. Fine for univariate measures.
mnames_good <- function() {
  c("CD3..Opal.570..Positive", "CD8..Opal.520..Positive",
    "FOXP3..Opal.620..Positive", "CD3..CD8.", "CD3..FOXP3.")
}

#' A marker pair that actually works for bivariate measures on the example data
#'
#' Do NOT use mnames_good()[1:2] (CD3 / CD8) for anything bivariate. Every
#' bivariate measure in this package drops cells positive for both markers of a
#' pair, and in this data CD8 is almost perfectly nested inside CD3 -- 16 of the 17
#' CD8+ cells in a 400-cell thinning are also CD3+, leaving ONE counted cell. That
#' is biologically correct (CD8 T cells are CD3+) but it makes the pair useless as a
#' fixture, and the failure looks like a bug in the function under test.
#'
#' CD8 / FOXP3 leaves 79 and 30 disjoint cells on the full `TMA3_[9,K].tif`, which
#' is the best real-data pair available. For tests that care about bookkeeping
#' rather than statistics, prefer [toy_spatial()], which guarantees disjoint counts.
mnames_bivariate <- function() {
  c("CD8..Opal.520..Positive", "FOXP3..Opal.620..Positive")
}

#' A sample whose markers sit in opposite corners of a large window
#'
#' Used to prove the window invariant: the convex hull of all cells is ~8x the area
#' of the hull of any one marker's cells, so a result computed against a
#' subset-derived window is off by that ratio and cannot be mistaken for noise.
#' Shared by test-window-invariant.R and the per-function tests.
corner_marker_mif <- function(seed = 6, n = 3000) {
  set.seed(seed)
  spat <- data.frame(
    deidentified_sample = "S1",
    XMin = runif(n, 0, 1000),
    YMin = runif(n, 0, 1000)
  )
  spat$XMax <- spat$XMin
  spat$YMax <- spat$YMin
  spat$cornerA <- as.integer(spat$XMin < 350 & spat$YMin < 350)
  spat$cornerB <- as.integer(spat$XMin > 650 & spat$YMin > 650)
  spat$spread  <- as.integer(seq_len(n) %% 7 == 0)
  # Two classifier levels, unbalanced, so subset_mif/marker_freq_diff can use it.
  spat$Classifier.Label <- ifelse(spat$XMin < 300, "Stroma", "Tumor")
  list(
    spat = spat,
    mif = create_mif(
      clinical_data = data.frame(deidentified_id = "p1", stringsAsFactors = FALSE),
      sample_data   = data.frame(deidentified_id = "p1", deidentified_sample = "S1",
                                 stringsAsFactors = FALSE),
      spatial_list  = list(S1 = spat),
      patient_id = "deidentified_id",
      sample_id  = "deidentified_sample"
    )
  )
}

#' Hand-built spatial data with full control over counts and levels
#'
#' For the cases the shipped data cannot express: a classifier level missing from
#' one sample, a level with <= 2 cells, more than two classifier levels, marker
#' names containing spaces, or several samples per patient.
toy_mif <- function(spatial_list,
                    patient_ids = NULL,
                    patient_id = "deidentified_id",
                    sample_id = "deidentified_sample") {
  ids <- names(spatial_list)
  if (is.null(patient_ids)) patient_ids <- as.character(seq_along(ids))
  create_mif(
    clinical_data = data.frame(x = patient_ids, stringsAsFactors = FALSE) |>
      stats::setNames(patient_id),
    sample_data = stats::setNames(
      data.frame(patient_ids, ids, stringsAsFactors = FALSE),
      c(patient_id, sample_id)
    ),
    spatial_list = spatial_list,
    patient_id = patient_id,
    sample_id  = sample_id
  )
}

#' A minimal spatial frame with exactly the requested marker counts
#'
#' @param n total cells.
#' @param markers named integer vector: how many cells are positive for each.
#' @param levels classifier levels to cycle through.
toy_spatial <- function(name = "S1", n = 200, markers = c(A = 40, B = 40),
                        levels = c("Tumor", "Stroma"), seed = 11) {
  set.seed(seed)
  d <- data.frame(
    deidentified_sample = name,
    XMin = runif(n, 0, 500), YMin = runif(n, 0, 500),
    stringsAsFactors = FALSE
  )
  d$XMax <- d$XMin
  d$YMax <- d$YMin
  d$Classifier.Label <- rep(levels, length.out = n)
  # Disjoint positive sets so bivariate measures have clean anchor/counted groups.
  taken <- integer(0)
  for (m in names(markers)) {
    k <- markers[[m]]
    pool <- setdiff(seq_len(n), taken)
    idx <- if (k > 0) sort(sample(pool, min(k, length(pool)))) else integer(0)
    taken <- c(taken, idx)
    d[[m]] <- 0L
    d[[m]][idx] <- 1L
  }
  d
}
