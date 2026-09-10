# Exact, memory-bounded engine for count-based (Ripley's K) statistics.
#
# Why this exists
# ---------------
# The naive way to compute K on a large sample is to build an n x n distance
# matrix and an n x n edge-correction matrix. For a whole-slide image that is
# fatal: at n = 200,000 an n x n double matrix is ~319 GB. Earlier versions of
# this package worked around that by tiling the distance matrix, which bounded
# memory but cost accuracy -- above a cell-count threshold the requested edge
# correction was silently replaced with "none", and the tiled branch hard-coded
# translation regardless of what the caller asked for.
#
# None of that is necessary. spatstat.geom::closepairs(X, rmax) returns only the
# pairs closer than rmax, and edge.Trans(paired = TRUE) / edge.Ripley(r = d)
# return one weight per pair rather than an n x n matrix. At n = 200,000 with a
# realistic rmax that is ~3.1M pairs (~75 MB) computed in a fraction of a second,
# and the result is exact -- verified equal to Kest() to ~1e-20.
#
# The second thing this buys us: for a random-labelling null, the translation and
# isotropic edge weight of a pair depends only on that pair's displacement and on
# the window, never on which other points are in the pattern. So the pair list and
# its weights can be computed ONCE per sample and reused for every marker, every
# marker pair, and every permutation -- each of which then costs only a logical
# mask plus a weighted histogram. Measured ~26x faster than re-calling Kcross per
# permutation, and equal to it to ~1e-17.
#
# Fidelity to spatstat
# --------------------
# These functions deliberately mirror the internals of spatstat.explore::Kest and
# ::Kmulti rather than re-deriving the estimator:
#
#   breaks <- handle.r.b.args(r, NULL, W, rmaxdefault)   # same bin edges
#   wh     <- whist(d, breaks$val, edgeweights)          # same weighted histogram
#   K      <- cumsum(wh) / lambda2area                   # same normalisation
#   K[r >= h] <- NA                                      # same rmax truncation
#
# Calling the same primitives in the same order is what makes agreement exact
# rather than approximate, and it means bin-edge conventions never have to be
# guessed at. tests/testthat/test-k-engine.R pins this against Kest/Kcross for
# every correction, on both continuous and integer-valued coordinates, because
# integer coordinates (which is what HALO and Vectra actually emit) are where
# binning conventions actually bite.
#
# Two spatstat behaviours worth knowing about, both discovered by testing rather
# than by reading:
#
# 1. Kest is not internally consistent about tied distances, and which answer you
#    get depends on what else you asked for. When "none" or "border" is the ONLY
#    correction requested and r is evenly spaced, Kest takes a fast C path
#    (will.do.fast) that bins right-closed, d <= r. Request "none" alongside any
#    non-fast correction and the same call routes through whist() instead, which
#    bins left-closed. Those disagree whenever a pair distance lands exactly on a
#    bin edge -- never on continuous coordinates, but constantly on the integer
#    coordinates HALO and Vectra emit. On a 41x41 integer grid the two Kest paths
#    give 390.48 vs 0 at the same radius.
#
#    This engine always uses the whist path, which means every value it produces
#    is bit-identical to Kest(correction = c("none", "translation"))$un -- a
#    legitimate spatstat number, verified to max|diff| = 0. The benefit over
#    chasing the fast path is that all four corrections then share one binning
#    convention, so "none" and "translation" from this package are directly
#    comparable to each other. Kest's fast path is not comparable to Kest's own
#    translation output on tied data.
#
# 2. Isotropic weights on a polygonal window have genuine discontinuities. The
#    weight is 1 / (fraction of the circle of radius d that lies inside W). If
#    that circle passes exactly through a vertex of the window polygon, the
#    fraction jumps. Because the observation window here is the convex hull of the
#    cells, its vertices ARE cells, so a pair whose distance happens to equal the
#    anchor's distance to a hull vertex sits exactly on such a jump. spatstat's
#    Kest uses closepairs() and its Kmulti uses crosspairs(), and those two
#    compute the same distance with up to 1 ulp (1.4e-17) of difference -- enough
#    to land on opposite sides of the jump. Consequence: univariate isotropic
#    agrees with Kest exactly (both use closepairs), while bivariate isotropic can
#    differ from Kcross by ~1e-6 relative at a single radius when such a
#    degenerate pair exists. Verified against a 2-million-point Monte Carlo
#    integration of the true arc fraction: spatstat's crosspairs-side value is the
#    accurate one, and the discrepancy is confined to that one degenerate pair.
#    Translation -- the package default -- is exact for both.


#' Bin edges matching those spatstat would use
#'
#' @param r_range numeric vector of radii, including 0.
#' @param win observation window (`owin`).
#' @param lambda intensity used for spatstat's default rmax rule.
#' @return the `breakpts` object produced by [spatstat.geom::handle.r.b.args()].
#' @keywords internal
#' @noRd
k_breaks <- function(r_range, win, lambda) {
  rmaxdefault <- spatstat.explore::rmax.rule("K", win, lambda)
  if (is.infinite(rmaxdefault)) rmaxdefault <- spatstat.geom::diameter(win)
  spatstat.geom::handle.r.b.args(r_range, NULL, win, rmaxdefault = rmaxdefault)
}


#' Close pairs and their edge-correction weights for one sample
#'
#' Computed once per sample from **all** cells in that sample, never from a
#' marker-positive subset. Marker selection happens afterwards as a logical mask
#' over the returned pair list (see [k_from_pairs()]), which is what keeps the
#' observation window, its area and the edge weights identical across every
#' marker and every permutation within a sample.
#'
#' @param pp `ppp` object containing **every** cell in the sample.
#' @param r_range numeric vector of radii, including 0.
#' @param edge_correction one of "translation", "isotropic", "none", "border".
#' @param block integer; if the pair list would exceed roughly this many pairs,
#'   the weights are computed in chunks to bound peak memory. Results are
#'   identical either way -- weighted-histogram counts are additive.
#' @return list with `i`, `j`, `d` (pair indices into `pp` and their distances),
#'   `w` (per-pair edge weight), `rmax_valid` (radius at or beyond which spatstat
#'   returns NA for this correction), and the `breaks` object.
#' @keywords internal
#' @noRd
k_pairs <- function(pp, r_range, edge_correction, block = 5e6) {
  edge_correction <- match_edge_correction(edge_correction)
  win  <- spatstat.geom::Window(pp)
  n    <- spatstat.geom::npoints(pp)
  area <- spatstat.geom::area(win)
  rmax <- max(r_range)

  breaks <- k_breaks(r_range, win, lambda = n / area)

  cp <- spatstat.geom::closepairs(pp, rmax = rmax, what = "ijd")
  npair <- length(cp$d)

  # Per-pair edge weights, chunked so peak memory stays bounded even when the
  # pair list itself is large. Chunking is exact: each pair's weight depends only
  # on that pair.
  w <- rep_len(1, npair)
  rmax_valid <- Inf

  if (edge_correction == "translation" && npair > 0L) {
    gW <- spatstat.geom::setcov(win)
    idx <- chunk_index(npair, block)
    for (k in idx) {
      ew <- spatstat.explore::edge.Trans(
        X = spatstat.geom::ppp(pp$x[cp$i[k]], pp$y[cp$i[k]], window = win, check = FALSE),
        Y = spatstat.geom::ppp(pp$x[cp$j[k]], pp$y[cp$j[k]], window = win, check = FALSE),
        paired = TRUE, gW = gW, give.rmax = TRUE
      )
      w[k] <- as.numeric(ew)
      rmax_valid <- attr(ew, "rmax")
    }
    if (is.null(rmax_valid)) rmax_valid <- Inf
  } else if (edge_correction == "isotropic" && npair > 0L) {
    idx <- chunk_index(npair, block)
    for (k in idx) {
      w[k] <- as.numeric(spatstat.explore::edge.Ripley(
        X = spatstat.geom::ppp(pp$x[cp$i[k]], pp$y[cp$i[k]], window = win, check = FALSE),
        r = cp$d[k]
      ))
    }
    rmax_valid <- spatstat.geom::boundingradius(win)
  }

  list(i = cp$i, j = cp$j, d = cp$d, w = w,
       n = n, area = area, win = win, breaks = breaks,
       r_range = r_range, edge_correction = edge_correction,
       rmax_valid = rmax_valid)
}


#' K from a precomputed pair list
#'
#' @param pairs value of [k_pairs()].
#' @param keep_i,keep_j logical vectors of length `pairs$n` selecting which cells
#'   act as the "from" (i) and "to" (j) type. For the univariate case pass the
#'   same vector twice.
#' @param n_i,n_j number of cells of each type; taken from `keep_i`/`keep_j` when
#'   not supplied.
#' @return numeric vector of length `length(pairs$r_range)`.
#' @keywords internal
#' @noRd
k_from_pairs <- function(pairs, keep_i, keep_j = keep_i,
                         n_i = sum(keep_i), n_j = sum(keep_j),
                         univariate = identical(keep_i, keep_j)) {
  r <- pairs$r_range
  # Denominator matches spatstat: npts*(npts-1)/area for Kest, nI*nJ/area for
  # Kmulti/Kcross.
  denom <- if (univariate) n_i * (n_i - 1) else n_i * n_j
  if (is.na(denom) || denom <= 0) return(rep(NA_real_, length(r)))

  sel <- keep_i[pairs$i] & keep_j[pairs$j]
  if (!any(sel)) return(rep(0, length(r)))

  k <- k_cumsum(pairs$d[sel], pairs$w[sel], pairs) / (denom / pairs$area)
  k[r >= pairs$rmax_valid] <- NA_real_
  k
}


#' Cumulative weighted pair counts, binned the way spatstat bins them
#'
#' One convention for every correction: `spatstat.univar::whist()` over the bin
#' edges `spatstat.geom::handle.r.b.args()` produces, which is the path
#' `spatstat.explore::Kest()` and `::Kmulti()` take for translation and isotropic,
#' and the path they take for "none" whenever any other correction is also
#' requested. See note 1 in the header comment for why we do not chase Kest's
#' fast-path binning instead.
#'
#' @param d,w pair distances and their edge weights.
#' @param pairs value of [k_pairs()], for `breaks`.
#' @return numeric vector of length `length(pairs$r_range)`.
#' @keywords internal
#' @noRd
k_cumsum <- function(d, w, pairs) {
  cumsum(spatstat.univar::whist(d, pairs$breaks$val, w))
}


#' Normalise the many spellings of an edge correction
#'
#' `Kest` accepts several spellings but names its output columns differently
#' again; earlier versions of this package accepted a spelling and then silently
#' did something else, or left `edge` undefined and failed later with an unhelpful
#' "object 'edge' not found". Fail loudly and early instead.
#'
#' @keywords internal
#' @noRd
match_edge_correction <- function(x) {
  if (length(x) != 1L || !is.character(x)) {
    stop("`edge_correction` must be a single string.", call. = FALSE)
  }
  switch(x,
    translation = , trans = "translation",
    isotropic = , iso = , Ripley = "isotropic",
    none = , un = "none",
    border = , bord = "border",
    stop("Unsupported `edge_correction`: \"", x, "\". ",
         "Use one of \"translation\", \"isotropic\", \"border\" or \"none\".",
         call. = FALSE)
  )
}


#' Split 1:n into chunks of at most `size`
#' @keywords internal
#' @noRd
chunk_index <- function(n, size) {
  if (n == 0L) return(list())
  if (n <= size) return(list(seq_len(n)))
  split(seq_len(n), ceiling(seq_len(n) / size))
}
