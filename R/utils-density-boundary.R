# Internals for split_tissue() / plot_tissue_split().
#
# The prototype in SplittingTissue/functions.R has three bugs this file exists to
# fix once, in one place, instead of at every call site:
#   1. contourLines(im$xcol, im$yrow, im$v) silently transposes the field -- an
#      `im` stores v with dim = c(length(yrow), length(xcol)) and contourLines()
#      has no dimension check. zero_contour() passes t(im$v).
#   2. Rescaling each class image to [0, 1] before differencing is not a shared
#      monotone transform (each image gets a different affine map), so it moves
#      and destroys contour components. compartment_diff() differences the raw
#      KDE intensities instead.
#   3. `dimyx` is exposed and can give non-square pixels on a non-square window.
#      density_pixel_size() derives a square eps from sigma instead.

#' @keywords internal
#' @noRd
DENSITY_EPS_DIVISOR <- 8L

#' @keywords internal
#' @noRd
DENSITY_MAX_PIXELS <- 4e7


#' Square pixel size for a given bandwidth
#'
#' Resolution is derived from `sigma` rather than exposed as `dimyx`: contour
#' *topology* is set by sigma, resolution only adds a converging length bias
#' (see `?split_tissue`), and deriving `eps` from `sigma` guarantees square
#' pixels regardless of the window's aspect ratio, which `dimyx` does not.
#'
#' @keywords internal
#' @noRd
density_pixel_size <- function(sigma) {
  sigma / DENSITY_EPS_DIVISOR
}


#' Error out before an OOM instead of after
#'
#' Hiding `dimyx` makes `eps` (and therefore `sigma`) the only lever on the
#' density grid's memory, and that relationship is quadratic and invisible to
#' the caller. Without this check a small `sigma` on a large window is an OOM
#' kill rather than an informative error.
#'
#' @keywords internal
#' @noRd
check_density_budget <- function(win, eps, sigma, label) {
  a <- spatstat.geom::area(spatstat.geom::Frame(win))
  n <- a / eps^2
  if (n > DENSITY_MAX_PIXELS) {
    sigma_min <- DENSITY_EPS_DIVISOR * sqrt(a / DENSITY_MAX_PIXELS)
    stop(sprintf(
      paste0("Sample \"%s\" needs a %s-pixel density grid at sigma = %g.\n",
             "  The internal pixel size is sigma/%d, so halving sigma quadruples memory.\n",
             "  Use sigma >= %.0f for this sample, or split the image first."),
      label, format(round(n), big.mark = ","), sigma, DENSITY_EPS_DIVISOR,
      ceiling(sigma_min)), call. = FALSE)
  }
  invisible(n)
}


#' Logical mask for one classifier level, NA-safe
#'
#' `values == level` propagates `NA`, and a logical subscript containing `NA`
#' makes `pp[keep]` error ("Index out of bounds in [.ppp"). Unclassified cells
#' are realistic, so they must be masked out rather than erroring.
#'
#' @keywords internal
#' @noRd
class_mask <- function(values, level) {
  !is.na(values) & values == level
}


#' Raw and (optionally) filtered class1-minus-class2 density difference
#'
#' The **filtered** difference drives the boundary (contour / length /
#' interface distances); the **unfiltered** difference drives each cell's sign.
#' That split means a `filter_density` that masks a tissue hole can stop it
#' inflating the boundary without ever orphaning a cell to `NA`.
#'
#' @param pp full-sample `ppp`, window = convex hull of every cell
#' @param keep1,keep2 logical masks (see [class_mask()]) selecting class1/class2
#' @param filter_density `NULL`, or a `function(im) im` applied to each class
#'   image before differencing
#' @return `list(raw = im, filtered = im)`
#' @keywords internal
#' @noRd
compartment_diff <- function(pp, keep1, keep2, sigma, eps, filter_density) {
  d1 <- spatstat.explore::density.ppp(pp[keep1], sigma = sigma, eps = eps)
  d2 <- spatstat.explore::density.ppp(pp[keep2], sigma = sigma, eps = eps)
  raw <- d1 - d2

  if (is.null(filter_density)) {
    return(list(raw = raw, filtered = raw))
  }

  f1 <- filter_density(d1)
  f2 <- filter_density(d2)
  validate_filtered_im(d1, f1)
  validate_filtered_im(d2, f2)
  list(raw = raw, filtered = f1 - f2)
}


#' Guard that `filter_density` did not change the pixel grid
#'
#' `d1 - d2` (and, downstream, `f1 - f2`) is only meaningful when both images
#' share a grid. A filter that returns a coarser/finer `im`, or something that
#' is not an `im` at all, must error here rather than misalign silently.
#'
#' @keywords internal
#' @noRd
validate_filtered_im <- function(before, after) {
  if (!inherits(after, "im")) {
    stop("`filter_density` must return an object of class \"im\", ",
         "got \"", class(after)[1], "\".", call. = FALSE)
  }
  if (!identical(dim(before$v), dim(after$v)) ||
      !isTRUE(all.equal(before$xcol, after$xcol)) ||
      !isTRUE(all.equal(before$yrow, after$yrow))) {
    stop("`filter_density` must return an image on the same pixel grid it ",
         "was given -- it may recolour pixels (e.g. set some to NA) but not ",
         "resize or resample the grid.", call. = FALSE)
  }
  invisible(TRUE)
}


#' Exact zero level set of a density-difference image
#'
#' `t(im$v)` is load-bearing -- see the file header. Pieces with fewer than 2
#' vertices (an isolated saddle point) carry no length and are dropped.
#'
#' @return `data.frame(<sample_id>, piece, x, y)`, 0 rows if no contour exists
#' @keywords internal
#' @noRd
zero_contour <- function(im, sample_id, label) {
  cs <- withCallingHandlers(
    grDevices::contourLines(x = im$xcol, y = im$yrow, z = t(im$v), levels = 0),
    warning = function(w) if (grepl("all z values are NA", conditionMessage(w)))
      invokeRestart("muffleWarning"))
  cs <- cs[vapply(cs, function(c) length(c$x) >= 2L, logical(1))]
  if (!length(cs)) return(empty_boundary_df(sample_id))
  out <- do.call(rbind, lapply(seq_along(cs), function(i)
    data.frame(label, piece = i, x = cs[[i]]$x, y = cs[[i]]$y,
               stringsAsFactors = FALSE)))
  names(out)[1] <- sample_id
  out$piece <- as.integer(out$piece)
  out
}


#' Zero-row boundary frame with correct column types
#'
#' @keywords internal
#' @noRd
empty_boundary_df <- function(sample_id) {
  out <- data.frame(character(0), integer(0), numeric(0), numeric(0),
                    stringsAsFactors = FALSE)
  names(out) <- c(sample_id, "piece", "x", "y")
  out
}


#' Boundary polyline as a spatstat line-segment pattern
#'
#' Each piece's vertices are already ordered along the polyline by
#' `contourLines()`, so consecutive rows within a piece become one segment.
#' 0-row input gives a 0-segment `psp` rather than an error.
#'
#' @keywords internal
#' @noRd
boundary_psp <- function(boundary_df, win) {
  if (!nrow(boundary_df)) {
    return(spatstat.geom::psp(numeric(0), numeric(0), numeric(0), numeric(0),
                              window = spatstat.geom::Frame(win)))
  }
  segs <- do.call(rbind, lapply(split(boundary_df, boundary_df$piece), function(p) {
    if (nrow(p) < 2L) return(NULL)
    data.frame(x0 = p$x[-nrow(p)], y0 = p$y[-nrow(p)],
               x1 = p$x[-1],       y1 = p$y[-1])
  }))
  if (is.null(segs) || !nrow(segs)) {
    return(spatstat.geom::psp(numeric(0), numeric(0), numeric(0), numeric(0),
                              window = spatstat.geom::Frame(win)))
  }
  spatstat.geom::psp(segs$x0, segs$y0, segs$x1, segs$y1,
                     window = spatstat.geom::Frame(win), check = FALSE)
}


#' Total boundary length
#'
#' @keywords internal
#' @noRd
boundary_length <- function(S) {
  sum(spatstat.geom::lengths_psp(S))
}


#' Density value at arbitrary points, with an edge-pixel fallback
#'
#' `lookup.im()` returns `NA` for cells whose pixel *centre* falls outside the
#' image mask even though the cell itself is inside the hull the mask was built
#' from (verified: 9-16 per example core). `nearest.valid.pixel()` resolves
#' those instead of leaving the cell's compartment as `NA`.
#'
#' @keywords internal
#' @noRd
im_value_at <- function(im, x, y) {
  v <- spatstat.geom::lookup.im(im, x, y, naok = TRUE)
  na <- is.na(v)
  if (any(na)) {
    nv <- spatstat.geom::nearest.valid.pixel(x[na], y[na], im)
    v[na] <- im$v[cbind(nv$row, nv$col)]
  }
  v
}


#' Retrieve `split_tissue()` provenance, from `settings` or from the slot's attribute
#'
#' `attr()` on a list survives `[[` access and `saveRDS()`, but is dropped by
#' `[` subsetting and by `dplyr::bind_rows()`. `settings` is the fallback for
#' whenever that attribute has gone missing.
#'
#' @keywords internal
#' @noRd
split_tissue_settings <- function(mif, settings = NULL) {
  required <- c("classifier", "class1", "class2", "sigma", "eps",
                "interface_width", "xloc", "yloc", "filter_density",
                "sample_id", "spatialTIME_version")

  if (!is.null(settings)) {
    missing <- setdiff(required, names(settings))
    if (length(missing)) {
      stop("`settings` is missing required field", if (length(missing) > 1) "s" else "",
           ": ", paste0("`", missing, "`", collapse = ", "), ".", call. = FALSE)
    }
    return(settings[required])
  }

  boundary <- mif$derived$density_boundary
  if (is.null(boundary)) {
    stop("`mif$derived$density_boundary` not found -- run `split_tissue()` first, ",
         "or pass `settings` explicitly.", call. = FALSE)
  }
  call_info <- attr(boundary, "call_info")
  if (is.null(call_info)) {
    stop("The provenance attached to `mif$derived$density_boundary` by ",
         "`split_tissue()` is missing (it does not survive `[` subsetting or ",
         "`merge_mifs()`). Re-run `split_tissue()`, or pass `settings` explicitly ",
         "-- see `?split_tissue_settings`.", call. = FALSE)
  }
  missing <- setdiff(required, names(call_info))
  if (length(missing)) {
    stop("`mif$derived$density_boundary`'s provenance is missing field",
         if (length(missing) > 1) "s" else "", ": ",
         paste0("`", missing, "`", collapse = ", "), ".", call. = FALSE)
  }
  call_info[required]
}
