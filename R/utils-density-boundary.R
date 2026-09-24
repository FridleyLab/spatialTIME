# Internals for split_tissue() / plot_tissue_split().
#
# Two bugs in the SplittingTissue/ prototype this file fixes once, in one place,
# instead of at every call site:
#   1. contourLines(im$xcol, im$yrow, im$v) silently transposes the field -- an
#      `im` stores v with dim = c(length(yrow), length(xcol)) and contourLines()
#      has no dimension check. zero_contour() passes t(im$v).
#   2. `dimyx` is exposed and can give non-square pixels on a non-square window.
#      density_pixel_size() derives a square eps from sigma instead.
#
# Per-image [0, 1] rescaling is NOT a bug -- it is the `rescale` argument, on by
# default -- but it is a real tradeoff and the measurement behind it belongs here.
# It is not a shared monotone transform: each image gets its own affine map, so
# the zero set moves rather than being reparameterised. Measured on
# example_spatial[["TMA3_[9,K].tif"]] at sigma 40, the raw difference gives 9
# contour pieces / length 6756.5 and the rescaled difference 7 / 6723.8 -- two
# components gone. Transforms that genuinely are shared and monotone leave the
# zero set alone (edge=FALSE -> 6756.6, log(l1)-log(l2) -> 6756.4, relative risk
# -> 6756.4). The consequence to keep in mind: under `rescale` the boundary sits
# where each sample's *rescaled* densities are equal, so it depends on that
# sample's own density range and `Boundary Length` stops being an absolute
# cross-sample criterion. `rescale = FALSE` restores the absolute one.

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


#' Rescale an image's values to [0, 1], leaving a zero-range image untouched
#'
#' The guard is required, not defensive. `scales::rescale()` maps a zero-range
#' input to `mean(to)`, i.e. **0.5** -- so an absent class, whose density image
#' is all zero, would come back as a constant 0.5 and `r1 - r2` would cross zero
#' wherever `r1 = 0.5`, manufacturing a boundary where there is none. Returning
#' such an image unchanged keeps it at 0, so the difference stays one-signed and
#' the contour stays empty.
#'
#' @keywords internal
#' @noRd
rescale01 <- function(v) {
  r <- range(v, na.rm = TRUE)
  if (!all(is.finite(r)) || diff(r) == 0) return(v)
  (v - r[1]) / diff(r)
}


#' Unfiltered and filtered class1-minus-class2 density difference
#'
#' The **filtered** difference drives the boundary (contour / length / interface
#' distances); the **unfiltered** difference drives each cell's sign and
#' `density_score`. That split means masking a tissue hole can stop it inflating
#' the boundary without ever orphaning a cell to `NA`.
#'
#' Order of operations is load-bearing:
#' \enumerate{
#'   \item kernel densities per class;
#'   \item `min_density` mask, on **raw** intensities -- "is there tissue here?"
#'     is a question about absolute cell density, so it must be asked before any
#'     rescaling destroys the scale;
#'   \item `rescale` each image to [0, 1];
#'   \item difference -> the unfiltered field;
#'   \item `filter_density` per image, **after** rescale, so a filter written
#'     against [0, 1] values behaves as written;
#'   \item `hard_threshold` on the filtered image **only** -- applying it to the
#'     unfiltered field would flatten `density_score` to +/-1 and the plot raster
#'     with it.
#' }
#'
#' @param pp full-sample `ppp`, window = convex hull of every cell
#' @param keep1,keep2 logical masks (see [class_mask()]) selecting class1/class2
#' @param filter_density `NULL`, or a `function(im) im` applied to each class
#'   image before differencing
#' @param rescale rescale each class image to [0, 1] before differencing
#' @param min_density `NULL`/`NA` for off, else a fraction of this sample's mean
#'   `lambda_class1 + lambda_class2`. Pixels below it become `NA` in the contour
#'   field only -- never in `raw`, so no cell is orphaned by masking
#' @param hard_threshold collapse the filtered field to `sign()` before contouring
#' @return `list(raw = im, filtered = im)`
#' @keywords internal
#' @noRd
compartment_diff <- function(pp, keep1, keep2, sigma, eps, filter_density,
                             rescale = TRUE, min_density = NULL,
                             hard_threshold = FALSE) {
  d1 <- spatstat.explore::density.ppp(pp[keep1], sigma = sigma, eps = eps)
  d2 <- spatstat.explore::density.ppp(pp[keep2], sigma = sigma, eps = eps)

  # "Is there tissue here?" is a question about ABSOLUTE density, so compute the
  # mask from the raw intensities -- before any rescaling destroys the scale --
  # relative to this sample's own mean, so the same fraction means the same thing
  # across samples and coordinate units. It is only *applied* further down, to the
  # contour field, never to the field that drives each cell's sign.
  drop <- NULL
  if (!is.null(min_density) && !is.na(min_density) && min_density > 0) {
    tot  <- d1$v + d2$v
    lbar <- (sum(keep1) + sum(keep2)) /
      spatstat.geom::area(spatstat.geom::Window(pp))
    drop <- is.na(tot) | tot < min_density * lbar
  }

  # Rescaling is a *relative* comparison, so it needs two populated classes. With
  # one absent, skipping it is not a nicety: rescaling the present class to [0, 1]
  # puts its minimum at exactly 0, which IS the contour level, so a boundary would
  # appear around the edge of its support where the unrescaled field (strictly
  # positive, a sum of Gaussian kernels) correctly has none.
  if (isTRUE(rescale) && any(keep1) && any(keep2)) {
    d1$v <- rescale01(d1$v)
    d2$v <- rescale01(d2$v)
  }

  # Unmasked and unfiltered: drives each cell's sign and `density_score`. Masking
  # must never orphan a cell to NA, which is the same reason `filter_density` is
  # kept off this branch.
  raw <- d1 - d2

  c1 <- d1
  c2 <- d2
  if (!is.null(drop)) {
    c1$v[drop] <- NA
    c2$v[drop] <- NA
  }

  if (is.null(filter_density)) {
    filtered <- c1 - c2
  } else {
    f1 <- filter_density(c1)
    f2 <- filter_density(c2)
    validate_filtered_im(c1, f1)
    validate_filtered_im(c2, f2)
    filtered <- f1 - f2
  }

  if (isTRUE(hard_threshold)) {
    v <- filtered$v
    v[!is.na(v) & v > 0] <-  1
    v[!is.na(v) & v < 0] <- -1
    filtered$v <- v
  }

  list(raw = raw, filtered = filtered)
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


#' Density value at arbitrary points: bilinear, with a nearest-pixel fallback
#'
#' `interp.im()` interpolates bilinearly between the four surrounding pixel
#' centres rather than snapping to the nearest one, which matters here for a
#' reason beyond accuracy: `contourLines()` places the boundary by linear
#' interpolation along grid edges, so interpolating the field the same way keeps
#' `density_score` **consistent with the drawn polyline**. At the polyline's own
#' vertices this field is 0 to machine precision (<= 4.3e-15 across the example
#' cores); nearest-pixel lookup reads up to 0.041 there, ~4% of the rescaled
#' field's range.
#'
#' The payoff is the score, not the labels. Switching from nearest-pixel changes
#' the sign for only 6-12 cells per example core, and because those all sit within
#' half a pixel of the contour -- hence inside any sensible interface band -- the
#' 3-level label moves for 0-1 cells per core. What does change materially is
#' `density_score` itself, by up to 0.3 near the boundary.
#'
#' It returns `NA` in two situations, both of which the fallback has to cover:
#' any of the four neighbours is outside the mask (a cell inside the hull but
#' near a masked hole), or the point is outside the image frame entirely. Verified
#' that `interp.im()` returns `NA` rather than erroring in both cases. Falling
#' back to `nearest.valid.pixel()` means no cell is ever left `NA` by a lookup --
#' 178 cells on that same sample -- so an `NA` compartment can only ever mean an
#' exactly-zero field, which is what it is documented to mean.
#'
#' @keywords internal
#' @noRd
im_value_at <- function(im, x, y) {
  v <- spatstat.geom::interp.im(im, x, y)
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
                "rescale", "min_density", "hard_threshold",
                "sample_id", "spatialTIME_version")

  # These three arrived after the others, so a mif saved before they existed has
  # none of them -- and both returns below subset to `required`, which would make
  # every such mif fail the missing-field check. The backfill values are not
  # arbitrary defaults: they are the behaviour that predated the arguments, so an
  # older mif replots *correctly* rather than merely not erroring.
  backfill <- function(x) {
    if (is.null(x$rescale))        x$rescale        <- FALSE
    if (is.null(x$hard_threshold)) x$hard_threshold <- FALSE
    if (is.null(x$min_density))    x$min_density    <- NA_real_
    x
  }

  if (!is.null(settings)) {
    settings <- backfill(settings)
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
  call_info <- backfill(call_info)
  missing <- setdiff(required, names(call_info))
  if (length(missing)) {
    stop("`mif$derived$density_boundary`'s provenance is missing field",
         if (length(missing) > 1) "s" else "", ": ",
         paste0("`", missing, "`", collapse = ", "), ".", call. = FALSE)
  }
  call_info[required]
}
