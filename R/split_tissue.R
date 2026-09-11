#' Split a sample into tissue compartments by density difference
#'
#' @param mif object of class `mif` created with `create_mif()`
#' @param classifier column in each spatial file giving the per-cell tissue
#'   classification (e.g. Tumor/Stroma/Lymph/Necrotic). May have more than two
#'   levels; only `class1` and `class2` are compared.
#' @param class1,class2 the two classifier levels to compare. The density
#'   difference is `class1` minus `class2`, so a cell on the positive side of
#'   the boundary is labelled `class1` and a cell on the negative side is
#'   labelled `class2`. Neither may be named `"Interface"` (reserved) and they
#'   must differ.
#' @param sigma kernel density bandwidth, in the same coordinate units as the
#'   spatial data. There is no default: the boundary's shape depends on it and
#'   values that make sense differ by imaging platform, so a silent default
#'   would make cores or cohorts processed with different (undocumented)
#'   defaults incomparable. Memory is quadratic in `1/sigma` -- see `check_density_budget()`
#'   in `R/utils-density-boundary.R` and the pixel-budget error below.
#' @param interface_width full width (diameter, not radius) of the interface
#'   band, in the same units as `sigma`. Cells within `interface_width / 2` of
#'   the boundary are labelled `"Interface"` in `refined_density_compartment`.
#' @param workers number of cores to use for calculations
#' @param overwrite whether to replace `density_compartment`,
#'   `refined_density_compartment`, `mif$derived$density_boundary` and
#'   `mif$sample`'s `Boundary Length` column if any already exist. There is no
#'   `Run`-based append here (unlike the metric functions): a cell can carry
#'   only one compartment label, so run `split_tissue()` twice under different
#'   `sigma`/`interface_width` on two separate mifs if you want to compare them.
#' @param xloc,yloc columns giving the cell centre. If left `NULL`, `XMin`, `XMax`,
#'   `YMin` and `YMax` must be present and the centre is their midpoint.
#' @param ... `filter_density`, an optional `function(im) im` applied to each
#'   class's density image before differencing, to suppress the boundary
#'   inflating/bouncing around near-zero-density holes in the tissue (a hole is
#'   close to zero in *both* compartments, not just one). It affects only the
#'   boundary geometry (contour, length, interface distances) -- the sign that
#'   drives `density_compartment`/`refined_density_compartment` always comes
#'   from the *unfiltered* difference, so a filtered-out hole never leaves a
#'   cell `NA`. Must return an `im` on the same pixel grid it was given.
#'   Anything else in `...` is an error naming the offender.
#'
#' @description
#' For each spatial sample, `split_tissue()` computes kernel density estimates
#' of `class1` and `class2`, takes their difference, and extracts the *exact*
#' zero level set of that difference with [grDevices::contourLines()] -- not a
#' thresholded band around zero. Every cell (not just `class1`/`class2` cells)
#' is then labelled by which side of that boundary it falls on.
#'
#' @section What is and is not kept:
#' The density images and point patterns used to compute the boundary are
#' discarded once the boundary and per-cell labels are derived -- keeping them
#' would multiply the size of the mif by the number of samples. Only the
#' boundary polyline (`mif$derived$density_boundary`) and the two new spatial
#' columns survive. To visualise the density alongside the boundary, call
#' [plot_tissue_split()], which recomputes the density on demand at plot time.
#'
#' @section Choosing sigma:
#' `sigma` has no default. Its units are the same as your spatial coordinates,
#' which differ by imaging platform, so there is no value that is reasonable
#' for every dataset. Memory is quadratic in `1/sigma` (pixel size is
#' `sigma / 8` internally, so halving `sigma` quadruples the density grid) --
#' too small a `sigma` errors naming a minimum viable value for that sample
#' rather than exhausting memory.
#'
#' @section Boundary length:
#' Pixel resolution (`dimyx` in the underlying [spatstat.explore::density.ppp()]
#' call) is not exposed. Measured on a real core (Tumor vs Stroma, sigma 40,
#' 1803 cells): the number of contour pieces is 9 at every resolution from
#' `eps = sigma` down to `eps = sigma/32` -- topology is set by `sigma`, not
#' resolution -- while boundary length converges from 6016.8 (`eps = sigma`,
#' -11.5%) to 6799.1 (`eps = sigma/32`, reference). `eps = sigma/8` (6756.5,
#' -0.6% bias) is used internally as a resolution/cost tradeoff; there is
#' nothing left for the user to tune. `Boundary Length` is not scale-free --
#' compare `Boundary Length / sqrt(area)` across cores of different size.
#'
#' @return object of class `mif`, with:
#' \item{spatial}{each sample gains `density_compartment` (factor, levels
#'   `c(class1, class2)`) and `refined_density_compartment` (factor, levels
#'   `c(class1, "Interface", class2)`)}
#' \item{derived$density_boundary}{a named list (by `names(mif$spatial)`), one
#'   data frame per sample with columns `<sample_id>`, `piece`, `x`, `y` (0 rows
#'   if no boundary exists), carrying a `"call_info"` attribute recording the
#'   settings used, consumed internally by [plot_tissue_split()]}
#' \item{sample}{gains a `Boundary Length` column (`NA` for samples with no
#'   matching spatial frame; `0`, not `NA`, for a sample with data but no
#'   contour)}
#'
#' @export
#' @examples
#' library(dplyr)
#' x <- create_mif(clinical_data = example_clinical %>%
#'   mutate(deidentified_id = as.character(deidentified_id)),
#'   sample_data = example_summary %>%
#'   mutate(deidentified_id = as.character(deidentified_id)),
#'   spatial_list = example_spatial[1],
#'   patient_id = "deidentified_id",
#'   sample_id = "deidentified_sample")
#' x <- split_tissue(x, classifier = "Classifier.Label",
#'   class1 = "Tumor", class2 = "Stroma",
#'   sigma = 40, interface_width = 100)
split_tissue <- function(mif, classifier, class1, class2, sigma, interface_width,
                         workers = 1, overwrite = FALSE, xloc = NULL, yloc = NULL, ...) {
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))))) {
    stop("`split_tissue()` received unnamed arguments in `...`. ",
         "All arguments must be named.", call. = FALSE)
  }
  dep  <- intersect(names(dots), names(deprecated_arg_map("split_tissue")))
  apply_deprecated_args(dots[dep], "split_tissue")
  rest <- dots[setdiff(names(dots), dep)]
  unknown <- setdiff(names(rest), "filter_density")
  if (length(unknown)) {
    stop(sprintf("Unknown argument%s passed to `split_tissue()`: %s.\n",
                 if (length(unknown) > 1) "s" else "",
                 paste0("`", unknown, "`", collapse = ", ")),
         "The only argument accepted through `...` is `filter_density`.",
         call. = FALSE)
  }
  filter_density <- rest[["filter_density"]]

  if (!inherits(mif, "mif")) {
    stop("mIF should be of class `mif` created with function `create_mif()`\n",
         "\tTo check use `inherits(mif, 'mif')`")
  }
  if (length(mif$spatial) == 1 && identical(mif$spatial, list(NA))) {
    stop("`mif` has no spatial data (`create_mif(spatial_list = NULL)`); ",
         "`split_tissue()` needs per-cell spatial data.", call. = FALSE)
  }
  if (missing(sigma)) {
    stop("`sigma` has no default -- see `?split_tissue` for why. ",
         "Choose a bandwidth in the same units as your spatial coordinates.",
         call. = FALSE)
  }
  if (missing(interface_width)) {
    stop("`interface_width` has no default -- choose the full width of the ",
         "interface band, in the same units as `sigma`.", call. = FALSE)
  }
  if (!is.numeric(sigma) || length(sigma) != 1 || !is.finite(sigma) || sigma <= 0) {
    stop("`sigma` must be a single finite value greater than 0.", call. = FALSE)
  }
  if (!is.numeric(interface_width) || length(interface_width) != 1 ||
      !is.finite(interface_width) || interface_width < 0) {
    stop("`interface_width` must be a single finite value >= 0.", call. = FALSE)
  }
  for (nm in c("classifier", "class1", "class2")) {
    val <- get(nm)
    if (!is.character(val) || length(val) != 1) {
      stop("`", nm, "` must be a single character string.", call. = FALSE)
    }
  }
  if (identical(class1, class2)) {
    stop("`class1` and `class2` must differ (the density difference would be ",
         "identically zero).", call. = FALSE)
  }
  if ("Interface" %in% c(class1, class2)) {
    stop("`class1`/`class2` cannot be named \"Interface\" -- that level is ",
         "reserved for `refined_density_compartment`.", call. = FALSE)
  }
  if (!is.null(filter_density) && !is.function(filter_density)) {
    stop("`filter_density` must be a function.", call. = FALSE)
  }

  clashes <- character(0)
  if (any(vapply(mif$spatial, function(s)
      "density_compartment" %in% colnames(s) || "refined_density_compartment" %in% colnames(s),
      logical(1)))) {
    clashes <- c(clashes, "`density_compartment`/`refined_density_compartment` in `mif$spatial`")
  }
  if (!is.null(mif$derived$density_boundary)) {
    clashes <- c(clashes, "`mif$derived$density_boundary`")
  }
  if ("Boundary Length" %in% colnames(mif$sample)) {
    clashes <- c(clashes, "`Boundary Length` in `mif$sample`")
  }
  if (length(clashes) && !overwrite) {
    stop("`split_tissue()` has already been run on this mif -- found: ",
         paste(clashes, collapse = "; "), ".\n",
         "Set `overwrite = TRUE` to replace, or keep two mifs to compare two settings.",
         call. = FALSE)
  }

  eps <- density_pixel_size(sigma)

  # No RNG anywhere in this function (unlike ripleys_k()/pair_correlation()), so
  # there is deliberately no per-sample seed here.
  res <- parallel::mclapply(seq_along(mif$spatial), function(sample_i) {
    spat  <- add_cell_centres(mif$spatial[[sample_i]], xloc, yloc)
    label <- as.character(spat[[mif$sample_id]][1])
    cls   <- spat[[classifier]]
    if (is.null(cls)) {
      stop("Column \"", classifier, "\" not found in sample \"", label, "\".", call. = FALSE)
    }

    # Window from EVERY cell in the sample -- pinned by test-window-invariant.R.
    win <- spatstat.geom::convexhull.xy(spat$xloc, spat$yloc)
    pp  <- spatstat.geom::ppp(spat$xloc, spat$yloc, window = win, check = FALSE)
    check_density_budget(win, eps, sigma, label)

    keep1 <- class_mask(cls, class1)
    keep2 <- class_mask(cls, class2)
    if (!any(keep1) || !any(keep2)) {
      warning(sprintf(
        "Sample \"%s\" has %d \"%s\" cell%s and %d \"%s\" cell%s; the boundary ",
        label, sum(keep1), class1, if (sum(keep1) == 1) "" else "s",
        sum(keep2), class2, if (sum(keep2) == 1) "" else "s"),
        "will be empty and every cell will fall on one side (or NA).", call. = FALSE)
    }

    d  <- compartment_diff(pp, keep1, keep2, sigma, eps, filter_density)
    bd <- zero_contour(d$filtered, mif$sample_id, label)
    S  <- boundary_psp(bd, win)
    v  <- im_value_at(d$raw, spat$xloc, spat$yloc)

    list(label     = label,
        code      = ifelse(v > 0, 1L, ifelse(v < 0, 2L, NA_integer_)),
        interface = spatstat.geom::nncross(pp, S)$dist <= interface_width / 2,
        boundary  = bd,
        length    = boundary_length(S))
  }, mc.cores = workers, mc.preschedule = FALSE)

  failed <- vapply(res, inherits, logical(1), "try-error")
  if (any(failed)) {
    stop("split_tissue() failed on sample", if (sum(failed) > 1) "s" else "", " ",
         paste0("\"", names(mif$spatial)[failed], "\"", collapse = ", "), ":\n",
         paste(vapply(res[failed], as.character, character(1)), collapse = "\n"),
         call. = FALSE)
  }

  lev  <- c(class1, class2)
  levr <- c(class1, "Interface", class2)
  for (i in seq_along(mif$spatial)) {
    lab <- lev[res[[i]]$code]
    mif$spatial[[i]][["density_compartment"]] <- factor(lab, levels = lev)
    # Overwriting `lab` in place guarantees the two columns cannot disagree --
    # a cell that is Interface in `refined` always has the sign-derived label
    # in `density_compartment` it would have had without the interface rule.
    lab[res[[i]]$interface] <- "Interface"
    mif$spatial[[i]][["refined_density_compartment"]] <- factor(lab, levels = levr)
  }

  boundary <- lapply(res, function(r) r$boundary)
  names(boundary) <- names(mif$spatial)
  mismatched <- vapply(seq_along(mif$spatial), function(i)
    !identical(names(mif$spatial)[i], res[[i]]$label), logical(1))
  if (any(mismatched)) {
    warning("`names(mif$spatial)` disagrees with `mif$sample_id`'s value for sample",
           if (sum(mismatched) > 1) "s " else " ",
           paste0("\"", names(mif$spatial)[mismatched], "\"", collapse = ", "),
           " -- the `Boundary Length` join below will drop or misassign these.",
           call. = FALSE)
  }
  attr(boundary, "call_info") <- list(
    classifier = classifier, class1 = class1, class2 = class2,
    sigma = sigma, eps = eps, interface_width = interface_width,
    xloc = xloc, yloc = yloc, filter_density = filter_density,
    sample_id = mif$sample_id,
    spatialTIME_version = as.character(utils::packageVersion("spatialTIME")))
  mif$derived$density_boundary <- boundary

  lengths <- data.frame(sample_id_ = names(mif$spatial),
                        `Boundary Length` = vapply(res, function(r) r$length, numeric(1)),
                        check.names = FALSE, stringsAsFactors = FALSE)
  names(lengths)[1] <- mif$sample_id

  if (anyDuplicated(names(mif$spatial))) {
    stop("`mif$spatial` has duplicate sample names -- a left join into ",
         "`mif$sample` would fan out its rows. Fix the duplicate names first.",
         call. = FALSE)
  }
  if (!identical(class(mif$sample[[mif$sample_id]]), class(lengths[[mif$sample_id]]))) {
    stop("`mif$sample[[\"", mif$sample_id, "\"]]` is class \"",
         class(mif$sample[[mif$sample_id]])[1], "\" but the spatial sample names are class \"",
         class(lengths[[mif$sample_id]])[1], "\" -- these must match for the ",
         "`Boundary Length` join, and will not be silently coerced.", call. = FALSE)
  }
  unmatched <- setdiff(lengths[[mif$sample_id]], mif$sample[[mif$sample_id]])
  if (length(unmatched)) {
    warning("Sample", if (length(unmatched) > 1) "s " else " ",
           paste0("\"", unmatched, "\"", collapse = ", "),
           " in `mif$spatial` have no matching row in `mif$sample`; their ",
           "`Boundary Length` is dropped.", call. = FALSE)
  }

  mif$sample[["Boundary Length"]] <- NULL   # drop before the join, else left_join makes .x/.y
  mif$sample <- dplyr::left_join(mif$sample, lengths, by = mif$sample_id)

  mif
}
