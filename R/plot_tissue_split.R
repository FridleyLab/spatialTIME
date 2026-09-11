#' Plot density-based tissue compartments alongside the density difference
#'
#' @param mif object of class `mif` that has been through [split_tissue()]
#' @param which samples to plot: `NULL` for all, or a numeric/character vector
#'   indexing/naming elements of `mif$spatial`
#' @param compartment which compartment column to colour cells by: the
#'   3-level `"refined_density_compartment"` (default) or the 2-level
#'   `"density_compartment"`
#' @param panels which panel(s) to draw: `"both"` (default, one faceted plot
#'   with the density-difference raster and the compartment scatter
#'   side by side), `"density"` alone, or `"compartment"` alone
#' @param colors named character vector of colours for `class1`, `class2` and
#'   (for `compartment = "refined_density_compartment"`) `"Interface"`. If
#'   `NULL`, defaults to `RColorBrewer::brewer.pal(3, "Set1")[1:2]` for the two
#'   classes and `"grey20"` for `"Interface"`; `NA`/unclassified cells are
#'   always drawn `"grey70"`.
#' @param point_size size passed to [ggplot2::geom_point()] for the compartment panel
#' @param raster_max_pixels the density raster is recomputed at plot time (see
#'   Details) at up to this many pixels; above the cap resolution is capped
#'   for display only -- the drawn boundary line is unaffected and stays exact.
#' @param workers number of cores to use, one sample per core
#' @param filename,path if `filename` is given, every requested plot is also
#'   written to a single PDF at `file.path(path, filename)` (`path` defaults
#'   to the working directory; a trailing `.pdf` in `filename` is not doubled).
#' @param settings provenance normally read from
#'   `attr(mif$derived$density_boundary, "call_info")`, which [split_tissue()]
#'   attaches. That attribute does not survive `[` subsetting or
#'   `dplyr::bind_rows()` (e.g. after [merge_mifs()] disagrees across mifs), so
#'   pass the settings list explicitly to recover from that: a list with
#'   `classifier`, `class1`, `class2`, `sigma`, `eps`, `interface_width`,
#'   `xloc`, `yloc`, `filter_density` and `sample_id`.
#' @param ... accepts no arguments; present only so that a mistyped named
#'   argument above produces an informative "Unknown argument" error instead
#'   of being silently absorbed.
#'
#' @description
#' [split_tissue()] discards the density images and point patterns it computes
#' the boundary from, to avoid inflating the mif. `plot_tissue_split()`
#' recomputes them on demand, at plot time, from the settings recorded by
#' [split_tissue()], and draws the *stored* boundary polyline (never a
#' recontoured one) on top of both a raster of the density difference and a
#' scatter of the per-cell compartment label.
#'
#' @return a named list of `ggplot` objects, one per requested sample, named
#'   like `mif$spatial` -- **not** the `mif`. Unlike [plot_immunoflo()], this
#'   deliberately does not attach the plots to `mif$derived`: each plot
#'   captures a raster data frame that can be as large as the sample itself,
#'   and `merge_mifs()`'s handling of list-valued derived slots is fragile
#'   enough already (see `NEWS.md`) without doubling the exposure.
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
#' plots <- plot_tissue_split(x)
plot_tissue_split <- function(mif, which = NULL,
                              compartment = c("refined_density_compartment",
                                              "density_compartment"),
                              panels = c("both", "density", "compartment"),
                              colors = NULL, point_size = 0.4,
                              raster_max_pixels = 5e5, workers = 1,
                              filename = NULL, path = NULL, settings = NULL, ...) {
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))))) {
    stop("`plot_tissue_split()` received unnamed arguments in `...`. ",
         "All arguments must be named.", call. = FALSE)
  }
  dep <- intersect(names(dots), names(deprecated_arg_map("plot_tissue_split")))
  apply_deprecated_args(dots[dep], "plot_tissue_split")
  unknown <- setdiff(names(dots), dep)
  if (length(unknown)) {
    stop(sprintf("Unknown argument%s passed to `plot_tissue_split()`: %s.\n",
                 if (length(unknown) > 1) "s" else "",
                 paste0("`", unknown, "`", collapse = ", ")),
         "`plot_tissue_split()` accepts no arguments through `...`.",
         call. = FALSE)
  }

  if (!inherits(mif, "mif")) {
    stop("mIF should be of class `mif` created with function `create_mif()`\n",
         "\tTo check use `inherits(mif, 'mif')`")
  }
  compartment <- match.arg(compartment)
  panels <- match.arg(panels)

  cfg <- split_tissue_settings(mif, settings)
  boundary_list <- mif$derived$density_boundary
  if (is.null(boundary_list)) {
    stop("`mif$derived$density_boundary` not found -- run `split_tissue()` first.",
         call. = FALSE)
  }

  idx <- if (is.null(which)) seq_along(mif$spatial) else which
  if (is.character(idx)) {
    resolved <- match(idx, names(mif$spatial))
    if (anyNA(resolved)) {
      stop("`which` names not found in `names(mif$spatial)`: ",
           paste0("\"", idx[is.na(resolved)], "\"", collapse = ", "), ".", call. = FALSE)
    }
    idx <- resolved
  }

  levs <- if (compartment == "refined_density_compartment") {
    c(cfg$class1, "Interface", cfg$class2)
  } else {
    c(cfg$class1, cfg$class2)
  }
  if (is.null(colors)) {
    base_cols <- RColorBrewer::brewer.pal(3, "Set1")
    colors <- stats::setNames(c(base_cols[1], "grey20", base_cols[2]),
                              c(cfg$class1, "Interface", cfg$class2))
  }
  colors <- colors[levs]

  plots <- parallel::mclapply(idx, function(i) {
    label <- names(mif$spatial)[i]
    spat  <- add_cell_centres(mif$spatial[[i]], cfg$xloc, cfg$yloc)
    cls   <- spat[[cfg$classifier]]

    win <- spatstat.geom::convexhull.xy(spat$xloc, spat$yloc)
    pp  <- spatstat.geom::ppp(spat$xloc, spat$yloc, window = win, check = FALSE)
    keep1 <- class_mask(cls, cfg$class1)
    keep2 <- class_mask(cls, cfg$class2)

    a <- spatstat.geom::area(spatstat.geom::Frame(win))
    eps_plot <- max(cfg$eps, sqrt(a / raster_max_pixels))
    d <- compartment_diff(pp, keep1, keep2, cfg$sigma, eps_plot, cfg$filter_density)
    diff_im <- d$filtered

    ras <- as.data.frame(diff_im)
    names(ras) <- c("x", "y", "value")
    ras <- ras[!is.na(ras$value), ]
    ras$panel <- "Density difference"
    m <- if (nrow(ras)) max(abs(ras$value)) else 1

    pts <- data.frame(x = spat$xloc, y = spat$yloc,
                      compartment = factor(as.character(spat[[compartment]]), levels = levs),
                      panel = "Compartment")

    bd <- boundary_list[[i]][, c("piece", "x", "y")]
    if (nrow(bd)) {
      bd2 <- rbind(cbind(bd, panel = "Density difference"),
                  cbind(bd, panel = "Compartment"))
    } else {
      bd2 <- data.frame(piece = integer(0), x = numeric(0), y = numeric(0),
                        panel = character(0))
    }

    show_density <- panels %in% c("both", "density")
    show_compartment <- panels %in% c("both", "compartment")
    ras <- if (show_density) ras else ras[0, ]
    pts <- if (show_compartment) pts else pts[0, ]
    bd2 <- bd2[bd2$panel %in% c(
      if (show_density) "Density difference",
      if (show_compartment) "Compartment"), , drop = FALSE]

    p <- ggplot2::ggplot()
    if (show_density) {
      p <- p +
        ggplot2::geom_raster(data = ras,
          ggplot2::aes(.data$x, .data$y, fill = .data$value)) +
        ggplot2::scale_fill_gradient2(
          name = sprintf("%s - %s\ndensity", cfg$class1, cfg$class2),
          low = "#2166AC", mid = "white", high = "#B2182B",
          midpoint = 0, limits = c(-m, m), oob = scales::squish)
    }
    if (show_compartment) {
      p <- p +
        ggplot2::geom_point(data = pts,
          ggplot2::aes(.data$x, .data$y, colour = .data$compartment),
          size = point_size) +
        ggplot2::scale_colour_manual(NULL, values = colors, drop = FALSE,
                                     na.value = "grey70")
    }
    p <- p +
      ggplot2::geom_path(data = bd2,
        ggplot2::aes(.data$x, .data$y, group = .data$piece),
        colour = "black", linewidth = 0.3) +
      ggplot2::coord_equal() +
      ggplot2::scale_y_reverse(breaks = scales::pretty_breaks(5)) +
      ggplot2::ggtitle(paste0("ID: ", label)) +
      ggplot2::theme_bw(base_size = 18) +
      ggplot2::theme(axis.title = ggplot2::element_blank())
    if (panels == "both") {
      p <- p + ggplot2::facet_wrap(~ panel)
    }
    p
  }, mc.cores = workers)
  names(plots) <- names(mif$spatial)[idx]

  failed <- vapply(plots, inherits, logical(1), "try-error")
  if (any(failed)) {
    stop("plot_tissue_split() failed on sample", if (sum(failed) > 1) "s" else "", " ",
         paste0("\"", names(plots)[failed], "\"", collapse = ", "), ":\n",
         paste(vapply(plots[failed], as.character, character(1)), collapse = "\n"),
         call. = FALSE)
  }

  if (!is.null(filename)) {
    out_dir  <- if (is.null(path)) "." else path
    out_stem <- sub("\\.pdf$", "", filename, ignore.case = TRUE)
    grDevices::pdf(file.path(out_dir, paste0(out_stem, ".pdf")), height = 8, width = 14)
    on.exit(grDevices::dev.off(), add = TRUE)
    invisible(lapply(plots, print))
  }

  plots
}
