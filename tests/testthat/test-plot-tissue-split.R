skip_if_not_installed("ggplot2")

# plot_tissue_split() recomputes the density difference at plot time and draws
# the STORED boundary (never a recontoured one) over it. Assertions go through
# ggplot_build()/layer_data() and the scale objects rather than snapshots, per
# test-plot-immunoflo.R's convention.

split_fixture <- function(interface_width = 100) {
  split_tissue(halfplane_mif(), classifier = "Classifier.Label",
              class1 = "Tumor", class2 = "Stroma",
              sigma = 40, interface_width = interface_width)
}

fill_scale <- function(p) {
  p$scales$scales[[which(vapply(p$scales$scales,
                                function(s) "fill" %in% s$aesthetics, logical(1)))]]
}

y_scales <- function(p) {
  vapply(p$scales$scales, function(s) "y" %in% s$aesthetics, logical(1))
}

test_that("plot_tissue_split returns a named list of ggplots, not the mif", {
  mif <- split_fixture()
  plots <- plot_tissue_split(mif)
  expect_type(plots, "list")
  expect_false(inherits(plots, "mif"))
  expect_named(plots, names(mif$spatial))
  for (p in plots) expect_s3_class(p, "ggplot")
  expect_identical(mif$derived, mif$derived)   # unchanged reference, sanity
})

test_that("which filters to the requested sample(s)", {
  mif <- split_fixture()
  plots <- plot_tissue_split(mif, which = 1)
  expect_length(plots, 1)
  expect_identical(names(plots), names(mif$spatial)[1])

  by_name <- plot_tissue_split(mif, which = names(mif$spatial)[1])
  expect_identical(names(by_name), names(mif$spatial)[1])

  expect_error(plot_tissue_split(mif, which = "nope"), "not found")
})

test_that("layer row counts match panels = \"both\"", {
  mif <- split_fixture()
  p <- plot_tissue_split(mif)[[1]]
  built <- expect_no_warning(ggplot2::ggplot_build(p))
  raster_rows <- nrow(built$data[[1]])
  point_rows  <- nrow(built$data[[2]])
  path_rows   <- nrow(built$data[[3]])
  expect_gt(raster_rows, 0)
  expect_equal(point_rows, nrow(mif$spatial[[1]]))
  expect_equal(path_rows, 2 * nrow(mif$derived$density_boundary[[1]]))
})

test_that("panels = \"density\" drops the point layer, \"compartment\" drops the raster", {
  mif <- split_fixture()
  bd_n <- nrow(mif$derived$density_boundary[[1]])

  dens <- plot_tissue_split(mif, panels = "density")[[1]]
  built_d <- expect_no_warning(ggplot2::ggplot_build(dens))
  expect_length(built_d$data, 2)   # raster, path
  expect_gt(nrow(built_d$data[[1]]), 0)
  expect_equal(nrow(built_d$data[[2]]), bd_n)

  comp <- plot_tissue_split(mif, panels = "compartment")[[1]]
  built_c <- expect_no_warning(ggplot2::ggplot_build(comp))
  expect_length(built_c$data, 2)   # point, path
  expect_equal(nrow(built_c$data[[1]]), nrow(mif$spatial[[1]]))
  expect_equal(nrow(built_c$data[[2]]), bd_n)
})

test_that("the drawn path is the STORED boundary, not a recomputed one", {
  mif <- split_fixture()
  moved <- mif
  moved$derived$density_boundary[[1]]$x <- moved$derived$density_boundary[[1]]$x + 50
  attr(moved$derived$density_boundary, "call_info") <- attr(mif$derived$density_boundary, "call_info")

  p1 <- plot_tissue_split(mif)[[1]]
  p2 <- plot_tissue_split(moved)[[1]]
  x1 <- ggplot2::ggplot_build(p1)$data[[3]]$x
  x2 <- ggplot2::ggplot_build(p2)$data[[3]]$x
  expect_equal(sort(x2), sort(x1) + 50, tolerance = 1e-9)
})

test_that("the fill scale is symmetric about zero", {
  mif <- split_fixture()
  p <- plot_tissue_split(mif)[[1]]
  lim <- fill_scale(p)$get_limits()
  expect_equal(lim[1], -lim[2], tolerance = 1e-9)
})

test_that("compartment = density_compartment gives 2 colour levels instead of 3", {
  mif <- split_fixture()
  refined <- plot_tissue_split(mif, compartment = "refined_density_compartment")[[1]]
  simple  <- plot_tissue_split(mif, compartment = "density_compartment")[[1]]
  colour_scale <- function(p) {
    built <- ggplot2::ggplot_build(p)$plot
    built$scales$scales[[which(vapply(built$scales$scales,
      function(s) "colour" %in% s$aesthetics, logical(1)))]]
  }
  expect_length(colour_scale(refined)$get_limits(), 3)
  expect_length(colour_scale(simple)$get_limits(), 2)
})

test_that("the y axis is reversed and present exactly once", {
  mif <- split_fixture()
  p <- plot_tissue_split(mif)[[1]]
  built <- ggplot2::ggplot_build(p)
  expect_true(all(built$data[[2]]$y <= 0))
  expect_equal(sum(y_scales(p)), 1)
})

test_that("a mif without split_tissue() errors pointing at split_tissue", {
  mif <- halfplane_mif()
  expect_error(plot_tissue_split(mif), "split_tissue")
})

test_that("stripped call_info errors mentioning settings, and settings= recovers", {
  mif <- split_fixture()
  ci <- attr(mif$derived$density_boundary, "call_info")
  stripped <- mif
  attr(stripped$derived$density_boundary, "call_info") <- NULL

  expect_error(plot_tissue_split(stripped), "settings")
  recovered <- plot_tissue_split(stripped, settings = ci)
  expect_s3_class(recovered[[1]], "ggplot")
})

test_that("partial settings errors naming the missing field", {
  mif <- split_fixture()
  ci <- attr(mif$derived$density_boundary, "call_info")
  partial <- ci[setdiff(names(ci), "sigma")]
  expect_error(plot_tissue_split(mif, settings = partial), "sigma")
})

test_that("filename/path writes one PDF, .pdf is not doubled, no device left open", {
  mif <- split_fixture()
  tmp <- withr::local_tempdir()
  before_n <- length(grDevices::dev.list())

  expect_no_error(plot_tissue_split(mif, filename = "tissue_split.pdf", path = tmp))
  expect_equal(length(grDevices::dev.list()), before_n)
  expect_true(file.exists(file.path(tmp, "tissue_split.pdf")))
  expect_false(file.exists(file.path(tmp, "tissue_split.pdf.pdf")))
})

test_that("every plot builds without warning", {
  mif <- split_fixture()
  for (p in plot_tissue_split(mif)) {
    expect_no_warning(ggplot2::ggplot_build(p))
  }
})

test_that("raster_max_pixels caps raster resolution but not the boundary path", {
  mif <- split_fixture()
  full   <- plot_tissue_split(mif)[[1]]
  capped <- plot_tissue_split(mif, raster_max_pixels = 100)[[1]]
  built_full   <- ggplot2::ggplot_build(full)
  built_capped <- ggplot2::ggplot_build(capped)
  expect_lt(nrow(built_capped$data[[1]]), nrow(built_full$data[[1]]))
  expect_identical(built_capped$data[[3]]$x, built_full$data[[3]]$x)
})

test_that("plot_tissue_split rejects a non-mif", {
  expect_error(plot_tissue_split(list()), "class `mif`")
})

test_that("an unrecognised argument in ... is an error naming it", {
  mif <- split_fixture()
  expect_error(plot_tissue_split(mif, bogus = TRUE), "Unknown argument.*bogus")
})
