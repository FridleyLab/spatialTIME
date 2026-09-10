skip_if_not_installed("ggplot2")

# Assertions go through ggplot_build()/layer_data() rather than snapshots, so
# nothing has to be rendered and the tests do not break when ggplot2 changes its
# default theme.

# Cores 1 and 2 are the SPARSEST (17 and 21 CD3+ cells out of ~3800), so thinning
# them leaves about one positive cell and the marker layer comes out empty --
# which looks like a plotting bug rather than a fixture problem. Use the two dense
# TMA3 cores.
plot_mif <- function() example_mif(which = c("TMA3_[9,K].tif", "TMA3_[8,U].tif"),
                                   n_cells = 400)
plot_markers <- function() mnames_good()[1:3]

test_that("plot_immunoflo returns the mif with one ggplot per sample", {
  out <- suppressWarnings(suppressMessages(
    plot_immunoflo(plot_mif(), plot_title = "deidentified_sample",
                   mnames = plot_markers())))
  expect_s3_class(out, "mif")
  expect_named(out$derived$spatial_plots, names(plot_mif()$spatial))
  for (p in out$derived$spatial_plots) expect_s3_class(p, "ggplot")
})

test_that("cell_type maps the column, giving one shape per level", {
  # THE regression: aes(shape = cell_type) mapped the STRING "Classifier.Label",
  # not the column, so one shape was drawn for every cell and the legend had a
  # single entry named after the column. The documented example passes
  # cell_type = "Classifier.Label" and rendered Tumor and Stroma identically.
  mif <- plot_mif()
  out <- suppressWarnings(suppressMessages(
    plot_immunoflo(mif, plot_title = "deidentified_sample", mnames = plot_markers(),
                   cell_type = "Classifier.Label")))
  p <- out$derived$spatial_plots[[1]]
  n_levels <- length(unique(mif$spatial[[1]]$Classifier.Label))
  expect_gt(n_levels, 1)   # guard the guard

  marker_layer <- suppressWarnings(ggplot2::ggplot_build(p))$data[[2]]
  expect_gt(nrow(marker_layer), 0)          # guard: an empty layer proves nothing
  expect_length(unique(marker_layer$shape), n_levels)
  # And the mapping is to the column, not the constant string "cell_type".
  expect_false(identical(deparse(rlang::quo_get_expr(p$mapping$shape)), "cell_type"))
})

test_that("a classifier with more than two levels still builds", {
  # scale_shape_manual(values = c(3, 16)) hard-coded exactly two shapes, so a third
  # level errored -- at PRINT time, so plot_immunoflo() itself appeared to succeed.
  mif <- plot_mif()
  mif$spatial[[1]]$Classifier.Label <- rep(c("A", "B", "C"),
                                          length.out = nrow(mif$spatial[[1]]))
  mif$spatial[[2]]$Classifier.Label <- rep(c("A", "B", "C"),
                                          length.out = nrow(mif$spatial[[2]]))
  out <- suppressWarnings(suppressMessages(
    plot_immunoflo(mif, plot_title = "deidentified_sample", mnames = plot_markers(),
                   cell_type = "Classifier.Label")))
  built <- suppressWarnings(ggplot2::ggplot_build(out$derived$spatial_plots[[1]]))
  expect_gt(nrow(built$data[[2]]), 0)
  expect_length(unique(built$data[[2]]$shape), 3)
})

test_that("no device is left open or closed behind the caller's back", {
  # on.exit(dev.off()) PLUS an explicit grDevices::dev.off() meant a second
  # dev.off() fired on exit: with no other device open the function errored after
  # writing the pdf, and with a user device open it silently closed the caller's.
  tmp <- withr::local_tempdir()
  before_n <- length(grDevices::dev.list())

  expect_no_error(suppressWarnings(suppressMessages(
    plot_immunoflo(plot_mif(), plot_title = "deidentified_sample",
                   mnames = plot_markers(), filename = "plots", path = tmp))))

  expect_equal(length(grDevices::dev.list()), before_n)
  expect_true(file.exists(file.path(tmp, "plots.pdf")))
})

test_that("path is honoured and .pdf is not doubled", {
  # `path` was documented but never referenced, so output always went to getwd().
  tmp <- withr::local_tempdir()
  suppressWarnings(suppressMessages(
    plot_immunoflo(plot_mif(), plot_title = "deidentified_sample",
                   mnames = plot_markers(), filename = "already.pdf", path = tmp)))
  expect_true(file.exists(file.path(tmp, "already.pdf")))
  expect_false(file.exists(file.path(tmp, "already.pdf.pdf")))
})

test_that("no scale-replacement message is emitted", {
  # scale_y_continuous() was immediately replaced by scale_y_reverse(), which both
  # discarded the pretty-breaks intent and warned once per sample.
  msgs <- capture.output(
    suppressWarnings(invisible(
      plot_immunoflo(plot_mif(), plot_title = "deidentified_sample",
                     mnames = plot_markers(), cell_type = "Classifier.Label"))),
    type = "message")
  expect_length(grep("already present", msgs), 0)
})

test_that("the y axis is reversed, as images require", {
  out <- suppressWarnings(suppressMessages(
    plot_immunoflo(plot_mif(), plot_title = "deidentified_sample",
                   mnames = plot_markers())))
  p <- out$derived$spatial_plots[[1]]
  # Checked behaviourally rather than by introspecting the scale object, whose
  # internals differ across ggplot2 versions: under scale_y_reverse the built y
  # values are the negation of the source coordinates.
  built <- suppressWarnings(ggplot2::ggplot_build(p))$data[[2]]
  expect_gt(nrow(built), 0)
  expect_true(all(built$y <= 0))
  # And a y scale is present exactly once, so nothing is silently replacing it.
  y_scales <- vapply(p$scales$scales, function(s) "y" %in% s$aesthetics, logical(1))
  expect_equal(sum(y_scales), 1)
})

test_that("plot_immunoflo rejects a non-mif and a missing mif", {
  expect_error(plot_immunoflo(mnames = "A", plot_title = "x"), "MIF is missing")
  expect_error(plot_immunoflo(list(), plot_title = "x", mnames = "A"),
               "Please use a mif object")
})
