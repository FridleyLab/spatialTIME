# split_tissue() computes a KDE difference, extracts its exact zero level set with
# contourLines(), and labels every cell by which side of that boundary it falls on.
# See plans/create-an-implementation-plan-squishy-turing.md for the numbers below --
# they are reproduced from real runs, not invented tolerances.

full_mif <- function() {
  create_mif(
    clinical_data = example_clinical %>%
      dplyr::mutate(deidentified_id = as.character(deidentified_id)),
    sample_data = example_summary %>%
      dplyr::mutate(deidentified_id = as.character(deidentified_id)),
    spatial_list = example_spatial,
    patient_id = "deidentified_id", sample_id = "deidentified_sample")
}

# name -> (cells, contour pieces, boundary length), measured at sigma = 40,
# interface_width = 100 on the shipped example data.
full_targets <- function() {
  list(
    `TMA1_[3,B].tif` = list(cells = 3803, pieces = 10, length = 3238.2),
    `TMA2_[3,B].tif` = list(cells = 3008, pieces = 4,  length = 5189.2),
    `TMA3_[7,B].tif` = list(cells = 1850, pieces = 6,  length = 2346.8),
    `TMA3_[9,K].tif` = list(cells = 1803, pieces = 9,  length = 6756.5),
    `TMA3_[8,U].tif` = list(cells = 2318, pieces = 8,  length = 4469.4)
  )
}

# ---- Structure --------------------------------------------------------------

test_that("split_tissue returns a mif with both columns on every frame", {
  out <- split_tissue(halfplane_mif(), classifier = "Classifier.Label",
                      class1 = "Tumor", class2 = "Stroma",
                      sigma = 40, interface_width = 100)
  expect_s3_class(out, "mif")
  for (nm in names(out$spatial)) {
    expect_true(all(c("density_compartment", "refined_density_compartment") %in%
                    colnames(out$spatial[[nm]])), info = nm)
  }
})

test_that("factor levels are exact and in order, and nothing else about spatial changes", {
  before <- halfplane_mif()
  out <- split_tissue(before, classifier = "Classifier.Label",
                      class1 = "Tumor", class2 = "Stroma",
                      sigma = 40, interface_width = 100)
  s <- out$spatial[[1]]
  expect_identical(levels(s$density_compartment), c("Tumor", "Stroma"))
  expect_identical(levels(s$refined_density_compartment), c("Tumor", "Interface", "Stroma"))
  expect_equal(nrow(s), nrow(before$spatial[[1]]))
  old_cols <- setdiff(colnames(before$spatial[[1]]), c("xloc", "yloc"))
  expect_identical(s[old_cols], before$spatial[[1]][old_cols])
})

test_that("derived$density_boundary is a named list of data frames, not a data frame", {
  out <- split_tissue(halfplane_mif(), classifier = "Classifier.Label",
                      class1 = "Tumor", class2 = "Stroma",
                      sigma = 40, interface_width = 100)
  b <- out$derived$density_boundary
  expect_type(b, "list")
  expect_false(is.data.frame(b))
  expect_identical(names(b), names(out$spatial))
  el <- b[[1]]
  expect_identical(colnames(el), c("deidentified_sample", "piece", "x", "y"))
  expect_type(el$deidentified_sample, "character")
  expect_type(el$piece, "integer")
  expect_type(el$x, "double")
  expect_type(el$y, "double")
})

# ---- halfplane_mif(): exact geometry ----------------------------------------

test_that("halfplane_mif recovers the exact half-plane boundary", {
  mif <- halfplane_mif()
  out <- split_tissue(mif, classifier = "Classifier.Label",
                      class1 = "Tumor", class2 = "Stroma",
                      sigma = 40, interface_width = 100)
  bd <- out$derived$density_boundary[[1]]
  eps <- 40 / 8
  expect_equal(length(unique(bd$piece)), 1)
  expect_equal(length(unique(bd$x)), 1)
  expect_equal(unique(bd$x), 500 + eps / 2, tolerance = 1e-6)
  height <- diff(range(mif$spatial[[1]]$YMin))
  expect_equal(out$sample$`Boundary Length`[1], height - eps, tolerance = 0.01)

  s <- out$spatial[[1]]
  expect_equal(sum(is.na(s$density_compartment)), 0)
  truth <- ifelse(s$XMin < 500, "Tumor", "Stroma")
  expect_true(all(as.character(s$density_compartment) == truth))
})

test_that("interface band on halfplane_mif matches interface_width", {
  mif <- halfplane_mif()
  out <- split_tissue(mif, classifier = "Classifier.Label",
                      class1 = "Tumor", class2 = "Stroma",
                      sigma = 40, interface_width = 100)
  s <- out$spatial[[1]]
  step <- 20   # lattice spacing in halfplane_mif()
  is_interface <- s$refined_density_compartment == "Interface"
  expect_true(all(abs(s$XMin[is_interface] - 500) <= 100 / 2 + step))
  expect_true(all(abs(s$XMin[!is_interface] - 500) >= 100 / 2 - step))
})

test_that("one_class_mif takes the no-boundary path with no Interface cells", {
  expect_warning(
    out <- split_tissue(one_class_mif(), classifier = "Classifier.Label",
                        class1 = "Tumor", class2 = "Stroma",
                        sigma = 40, interface_width = 100),
    "will be empty"
  )
  expect_equal(nrow(out$derived$density_boundary[[1]]), 0)
  expect_equal(out$sample$`Boundary Length`[1], 0)
  expect_false(any(out$spatial[[1]]$refined_density_compartment == "Interface"))
})

# ---- Real data ---------------------------------------------------------------

test_that("reported Boundary Length reproduces the target table", {
  out <- split_tissue(full_mif(), classifier = "Classifier.Label",
                      class1 = "Tumor", class2 = "Stroma",
                      sigma = 40, interface_width = 100, overwrite = TRUE)
  targets <- full_targets()
  for (nm in names(targets)) {
    len <- out$sample$`Boundary Length`[out$sample$deidentified_sample == nm]
    pieces <- length(unique(out$derived$density_boundary[[nm]]$piece))
    expect_equal(len, targets[[nm]]$length, tolerance = 0.05, info = nm)
    expect_equal(pieces, targets[[nm]]$pieces, info = nm)
    expect_equal(nrow(out$spatial[[nm]]), targets[[nm]]$cells, info = nm)
  }
})

test_that("per-piece polyline length sums to the reported Boundary Length", {
  out <- split_tissue(full_mif(), classifier = "Classifier.Label",
                      class1 = "Tumor", class2 = "Stroma",
                      sigma = 40, interface_width = 100, overwrite = TRUE)
  for (nm in names(out$derived$density_boundary)) {
    b <- out$derived$density_boundary[[nm]]
    naive <- if (!nrow(b)) 0 else sum(vapply(split(b, b$piece), function(p)
      sum(sqrt(diff(p$x)^2 + diff(p$y)^2)), numeric(1)))
    reported <- out$sample$`Boundary Length`[out$sample$deidentified_sample == nm]
    expect_equal(naive, reported, tolerance = 1e-9, info = nm)
  }
})

test_that("the interface flag reproduces independently from the stored boundary", {
  mif <- full_mif()
  out <- split_tissue(mif, classifier = "Classifier.Label",
                      class1 = "Tumor", class2 = "Stroma",
                      sigma = 40, interface_width = 100, overwrite = TRUE)
  nm <- "TMA3_[9,K].tif"
  s <- add_cell_centres(out$spatial[[nm]], NULL, NULL)
  win <- spatstat.geom::convexhull.xy(s$xloc, s$yloc)
  pp  <- spatstat.geom::ppp(s$xloc, s$yloc, window = win, check = FALSE)
  S   <- boundary_psp(out$derived$density_boundary[[nm]], win)
  dist <- spatstat.geom::nncross(pp, S)$dist
  expect_identical(dist <= 100 / 2, s$refined_density_compartment == "Interface")
})

test_that("refined_density_compartment agrees with density_compartment outside Interface", {
  out <- split_tissue(full_mif(), classifier = "Classifier.Label",
                      class1 = "Tumor", class2 = "Stroma",
                      sigma = 40, interface_width = 100, overwrite = TRUE)
  for (nm in names(out$spatial)) {
    s <- out$spatial[[nm]]
    not_interface <- s$refined_density_compartment != "Interface"
    expect_identical(as.character(s$density_compartment)[not_interface],
                     as.character(s$refined_density_compartment)[not_interface], info = nm)
  }
})

test_that("Interface count is monotone non-decreasing in interface_width", {
  mif <- full_mif()
  widths <- c(0, 25, 50, 100, 200)
  counts <- vapply(widths, function(w) {
    out <- split_tissue(mif, classifier = "Classifier.Label",
                        class1 = "Tumor", class2 = "Stroma",
                        sigma = 40, interface_width = w, overwrite = TRUE)
    sum(out$spatial[["TMA3_[9,K].tif"]]$refined_density_compartment == "Interface")
  }, integer(1))
  expect_true(all(diff(counts) >= 0))
})

test_that("swapping class1/class2 leaves geometry and per-cell labels unchanged", {
  # The sign of the field AND the level order both flip on a swap, so the label
  # actually assigned to a given cell is invariant to argument order -- which
  # class the caller wrote first must not change the answer.
  mif <- full_mif()
  a <- split_tissue(mif, classifier = "Classifier.Label",
                    class1 = "Tumor", class2 = "Stroma",
                    sigma = 40, interface_width = 100, overwrite = TRUE)
  b <- split_tissue(mif, classifier = "Classifier.Label",
                    class1 = "Stroma", class2 = "Tumor",
                    sigma = 40, interface_width = 100, overwrite = TRUE)
  expect_equal(a$sample$`Boundary Length`, b$sample$`Boundary Length`, tolerance = 1e-9)
  for (nm in names(a$spatial)) {
    sa <- a$spatial[[nm]]$density_compartment
    sb <- b$spatial[[nm]]$density_compartment
    expect_identical(as.character(sa), as.character(sb), info = nm)
  }
})

test_that("no cell is left NA on real data (lookup-NAs are resolved)", {
  out <- split_tissue(full_mif(), classifier = "Classifier.Label",
                      class1 = "Tumor", class2 = "Stroma",
                      sigma = 40, interface_width = 100, overwrite = TRUE)
  for (nm in names(out$spatial)) {
    expect_equal(sum(is.na(out$spatial[[nm]]$density_compartment)), 0, info = nm)
  }
})

# ---- Boundary Length join -----------------------------------------------------

test_that("Boundary Length has exactly 5 non-NA values of 229 rows", {
  out <- split_tissue(full_mif(), classifier = "Classifier.Label",
                      class1 = "Tumor", class2 = "Stroma",
                      sigma = 40, interface_width = 100, overwrite = TRUE)
  expect_equal(nrow(out$sample), 229)
  expect_equal(sum(!is.na(out$sample$`Boundary Length`)), 5)
})

test_that("a sample with no contour gets 0, not NA, in Boundary Length", {
  out <- split_tissue(one_class_mif(), classifier = "Classifier.Label",
                      class1 = "Tumor", class2 = "Stroma",
                      sigma = 40, interface_width = 100) |> suppressWarnings()
  expect_identical(out$sample$`Boundary Length`, 0)
})

test_that("repeated overwrite = TRUE never leaves .x/.y suffixed columns", {
  mif <- halfplane_mif()
  out1 <- split_tissue(mif, classifier = "Classifier.Label",
                       class1 = "Tumor", class2 = "Stroma",
                       sigma = 40, interface_width = 100)
  out2 <- split_tissue(out1, classifier = "Classifier.Label",
                       class1 = "Tumor", class2 = "Stroma",
                       sigma = 40, interface_width = 100, overwrite = TRUE)
  expect_false(any(grepl("Boundary Length\\.[xy]$", colnames(out2$sample))))
  expect_true("Boundary Length" %in% colnames(out2$sample))
})

test_that("mismatched sample-id classes between spatial names and mif$sample error", {
  mif <- halfplane_mif()
  mif$sample$deidentified_sample <- 1L
  expect_error(
    split_tissue(mif, classifier = "Classifier.Label", class1 = "Tumor", class2 = "Stroma",
                sigma = 40, interface_width = 100),
    "class"
  )
})

test_that("duplicate spatial sample names error rather than fanning out mif$sample", {
  mif <- halfplane_mif()
  mif$spatial[[2]] <- mif$spatial[[1]]
  names(mif$spatial) <- rep(names(mif$spatial)[1], 2)
  expect_error(
    split_tissue(mif, classifier = "Classifier.Label", class1 = "Tumor", class2 = "Stroma",
                sigma = 40, interface_width = 100),
    "duplicate sample names"
  )
})

# ---- Provenance ---------------------------------------------------------------

test_that("call_info round-trips every setting, including xloc/yloc and the closure", {
  f <- function(im) im
  out <- split_tissue(halfplane_mif(), classifier = "Classifier.Label",
                      class1 = "Tumor", class2 = "Stroma",
                      sigma = 40, interface_width = 100,
                      xloc = "XMin", yloc = "YMin", filter_density = f)
  ci <- attr(out$derived$density_boundary, "call_info")
  expect_identical(ci$classifier, "Classifier.Label")
  expect_identical(ci$class1, "Tumor")
  expect_identical(ci$class2, "Stroma")
  expect_identical(ci$sigma, 40)
  expect_identical(ci$eps, 40 / 8)
  expect_identical(ci$interface_width, 100)
  expect_identical(ci$xloc, "XMin")
  expect_identical(ci$yloc, "YMin")
  expect_identical(ci$sample_id, "deidentified_sample")
  expect_true(is.function(ci$filter_density))
  expect_identical(ci$filter_density, f)
  expect_identical(ci$spatialTIME_version, as.character(utils::packageVersion("spatialTIME")))
})

test_that("split_tissue_settings validates a partial settings list", {
  out <- split_tissue(halfplane_mif(), classifier = "Classifier.Label",
                      class1 = "Tumor", class2 = "Stroma",
                      sigma = 40, interface_width = 100)
  full <- attr(out$derived$density_boundary, "call_info")
  expect_identical(split_tissue_settings(out, NULL), full[names(split_tissue_settings(out, NULL))])
  partial <- full[setdiff(names(full), "sigma")]
  expect_error(split_tissue_settings(out, partial), "missing required field")
})

# ---- overwrite ------------------------------------------------------------------

test_that("overwrite = FALSE on a populated mif errors naming all clashing locations", {
  out <- split_tissue(halfplane_mif(), classifier = "Classifier.Label",
                      class1 = "Tumor", class2 = "Stroma",
                      sigma = 40, interface_width = 100)
  err <- tryCatch(
    split_tissue(out, classifier = "Classifier.Label", class1 = "Tumor", class2 = "Stroma",
                sigma = 40, interface_width = 100),
    error = function(e) conditionMessage(e))
  expect_match(err, "density_compartment")
  expect_match(err, "density_boundary")
  expect_match(err, "Boundary Length")
})

test_that("overwrite = TRUE replaces cleanly, same nrow and column count", {
  mif <- halfplane_mif()
  out1 <- split_tissue(mif, classifier = "Classifier.Label",
                       class1 = "Tumor", class2 = "Stroma",
                       sigma = 40, interface_width = 100)
  out2 <- split_tissue(out1, classifier = "Classifier.Label",
                       class1 = "Tumor", class2 = "Stroma",
                       sigma = 40, interface_width = 100, overwrite = TRUE)
  expect_equal(nrow(out2$spatial[[1]]), nrow(out1$spatial[[1]]))
  expect_equal(ncol(out2$spatial[[1]]), ncol(out1$spatial[[1]]))
  expect_equal(ncol(out2$sample), ncol(out1$sample))
})

test_that("clash detection fires even with only a partial hand-built state", {
  mif <- halfplane_mif()
  mif$spatial[[1]]$density_compartment <- "Tumor"
  expect_error(
    split_tissue(mif, classifier = "Classifier.Label", class1 = "Tumor", class2 = "Stroma",
                sigma = 40, interface_width = 100),
    "already been run"
  )
})

# ---- Errors and edges -----------------------------------------------------------

test_that("sigma and interface_width have no default and error when missing", {
  mif <- halfplane_mif()
  expect_error(
    split_tissue(mif, classifier = "Classifier.Label", class1 = "Tumor", class2 = "Stroma",
                interface_width = 100),
    "no default"
  )
  expect_error(
    split_tissue(mif, classifier = "Classifier.Label", class1 = "Tumor", class2 = "Stroma",
                sigma = 40),
    "interface_width"
  )
})

test_that("class1 == class2 and an \"Interface\"-named class both error", {
  mif <- halfplane_mif()
  expect_error(
    split_tissue(mif, classifier = "Classifier.Label", class1 = "Tumor", class2 = "Tumor",
                sigma = 40, interface_width = 100),
    "must differ"
  )
  expect_error(
    split_tissue(mif, classifier = "Classifier.Label", class1 = "Interface", class2 = "Stroma",
                sigma = 40, interface_width = 100),
    "reserved"
  )
})

test_that("a missing classifier column errors naming the sample", {
  mif <- halfplane_mif()
  expect_error(
    split_tissue(mif, classifier = "NoSuchColumn", class1 = "Tumor", class2 = "Stroma",
                sigma = 40, interface_width = 100),
    "not found"
  )
})

test_that("NA in the classifier column does not error (the pp[keep] index trap)", {
  mif <- halfplane_mif()
  mif$spatial[[1]]$Classifier.Label[1:5] <- NA
  out <- split_tissue(mif, classifier = "Classifier.Label", class1 = "Tumor", class2 = "Stroma",
                      sigma = 40, interface_width = 100)
  expect_s3_class(out, "mif")
})

test_that("both classes absent warns and gives an empty boundary", {
  mif <- halfplane_mif()
  mif$spatial[[1]]$Classifier.Label <- "Other"
  expect_warning(
    out <- split_tissue(mif, classifier = "Classifier.Label", class1 = "Tumor", class2 = "Stroma",
                        sigma = 40, interface_width = 100),
    "will be empty"
  )
  expect_equal(nrow(out$derived$density_boundary[[1]]), 0)
  expect_equal(out$sample$`Boundary Length`[1], 0)
})

test_that("filter_density = identity is bit-identical to omitting it", {
  mif <- halfplane_mif()
  a <- split_tissue(mif, classifier = "Classifier.Label", class1 = "Tumor", class2 = "Stroma",
                    sigma = 40, interface_width = 100)
  b <- split_tissue(mif, classifier = "Classifier.Label", class1 = "Tumor", class2 = "Stroma",
                    sigma = 40, interface_width = 100, overwrite = TRUE,
                    filter_density = function(im) im)
  expect_identical(a$derived$density_boundary[[1]][, c("piece", "x", "y")],
                   b$derived$density_boundary[[1]][, c("piece", "x", "y")])
  expect_identical(a$sample$`Boundary Length`, b$sample$`Boundary Length`)
})

test_that("filter_density affects only geometry, never the per-cell sign", {
  mif <- halfplane_mif()
  unfiltered <- split_tissue(mif, classifier = "Classifier.Label",
                             class1 = "Tumor", class2 = "Stroma",
                             sigma = 40, interface_width = 100)
  mask_half <- function(im) {
    im$v[im$xcol < 500, ] <- NA
    im
  }
  filtered <- split_tissue(mif, classifier = "Classifier.Label",
                           class1 = "Tumor", class2 = "Stroma",
                           sigma = 40, interface_width = 100, overwrite = TRUE,
                           filter_density = mask_half)
  expect_identical(unfiltered$spatial[[1]]$density_compartment,
                   filtered$spatial[[1]]$density_compartment)
  expect_lt(filtered$sample$`Boundary Length`[1], unfiltered$sample$`Boundary Length`[1])
})

test_that("filter_density returning a different grid or a non-im errors", {
  mif <- halfplane_mif()
  wrong_grid <- function(im) spatstat.geom::im(im$v[-1, -1], xcol = im$xcol[-1], yrow = im$yrow[-1])
  expect_error(
    split_tissue(mif, classifier = "Classifier.Label", class1 = "Tumor", class2 = "Stroma",
                sigma = 40, interface_width = 100, filter_density = wrong_grid),
    "same pixel grid"
  )
  not_im <- function(im) im$v
  expect_error(
    split_tissue(mif, classifier = "Classifier.Label", class1 = "Tumor", class2 = "Stroma",
                sigma = 40, interface_width = 100, filter_density = not_im),
    "class \"im\""
  )
})

test_that("an unrecognised argument in ... is an error naming it", {
  mif <- halfplane_mif()
  expect_error(
    split_tissue(mif, classifier = "Classifier.Label", class1 = "Tumor", class2 = "Stroma",
                sigma = 40, interface_width = 100, workerss = 4),
    "Unknown argument.*workerss"
  )
})

test_that("too small a sigma hits the pixel budget and names a minimum sigma", {
  mif <- halfplane_mif()
  expect_error(
    split_tissue(mif, classifier = "Classifier.Label", class1 = "Tumor", class2 = "Stroma",
                sigma = 1, interface_width = 10),
    "pixel density grid"
  )
})

test_that("explicit xloc/yloc naming XMin/YMin match the XMax/XMin midpoint default", {
  mif <- halfplane_mif()   # XMin == XMax, YMin == YMax by construction
  default <- split_tissue(mif, classifier = "Classifier.Label",
                          class1 = "Tumor", class2 = "Stroma",
                          sigma = 40, interface_width = 100)
  explicit <- split_tissue(mif, classifier = "Classifier.Label",
                           class1 = "Tumor", class2 = "Stroma",
                           sigma = 40, interface_width = 100, overwrite = TRUE,
                           xloc = "XMin", yloc = "YMin")
  expect_identical(default$derived$density_boundary[[1]][, c("piece", "x", "y")],
                   explicit$derived$density_boundary[[1]][, c("piece", "x", "y")])
})

test_that("workers = 2 gives identical results to workers = 1 (no RNG involved)", {
  mif <- halfplane_mif()
  one <- split_tissue(mif, classifier = "Classifier.Label", class1 = "Tumor", class2 = "Stroma",
                      sigma = 40, interface_width = 100, workers = 1)
  two <- split_tissue(mif, classifier = "Classifier.Label", class1 = "Tumor", class2 = "Stroma",
                      sigma = 40, interface_width = 100, workers = 2, overwrite = TRUE)
  expect_identical(one$spatial, two$spatial)
})

test_that("a mif with no spatial data errors clearly, not with a subscript failure", {
  # create_mif(spatial_list = NULL) itself rejects NULL before it can reach the
  # list(NA) fallback (is.list(NULL) is FALSE) -- construct that shape directly
  # to test split_tissue()'s guard against it.
  mif <- structure(
    list(clinical = data.frame(deidentified_id = "1"),
        sample = data.frame(deidentified_id = "1", deidentified_sample = "S1"),
        spatial = list(NA), derived = list(),
        patient_id = "deidentified_id", sample_id = "deidentified_sample"),
    class = "mif")
  expect_error(
    split_tissue(mif, classifier = "Classifier.Label", class1 = "Tumor", class2 = "Stroma",
                sigma = 40, interface_width = 100),
    "no spatial data"
  )
})

test_that("split_tissue rejects a non-mif", {
  expect_error(
    split_tissue(list(), classifier = "x", class1 = "a", class2 = "b",
                sigma = 40, interface_width = 100),
    "class `mif`"
  )
})
