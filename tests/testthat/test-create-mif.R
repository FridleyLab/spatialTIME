# Replaces test_creatMIF.R (misspelled, and its four assertions were
# near-tautological: they compared mif$clinical / mif$sample / mif$spatial to the
# unmodified inputs, which passes for any implementation that stores its arguments).

test_that("create_mif returns the six documented slots", {
  mif <- example_mif()
  expect_s3_class(mif, "mif")
  expect_named(mif, c("clinical", "sample", "spatial", "derived",
                      "patient_id", "sample_id"))
  expect_identical(mif$patient_id, "deidentified_id")
  expect_identical(mif$sample_id, "deidentified_sample")
  expect_identical(mif$derived, list())
  expect_true(is.list(mif$spatial))
})

test_that("create_mif accepts integer patient ids", {
  # Before 2.0.0 a full_join was performed on the ids and then discarded, so
  # mismatched integer/character ids could make create_mif fail over a result it
  # never used. That is what every example's as.character() call was working
  # around.
  sp <- list(S1 = toy_spatial("S1"))
  expect_s3_class(
    create_mif(
      clinical_data = data.frame(deidentified_id = 1L),
      sample_data   = data.frame(deidentified_id = 1L, deidentified_sample = "S1"),
      spatial_list  = sp,
      patient_id = "deidentified_id", sample_id = "deidentified_sample"),
    "mif"
  )
})

test_that("create_mif recovers spatial names from the sample_id column", {
  sp <- list(toy_spatial("A"), toy_spatial("B"))   # deliberately unnamed
  mif <- create_mif(
    clinical_data = data.frame(deidentified_id = c("1", "2")),
    sample_data   = data.frame(deidentified_id = c("1", "2"),
                               deidentified_sample = c("A", "B")),
    spatial_list  = sp,
    patient_id = "deidentified_id", sample_id = "deidentified_sample")
  expect_identical(names(mif$spatial), c("A", "B"))
})

test_that("create_mif validates its inputs", {
  ok_clin <- data.frame(deidentified_id = "1")
  ok_samp <- data.frame(deidentified_id = "1", deidentified_sample = "S1")
  ok_sp   <- list(S1 = toy_spatial("S1"))
  # Note: not utils::modifyList() -- it RECURSES into list arguments rather than
  # replacing them, so an unnamed replacement for spatial_list silently no-ops.
  mk <- function(...) {
    args <- list(clinical_data = ok_clin, sample_data = ok_samp, spatial_list = ok_sp,
                 patient_id = "deidentified_id", sample_id = "deidentified_sample")
    over <- list(...)
    args[names(over)] <- over
    do.call(create_mif, args)
  }
  expect_error(mk(clinical_data = "nope"), "clinical_data must be a data frame")
  expect_error(mk(sample_data = "nope"), "sample_data must be a data frame")
  expect_error(mk(spatial_list = list("nope")), "must be a data frame")
  expect_error(mk(patient_id = 1), "patient_id must be a character")
  expect_error(mk(patient_id = "absent"), "could not be found in 'clinical_data'")
  expect_error(mk(sample_id = "absent"), "could not be found in 'sample_data'")
})

test_that("create_mif rejects a partially named spatial list", {
  # The original check `all(!sapply(names(spatial_list), is.null))` could never
  # fail: names() returns a character vector, so is.null() per element is always
  # FALSE, and when names() is NULL sapply gives list() and all(logical(0)) is TRUE.
  sp <- list(A = toy_spatial("A"), toy_spatial("B"))
  names(sp) <- c("A", "")
  expect_error(
    create_mif(clinical_data = data.frame(deidentified_id = "1"),
               sample_data = data.frame(deidentified_id = "1", deidentified_sample = "A"),
               spatial_list = sp,
               patient_id = "deidentified_id", sample_id = "deidentified_sample"),
    "must be named"
  )
})

test_that("print.mif reports counts, returns its argument invisibly", {
  mif <- example_mif()
  withr::local_options(crayon.enabled = FALSE)
  expect_output(print(mif), "patients spanning")
  expect_output(print(mif), "spatial data frames were found")
  # A print method must return invisible(x); before 2.0.0 this returned NULL, so
  # `y <- print(mif)` gave NULL and print() could not sit in a pipe.
  res <- NULL
  invisible(capture.output(res <- withVisible(print(mif))))
  expect_false(res$visible)
  expect_identical(res$value, mif)
})

test_that("print.mif survives NULL and empty slots", {
  # NULL[["id"]] is an error in R, not NULL, so a partially built mif used to be
  # impossible even to echo at the console.
  withr::local_options(crayon.enabled = FALSE)
  mif <- example_mif()

  for (slot in c("clinical", "sample", "patient_id", "sample_id")) {
    broken <- mif
    broken[[slot]] <- NULL
    expect_output(print(broken), "spatial data frames were found",
                  info = slot)
  }
  empty <- mif
  empty$spatial <- list()
  expect_output(print(empty), "0 spatial data frames")
})

test_that("print.mif lists derived slots once they exist", {
  withr::local_options(crayon.enabled = FALSE)
  mif <- example_mif()
  expect_output(print(mif), "spatial data frames were found")
  filled <- ripleys_k(mif, mnames = mnames_good()[1], r_range = seq(0, 40, 10),
                      permute = FALSE, workers = 1, overwrite = TRUE)
  expect_output(print(filled), "univariate_Count")
})
