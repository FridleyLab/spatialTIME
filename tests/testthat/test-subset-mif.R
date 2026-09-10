test_that("subset_mif keeps only cells at the requested level", {
  mif <- example_mif(which = 1:3, n_cells = 300)
  out <- subset_mif(mif, classifier = "Classifier.Label", level = "Tumor",
                    markers = mnames_good()[1:2])
  expect_s3_class(out, "mif")
  for (nm in names(out$spatial)) {
    expect_true(all(out$spatial[[nm]]$Classifier.Label == "Tumor"), info = nm)
  }
})

test_that("summary rows stay aligned with retained spatial frames", {
  # THE regression. `out` was assigned only inside if(nrow(tmp) > 2) while
  # rbind.data.frame(summary, t(out)) ran unconditionally, so a sample with <= 2
  # cells at the requested level either killed the function (if it was first) or
  # silently re-appended the PREVIOUS sample's row (if it was later), leaving a
  # summary row that matched no retained spatial frame.
  mk <- c("A", "B")

  # first sample degenerate
  sp <- list(S1 = toy_spatial("S1", n = 100, markers = c(A = 20, B = 20)),
             S2 = toy_spatial("S2", n = 100, markers = c(A = 20, B = 20), seed = 12),
             S3 = toy_spatial("S3", n = 100, markers = c(A = 20, B = 20), seed = 13))
  sp$S1$Classifier.Label <- "Stroma"
  sp$S1$Classifier.Label[1] <- "Tumor"           # exactly 1 Tumor cell
  out <- subset_mif(toy_mif(sp), classifier = "Classifier.Label", level = "Tumor",
                    markers = mk)
  expect_equal(length(out$spatial), 2)
  expect_equal(nrow(out$sample), length(out$spatial))
  expect_setequal(out$sample$deidentified_sample, names(out$spatial))

  # middle sample degenerate -- this is the case that used to corrupt silently
  sp2 <- sp
  sp2$S1 <- toy_spatial("S1", n = 100, markers = c(A = 20, B = 20))
  sp2$S2$Classifier.Label <- "Stroma"
  sp2$S2$Classifier.Label[1:2] <- "Tumor"        # exactly 2 Tumor cells
  out2 <- subset_mif(toy_mif(sp2), classifier = "Classifier.Label", level = "Tumor",
                     markers = mk)
  expect_equal(nrow(out2$sample), length(out2$spatial))
  expect_false(any(duplicated(out2$sample$deidentified_sample)))
  expect_false("S2" %in% out2$sample$deidentified_sample)
})

test_that("every sample degenerate yields an empty but valid mif", {
  sp <- list(S1 = toy_spatial("S1", n = 50, markers = c(A = 10, B = 10)))
  sp$S1$Classifier.Label <- "Stroma"
  expect_warning(
    out <- subset_mif(toy_mif(sp), classifier = "Classifier.Label", level = "Tumor",
                      markers = c("A", "B")),
    "No sample had more than 2 cells"
  )
  expect_s3_class(out, "mif")
  expect_length(out$spatial, 0)
  expect_equal(nrow(out$sample), 0)
  # The zero-row summary must still carry the id columns, or create_mif rejects it.
  expect_true(all(c("deidentified_id", "deidentified_sample") %in% names(out$sample)))
})

test_that("the summary table is typed, not all character", {
  # c(patient, id, unlist(counts), unlist(percent)) coerced everything to character
  # at the c(), so counts shipped as "3611" and proportions as
  # "0.00249238438105788".
  mif <- example_mif(which = 1:2, n_cells = 300)
  out <- subset_mif(mif, classifier = "Classifier.Label", level = "Tumor",
                    markers = mnames_good()[1:2])
  value_cols <- setdiff(names(out$sample), c("deidentified_id", "deidentified_sample"))
  expect_true(all(vapply(out$sample[value_cols], is.numeric, logical(1))))
})

test_that("the % columns are percentages, matching marker_freq_diff", {
  # Before 2.0.0 these held sum/nrow -- a proportion -- while marker_freq_diff put
  # true percentages in its "%" columns, so the two exports disagreed on what "%"
  # meant. This asserts the 100x fix.
  mk <- "A"
  sp <- list(S1 = toy_spatial("S1", n = 200, markers = c(A = 50)))
  # Every cell is Tumor so the arithmetic is unambiguous.
  sp$S1$Classifier.Label <- "Tumor"
  out <- subset_mif(toy_mif(sp), classifier = "Classifier.Label", level = "Tumor",
                    markers = mk)
  expect_equal(out$sample[["Tumor: A"]], 50)
  expect_equal(out$sample[["Tumor: Total Cells"]], 200)
  expect_equal(out$sample[["Tumor: % A"]], 25)     # not 0.25
})

test_that("a duplicated sample id does not shift the row", {
  # `patient` came from an unguarded == lookup, so more than one matching row in
  # mif$sample made it length > 1 and silently pushed every later value along.
  sp <- list(S1 = toy_spatial("S1", n = 100, markers = c(A = 20, B = 20)))
  mif <- toy_mif(sp)
  mif$sample <- rbind(mif$sample, mif$sample)       # duplicate the row
  out <- subset_mif(mif, classifier = "Classifier.Label", level = "Tumor",
                    markers = c("A", "B"))
  expect_equal(nrow(out$sample), 1)
  expect_identical(out$sample$deidentified_sample, "S1")
  expect_true(is.numeric(out$sample[["Tumor: Total Cells"]]))
})

test_that("subset_mif rejects a non-mif", {
  expect_error(subset_mif(list(), classifier = "x", level = "y", markers = "z"),
               "class `mif`")
})
