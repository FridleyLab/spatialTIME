two_mifs <- function() {
  list(
    toy_mif(list(S1 = toy_spatial("S1"), S2 = toy_spatial("S2", seed = 12))),
    toy_mif(list(S3 = toy_spatial("S3", seed = 13), S4 = toy_spatial("S4", seed = 14)),
            patient_ids = c("3", "4"))
  )
}

test_that("merging disjoint mifs unions every slot", {
  m <- two_mifs()
  out <- merge_mifs(m)
  expect_s3_class(out, "mif")
  expect_length(out$spatial, 4)
  expect_setequal(names(out$spatial), c("S1", "S2", "S3", "S4"))
  expect_equal(nrow(out$clinical), 4)
  expect_identical(out$patient_id, "deidentified_id")
  expect_identical(out$sample_id, "deidentified_sample")
})

test_that("fewer than two mifs is rejected", {
  # length(mifs) == 1 used to be the only guard, so an empty list fell through to
  # names(sizes) = seq(length(sizes)) -> seq(0) is c(1, 0), length 2, and died on
  # "'names' attribute [2] must be the same length as the vector [0]".
  expect_error(merge_mifs(list()), "at least 2 MIF objects")
  expect_error(merge_mifs(NULL), "at least 2 MIF objects")
  expect_error(merge_mifs(two_mifs()[1]), "at least 2 MIF objects")
})

test_that("non-mif elements are rejected", {
  expect_error(merge_mifs(list(two_mifs()[[1]], list(a = 1))), "must be a mif object")
})

test_that("duplicate spatial names and clashing ids are caught by check.names", {
  m <- two_mifs()
  expect_error(merge_mifs(list(m[[1]], m[[1]])), "same name")

  clash <- m[[2]]
  clash$patient_id <- "other_id"
  clash$sample[["other_id"]] <- clash$sample[["deidentified_id"]]
  clash$clinical[["other_id"]] <- clash$clinical[["deidentified_id"]]
  expect_error(merge_mifs(list(m[[1]], clash)), "different patient_id")
})

test_that("derived slots with differing columns merge without dropping data", {
  m <- two_mifs()
  a <- ripleys_k(m[[1]], mnames = "A", r_range = seq(0, 40, 10), permute = FALSE,
                 workers = 1, overwrite = TRUE)
  b <- NN_G(m[[2]], mnames = "A", r_range = seq(0, 40, 10), num_permutations = 2,
            workers = 1, overwrite = TRUE)
  out <- merge_mifs(list(a, b))
  # Different metrics -> different slots, both preserved.
  expect_true(all(c("univariate_Count", "univariate_NN") %in% names(out$derived)))
  expect_equal(nrow(out$derived$univariate_Count), nrow(a$derived$univariate_Count))
  expect_equal(nrow(out$derived$univariate_NN), nrow(b$derived$univariate_NN))
})

test_that("the same derived slot from two mifs is row-bound", {
  m <- two_mifs()
  a <- ripleys_k(m[[1]], mnames = "A", r_range = seq(0, 40, 10), permute = FALSE,
                 workers = 1, overwrite = TRUE)
  b <- ripleys_k(m[[2]], mnames = "A", r_range = seq(0, 40, 10), permute = FALSE,
                 workers = 1, overwrite = TRUE)
  out <- merge_mifs(list(a, b))
  expect_equal(nrow(out$derived$univariate_Count),
               nrow(a$derived$univariate_Count) + nrow(b$derived$univariate_Count))
  expect_setequal(unique(out$derived$univariate_Count$deidentified_sample),
                  c("S1", "S2", "S3", "S4"))
})

test_that("merging mifs with no derived data reports it and still works", {
  expect_message(out <- merge_mifs(two_mifs()), "No variables have been derived")
  expect_identical(out$derived, stats::setNames(list(), character(0)))
})
