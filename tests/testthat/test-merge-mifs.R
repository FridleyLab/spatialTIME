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

test_that("list-valued derived slots (density_boundary) concatenate, not bind_rows", {
  # Before the type-check guard, bind_rows() on a list keyed by sample name
  # silently collapsed it into a nameless table. split_tissue()'s
  # density_boundary must instead end up as a list of 4 named, per-sample frames.
  m <- two_mifs()
  a <- split_tissue(m[[1]], classifier = "Classifier.Label", class1 = "Tumor",
                    class2 = "Stroma", sigma = 40, interface_width = 50)
  b <- split_tissue(m[[2]], classifier = "Classifier.Label", class1 = "Tumor",
                    class2 = "Stroma", sigma = 40, interface_width = 50)
  out <- merge_mifs(list(a, b))
  bd <- out$derived$density_boundary
  expect_type(bd, "list")
  expect_false(is.data.frame(bd))
  expect_setequal(names(bd), c("S1", "S2", "S3", "S4"))
  for (nm in names(bd)) {
    expect_true(is.data.frame(bd[[nm]]), info = nm)
  }
})

test_that("duplicate sample names in a list-valued derived slot error under check.names", {
  # Spatial/patient names must stay distinct (or check.names would already
  # reject the merge for that reason); force the clash inside density_boundary
  # itself so the derived-slot guard is what actually fires.
  m <- two_mifs()
  a <- split_tissue(m[[1]], classifier = "Classifier.Label", class1 = "Tumor",
                    class2 = "Stroma", sigma = 40, interface_width = 50)
  b <- split_tissue(m[[2]], classifier = "Classifier.Label", class1 = "Tumor",
                    class2 = "Stroma", sigma = 40, interface_width = 50)
  names(b$derived$density_boundary)[1] <- names(a$derived$density_boundary)[1]
  expect_error(merge_mifs(list(a, b)), "share a sample name")
})

test_that("differing call_info on a list-valued derived slot warns and keeps the first", {
  m <- two_mifs()
  a <- split_tissue(m[[1]], classifier = "Classifier.Label", class1 = "Tumor",
                    class2 = "Stroma", sigma = 40, interface_width = 50)
  b <- split_tissue(m[[2]], classifier = "Classifier.Label", class1 = "Tumor",
                    class2 = "Stroma", sigma = 80, interface_width = 50)
  expect_warning(out <- merge_mifs(list(a, b)), "different settings")
  merged_ci <- attr(out$derived$density_boundary, "call_info")
  expect_equal(merged_ci$sigma, attr(a$derived$density_boundary, "call_info")$sigma)
})

test_that("spatial_plots, a pre-existing list-valued derived slot, also merges cleanly", {
  # merge_mifs() used to bind_rows() this too; the same type-check guard that
  # fixes density_boundary fixes this independently of split_tissue.
  m <- two_mifs()
  a <- m[[1]]; b <- m[[2]]
  a$derived$spatial_plots <- list(S1 = ggplot2::ggplot(), S2 = ggplot2::ggplot())
  b$derived$spatial_plots <- list(S3 = ggplot2::ggplot(), S4 = ggplot2::ggplot())
  out <- merge_mifs(list(a, b))
  expect_type(out$derived$spatial_plots, "list")
  expect_false(is.data.frame(out$derived$spatial_plots))
  expect_setequal(names(out$derived$spatial_plots), c("S1", "S2", "S3", "S4"))
})
