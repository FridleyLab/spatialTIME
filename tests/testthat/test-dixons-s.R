skip_if_not_installed("dixon")

# dixons_s() had no test, which is how overwrite = FALSE -- its DEFAULT -- came to
# always error. Run-semantics live in test-run-semantics.R; this file covers the
# table schemas and the guards.

dix_mif <- function() {
  toy_mif(list(S1 = toy_spatial("S1", n = 200, markers = c(A = 40, B = 40)),
               S2 = toy_spatial("S2", n = 200, markers = c(A = 35, B = 45), seed = 12)))
}

test_that("dixons_s fills both Dixon slots and returns a mif", {
  out <- dixons_s(dix_mif(), mnames = c("A", "B"), num_permutations = 5,
                  workers = 1, overwrite = TRUE)
  expect_s3_class(out, "mif")
  expect_true(all(c("Dixon_Z", "Dixon_C") %in% names(out$derived)))
  expect_gt(nrow(out$derived$Dixon_Z), 0)
  expect_gt(nrow(out$derived$Dixon_C), 0)
})

test_that("type selects which slots are written", {
  z <- dixons_s(dix_mif(), mnames = c("A", "B"), num_permutations = 5, workers = 1,
                type = "Z", overwrite = TRUE)
  expect_true("Dixon_Z" %in% names(z$derived))
  expect_false("Dixon_C" %in% names(z$derived))

  c_only <- dixons_s(dix_mif(), mnames = c("A", "B"), num_permutations = 5, workers = 1,
                     type = "C", overwrite = TRUE)
  expect_true("Dixon_C" %in% names(c_only$derived))
  expect_false("Dixon_Z" %in% names(c_only$derived))
})

test_that("an invalid type is rejected instead of silently doing nothing", {
  # `type` was never validated, so type = "z" ran the whole permutation loop and
  # returned the mif unchanged with no message.
  expect_error(dixons_s(dix_mif(), mnames = c("A", "B"), type = "z"), "must be")
  expect_error(dixons_s(dix_mif(), mnames = c("A", "B"), type = character(0)), "must be")
})

test_that("column names are clean and pair-identifiable", {
  # dixon::dixon()'s own tablaC names carry padding spaces ("  df ", "  P.rand") and
  # only tablaZ used to be de-spaced, so binding a success row to an early-return row
  # produced BOTH "  P.rand" and "P.rand".
  out <- dixons_s(dix_mif(), mnames = c("A", "B"), num_permutations = 5,
                  workers = 1, overwrite = TRUE)
  for (slot in c("Dixon_Z", "Dixon_C")) {
    nms <- names(out$derived[[slot]])
    expect_equal(length(grep("^\\s|\\s$", nms)), 0, info = slot)
    expect_false(any(duplicated(nms)), info = slot)
    # Per-marker count columns are pair-specific, so a pair identifier is required
    # for the table to be interpretable at all.
    expect_true(all(c("Anchor", "Counted") %in% nms), info = slot)
  }
  expect_true("df" %in% names(out$derived$Dixon_C))
})

test_that("the schema does not depend on whether a pair hit the sparse guard", {
  # Marker pairs with fewer than 3 cells of either type return NA statistics rather
  # than a differently shaped frame, so bind_rows across pairs cannot go ragged.
  dense  <- toy_mif(list(S1 = toy_spatial("S1", n = 200, markers = c(A = 40, B = 40))))
  sparse <- toy_mif(list(S1 = toy_spatial("S1", n = 200, markers = c(A = 40, B = 1))))

  a <- dixons_s(dense,  mnames = c("A", "B"), num_permutations = 5, workers = 1, overwrite = TRUE)
  b <- dixons_s(sparse, mnames = c("A", "B"), num_permutations = 5, workers = 1, overwrite = TRUE)

  expect_identical(names(a$derived$Dixon_Z), names(b$derived$Dixon_Z))
  expect_identical(names(a$derived$Dixon_C), names(b$derived$Dixon_C))
  # The sparse pair reports NA rather than being dropped.
  expect_true(all(is.na(b$derived$Dixon_Z$Z)))
})

test_that("a mixed dense/sparse run gives one consistent table", {
  # Three markers: A/B is estimable, anything with C is not. This is the case that
  # used to yield both "P.rand" and "  P.rand" columns.
  mif <- toy_mif(list(S1 = toy_spatial("S1", n = 240,
                                       markers = c(A = 40, B = 40, C = 1))))
  out <- dixons_s(mif, mnames = c("A", "B", "C"), num_permutations = 5,
                  workers = 1, overwrite = TRUE)
  z <- out$derived$Dixon_Z
  expect_false(any(duplicated(names(z))))
  # Both an estimable and a non-estimable pair are present.
  pairs <- unique(paste(z$Anchor, z$Counted))
  expect_gt(length(pairs), 1)
  expect_true(any(!is.na(z$Z)))
  expect_true(any(is.na(z$Z)))
})

test_that("dixons_s validates its inputs", {
  expect_error(dixons_s(list(), mnames = c("A", "B")), "class `mif`")
  expect_error(dixons_s(dix_mif(), mnames = 1), "vector of marker names")
  # 1:nrow(mnames) on a zero-row combination frame used to iterate c(1, 0) and fail
  # cryptically inside the loop.
  expect_error(dixons_s(dix_mif(), mnames = "A"), "at least 2 markers")
  expect_error(dixons_s(dix_mif(), mnames = c("A", "NotAColumn"), num_permutations = 5),
               "not found in spatial data")
})

test_that("dixons_s prints nothing to stdout", {
  # dixon::dixon() writes a permutation counter to stdout; this package strips its
  # own progress output and should not leak a dependency's either.
  out <- capture.output(invisible(
    dixons_s(dix_mif(), mnames = c("A", "B"), num_permutations = 5,
             workers = 1, overwrite = TRUE)))
  expect_length(out, 0)
})
