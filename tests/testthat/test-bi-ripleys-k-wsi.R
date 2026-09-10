test_that("bi_ripleys_k_WSI errors with a pointer to bi_ripleys_k", {
  # Kept exported rather than removed so the name still resolves and the error can
  # explain the migration, instead of "could not find function".
  mif <- example_mif()
  expect_error(
    bi_ripleys_k_WSI(mif, mnames = mnames_bivariate(), r_range = seq(0, 40, 10)),
    "removed in spatialTIME 2\\.0\\.0"
  )
  expect_error(
    bi_ripleys_k_WSI(mif, mnames = mnames_bivariate()),
    "bi_ripleys_k\\(\\)"
  )
  # The old size arguments are accepted by the signature -- so a legacy call reaches
  # the explanatory error rather than an argument-matching failure.
  expect_error(
    bi_ripleys_k_WSI(mif, mnames = mnames_bivariate(), big = 1000, nlarge = 1000),
    "removed in spatialTIME"
  )
})

test_that("bi_ripleys_k handles what the WSI variant used to be needed for", {
  # The whole reason WSI existed was that bi_ripleys_k built an n x n matrix. It now
  # materialises only pairs within max(r_range), so a cell count that would have
  # required the tiled path works directly -- with the requested edge correction,
  # which the tiled path silently replaced with "none".
  set.seed(4)
  n <- 20000
  spat <- data.frame(deidentified_sample = "WSI",
                     XMin = runif(n, 0, 5000), YMin = runif(n, 0, 5000))
  spat$XMax <- spat$XMin; spat$YMax <- spat$YMin
  set.seed(5)
  spat$A <- rbinom(n, 1, 0.05); spat$B <- rbinom(n, 1, 0.05)
  spat$A[spat$B == 1] <- 0L                       # keep the sets disjoint
  mif <- toy_mif(list(WSI = spat))

  out <- bi_ripleys_k(mif, mnames = c("A", "B"), r_range = seq(0, 50, 10),
                      edge_correction = "translation", permute = FALSE,
                      workers = 1, big = 5000, overwrite = TRUE)
  tbl <- out$derived$bivariate_Count
  expect_true(any(!is.na(tbl$`Observed K`)))
  # And chunking is a memory knob only: it must not change the numbers.
  unchunked <- bi_ripleys_k(mif, mnames = c("A", "B"), r_range = seq(0, 50, 10),
                            edge_correction = "translation", permute = FALSE,
                            workers = 1, big = 1e9, overwrite = TRUE)
  expect_equal(tbl, unchunked$derived$bivariate_Count)
})
