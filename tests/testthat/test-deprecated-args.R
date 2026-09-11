# Renamed arguments keep working through the 2.x line, with a warning naming the
# replacement. Equally important: adding `...` to a user-facing function must not
# turn typos into silent no-ops, so anything in `...` that is not a recognised
# former argument is an error.

dep_fixture <- function() {
  mif <- create_mif(
    clinical_data = example_clinical %>%
      dplyr::mutate(deidentified_id = as.character(deidentified_id)),
    sample_data = example_summary %>%
      dplyr::mutate(deidentified_id = as.character(deidentified_id)),
    spatial_list = example_spatial["TMA3_[9,K].tif"],
    patient_id = "deidentified_id", sample_id = "deidentified_sample")
  list(mif = mif, mnames = c("CD3..Opal.570..Positive", "CD8..Opal.520..Positive"),
       r = seq(0, 60, 10))
}

test_that("keep_perm_dis warns and forwards to keep_permutation_distribution", {
  f <- dep_fixture()
  expect_warning(
    out <- NN_G(f$mif, mnames = f$mnames[1], r_range = f$r, num_permutations = 4,
                workers = 1, overwrite = TRUE, keep_perm_dis = TRUE),
    "`keep_perm_dis` is deprecated"
  )
  # Forwarded, not just warned about: TRUE means one row per permutation per radius.
  expect_equal(nrow(out$derived$univariate_NN), 4 * length(f$r))

  expect_warning(
    out2 <- NN_G(f$mif, mnames = f$mnames[1], r_range = f$r, num_permutations = 4,
                 workers = 1, overwrite = TRUE, keep_perm_dis = FALSE),
    "keep_permutation_distribution"
  )
  expect_equal(nrow(out2$derived$univariate_NN), length(f$r))
})

test_that("keep_perm_dis is accepted by every metric that used to have it", {
  f <- dep_fixture()
  calls <- list(
    function() ripleys_k(f$mif, mnames = f$mnames[1], r_range = f$r, permute = TRUE,
                         num_permutations = 2, workers = 1, overwrite = TRUE,
                         keep_perm_dis = FALSE),
    function() bi_ripleys_k(f$mif, mnames = f$mnames, r_range = f$r, permute = TRUE,
                            num_permutations = 2, workers = 1, keep_perm_dis = FALSE),
    function() NN_G(f$mif, mnames = f$mnames[1], r_range = f$r, num_permutations = 2,
                    workers = 1, overwrite = TRUE, keep_perm_dis = FALSE),
    function() bi_NN_G(f$mif, mnames = f$mnames, r_range = f$r, num_permutations = 2,
                       workers = 1, overwrite = TRUE, keep_perm_dis = FALSE)
  )
  for (i in seq_along(calls)) {
    expect_warning(res <- calls[[i]](), "keep_perm_dis", info = i)
    expect_s3_class(res, "mif")
  }
})

test_that("defunct arguments warn, are ignored, and do not change results", {
  f <- dep_fixture()
  clean <- ripleys_k(f$mif, mnames = f$mnames[1], r_range = f$r, permute = FALSE,
                     workers = 1, overwrite = TRUE)$derived$univariate_Count
  expect_warning(
    withmethod <- ripleys_k(f$mif, mnames = f$mnames[1], r_range = f$r, permute = FALSE,
                            workers = 1, overwrite = TRUE, method = "K"),
    "`method` is defunct"
  )
  expect_equal(withmethod$derived$univariate_Count, clean)

  # `force` existed only to get past a 10,000-cell stop that no longer exists.
  expect_warning(
    bi_ripleys_k(f$mif, mnames = f$mnames, r_range = f$r, permute = FALSE,
                 workers = 1, force = TRUE),
    "`force` is defunct"
  )
})

test_that("nlarge forwards to big for bi_ripleys_k", {
  f <- dep_fixture()
  expect_warning(
    out <- bi_ripleys_k(f$mif, mnames = f$mnames, r_range = f$r, permute = FALSE,
                        workers = 1, nlarge = 500),
    "`nlarge` is deprecated as of spatialTIME 2\\.0\\.0; use `big`"
  )
  # big only bounds memory, so the numbers must equal an unchunked run.
  ref <- bi_ripleys_k(f$mif, mnames = f$mnames, r_range = f$r, permute = FALSE,
                      workers = 1, big = 1e9)
  expect_equal(out$derived$bivariate_Count, ref$derived$bivariate_Count)
})

test_that("unknown arguments in ... are an error, not silently ignored", {
  # The whole risk of adding `...`. A typo must fail loudly.
  f <- dep_fixture()
  expect_error(
    ripleys_k(f$mif, mnames = f$mnames[1], r_range = f$r, workers = 1, workerss = 2),
    "Unknown argument.*workerss"
  )
  expect_error(
    NN_G(f$mif, mnames = f$mnames[1], r_range = f$r, workers = 1, edge_corection = "rs"),
    "Unknown argument.*edge_corection"
  )
  # An unnamed extra only reaches `...` once every earlier formal is named --
  # otherwise R matches it positionally to the first free formal.
  expect_error(
    bi_NN_G(mif = f$mif, mnames = f$mnames, r_range = f$r, num_permutations = 2,
            edge_correction = "rs", keep_permutation_distribution = FALSE,
            workers = 1, overwrite = TRUE, xloc = NULL, yloc = NULL, TRUE),
    "unnamed arguments"
  )
})

test_that("pcf functions still forward genuine extra arguments to spatstat", {
  # pair_correlation()/bi_pair_correlation() use `...` for BOTH deprecated names
  # and real pcf arguments, so the split between the two must work.
  f <- dep_fixture()
  a <- suppressWarnings(pair_correlation(f$mif, mnames = f$mnames[1], r_range = f$r,
                        num_permutations = 2, workers = 1, overwrite = TRUE))
  # "epanechnikov" is the default in spatstat.explore >= 3.8-1, so pick another
  # that genuinely differs from it.
  b <- suppressWarnings(pair_correlation(f$mif, mnames = f$mnames[1], r_range = f$r,
                        num_permutations = 2, workers = 1, overwrite = TRUE,
                        kernel = "gaussian"))
  expect_s3_class(b, "mif")
  # A different smoothing kernel must actually change the estimate.
  expect_false(isTRUE(all.equal(a$derived$univariate_pair_correlation$`Observed g`,
                                b$derived$univariate_pair_correlation$`Observed g`)))
})

test_that("split_tissue/plot_tissue_split map to no deprecated args, so keep_perm_dis errors", {
  # Both map to character(0) in deprecated_arg_map(), by explicit branches that
  # must come BEFORE the fall-through `common` default -- otherwise
  # keep_perm_dis would be silently accepted here too.
  expect_identical(deprecated_arg_map("split_tissue"), character(0))
  expect_identical(deprecated_arg_map("plot_tissue_split"), character(0))

  mif <- toy_mif(list(S1 = toy_spatial("S1", n = 100, markers = c(A = 20, B = 20))))
  expect_error(
    split_tissue(mif, classifier = "Classifier.Label", class1 = "Tumor",
                class2 = "Stroma", sigma = 40, interface_width = 50,
                keep_perm_dis = TRUE),
    "Unknown argument.*keep_perm_dis"
  )

  split <- split_tissue(mif, classifier = "Classifier.Label", class1 = "Tumor",
                        class2 = "Stroma", sigma = 40, interface_width = 50)
  expect_error(plot_tissue_split(split, keep_perm_dis = TRUE),
              "Unknown argument.*keep_perm_dis")
})

test_that("the deprecation map covers every function that takes dots", {
  for (fn in c("ripleys_k", "bi_ripleys_k", "NN_G", "bi_NN_G",
               "pair_correlation", "bi_pair_correlation", "interaction_variable")) {
    m <- deprecated_arg_map(fn)
    expect_true("keep_perm_dis" %in% names(m), info = fn)
    expect_identical(m[["keep_perm_dis"]], "keep_permutation_distribution", info = fn)
  }
})
