skip_if_not_installed("spatstat.geom")
skip_if_not_installed("spatstat.explore")

# spatstat's Gest/Gcross output column name for each correction. Note that
# Gcross(correction = "rs") returns BOTH an `rs` and a `km` column, so selecting
# the estimate by column POSITION returns the wrong estimator for "km". This map
# exists so the tests, like the package, select by name.
GCOL <- c(rs = "rs", km = "km", han = "han", none = "raw")

g_fixture <- function() {
  spat <- example_spatial[["TMA3_[9,K].tif"]]
  spat$xloc <- (spat$XMax + spat$XMin) / 2
  spat$yloc <- (spat$YMax + spat$YMin) / 2
  list(
    spat = spat,
    win  = spatstat.geom::convexhull.xy(spat$xloc, spat$yloc),
    mif  = create_mif(
      clinical_data = example_clinical %>%
        dplyr::mutate(deidentified_id = as.character(deidentified_id)),
      sample_data = example_summary %>%
        dplyr::mutate(deidentified_id = as.character(deidentified_id)),
      spatial_list = example_spatial["TMA3_[9,K].tif"],
      patient_id = "deidentified_id",
      sample_id = "deidentified_sample"
    ),
    mnames = c("CD3..Opal.570..Positive", "CD8..Opal.520..Positive"),
    r = seq(0, 60, 10)
  )
}

test_that("NN_G matches spatstat Gest for every correction", {
  f <- g_fixture()
  for (ec in names(GCOL)) {
    got <- NN_G(f$mif, mnames = f$mnames, r_range = f$r, num_permutations = 2,
                workers = 1, edge_correction = ec,
                overwrite = TRUE)$derived$univariate_NN
    for (m in f$mnames) {
      pos <- f$spat[[m]] == 1
      ref <- as.data.frame(spatstat.explore::Gest(
        spatstat.geom::ppp(f$spat$xloc[pos], f$spat$yloc[pos],
                           window = f$win, check = FALSE),
        r = f$r, correction = ec))[[GCOL[[ec]]]]
      mine <- got[got$Marker == m, ]
      mine <- mine[order(mine$r), ]
      expect_equal(mine$`Observed G`, ref, tolerance = 1e-12,
                   info = paste(ec, m))
    }
  }
})

test_that("bi_NN_G matches spatstat Gcross for every correction", {
  # This is what justifies deleting ~150 lines of hand-rolled rs/han estimators
  # built on a full as.matrix(dist(...)).
  f <- g_fixture()
  for (ec in names(GCOL)) {
    got <- bi_NN_G(f$mif, mnames = f$mnames, r_range = f$r, num_permutations = 2,
                   workers = 1, edge_correction = ec,
                   overwrite = TRUE)$derived$bivariate_NN
    for (i in 1:2) {
      a <- f$mnames[i]; c2 <- f$mnames[3 - i]
      pa <- f$spat[[a]] == 1 & f$spat[[c2]] != 1
      pc <- f$spat[[c2]] == 1 & f$spat[[a]] != 1
      kp <- pa | pc
      Y <- spatstat.geom::ppp(f$spat$xloc[kp], f$spat$yloc[kp],
                              window = f$win, check = FALSE)
      spatstat.geom::marks(Y) <- factor(ifelse(pa[kp], "i", "j"), levels = c("i", "j"))
      ref <- as.data.frame(spatstat.explore::Gcross(
        Y, "i", "j", r = f$r, correction = ec))[[GCOL[[ec]]]]
      mine <- got[got$Anchor == a & got$Counted == c2, ]
      mine <- mine[order(mine$r), ]
      expect_equal(mine$`Observed G`, ref, tolerance = 1e-12, info = paste(ec, a))
    }
  }
})

test_that("bivariate km correction returns km, not rs", {
  # Regression guard for selecting spatstat's estimate by position: for
  # correction = "km", Gcross's third column is `rs`.
  f <- g_fixture()
  pa <- f$spat[[f$mnames[1]]] == 1 & f$spat[[f$mnames[2]]] != 1
  pc <- f$spat[[f$mnames[2]]] == 1 & f$spat[[f$mnames[1]]] != 1
  kp <- pa | pc
  Y <- spatstat.geom::ppp(f$spat$xloc[kp], f$spat$yloc[kp], window = f$win, check = FALSE)
  spatstat.geom::marks(Y) <- factor(ifelse(pa[kp], "i", "j"), levels = c("i", "j"))
  gg <- as.data.frame(spatstat.explore::Gcross(Y, "i", "j", r = f$r, correction = "km"))
  expect_identical(names(gg)[3], "rs")   # the trap this test guards

  got <- bi_NN_G(f$mif, mnames = f$mnames, r_range = f$r, num_permutations = 1,
                 workers = 1, edge_correction = "km",
                 overwrite = TRUE)$derived$bivariate_NN
  mine <- got[got$Anchor == f$mnames[1] & got$Counted == f$mnames[2], ]
  mine <- mine[order(mine$r), ]
  expect_equal(mine$`Observed G`, gg$km, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(mine$`Observed G`, gg$rs)))
})

test_that("G output schema matches Ripley's K except for the Observed column", {
  # Requirement: NN G output columns should line up with Ripley's K so results from
  # different metrics can be handled by the same downstream code.
  f <- g_fixture()
  kc <- names(ripleys_k(f$mif, mnames = f$mnames, r_range = f$r, permute = FALSE,
                        workers = 1, overwrite = TRUE)$derived$univariate_Count)
  gc <- names(NN_G(f$mif, mnames = f$mnames, r_range = f$r, num_permutations = 1,
                   workers = 1, overwrite = TRUE)$derived$univariate_NN)
  expect_identical(setdiff(kc, gc), "Observed K")
  expect_identical(setdiff(gc, kc), "Observed G")
  # Same for the bivariate pair.
  bkc <- names(bi_ripleys_k(f$mif, mnames = f$mnames, r_range = f$r, permute = FALSE,
                            workers = 1, overwrite = TRUE)$derived$bivariate_Count)
  bgc <- names(bi_NN_G(f$mif, mnames = f$mnames, r_range = f$r, num_permutations = 1,
                       workers = 1, overwrite = TRUE)$derived$bivariate_NN)
  expect_identical(setdiff(bkc, bgc), "Observed K")
  expect_identical(setdiff(bgc, bkc), "Observed G")
})

test_that("Exact CSR is present but always NA for G", {
  # G is intensity-dependent, so unlike K there is no closed form for the expected
  # value under random labelling. The column exists only for schema consistency.
  f <- g_fixture()
  u <- NN_G(f$mif, mnames = f$mnames, r_range = f$r, num_permutations = 2,
            workers = 1, overwrite = TRUE)$derived$univariate_NN
  b <- bi_NN_G(f$mif, mnames = f$mnames, r_range = f$r, num_permutations = 2,
               workers = 1, overwrite = TRUE)$derived$bivariate_NN
  expect_true("Exact CSR" %in% names(u) && all(is.na(u$`Exact CSR`)))
  expect_true("Exact CSR" %in% names(b) && all(is.na(b$`Exact CSR`)))

  # And demonstrate WHY: Gest over all cells is nothing like the permuted mean.
  pp_all <- spatstat.geom::ppp(f$spat$xloc, f$spat$yloc, window = f$win, check = FALSE)
  g_all <- as.data.frame(spatstat.explore::Gest(pp_all, r = f$r, correction = "rs"))$rs
  perm <- u[u$Marker == f$mnames[1], ]
  perm <- perm[order(perm$r), ]
  mid <- which(f$r == 30)
  expect_gt(abs(g_all[mid] - perm$`Permuted CSR`[mid]), 0.2)
})

test_that("G permutations are reproducible and independent of workers", {
  f <- g_fixture()
  runp <- function(fn, w) {
    set.seed(7)
    d <- fn(f$mif, mnames = f$mnames, r_range = f$r, num_permutations = 4,
            keep_permutation_distribution = TRUE, workers = w,
            overwrite = TRUE)$derived
    d[[grep("NN", names(d))]]$`Permuted CSR`
  }
  expect_equal(runp(NN_G, 1), runp(NN_G, 1))
  expect_equal(runp(NN_G, 1), runp(NN_G, 3))
  expect_equal(runp(bi_NN_G, 1), runp(bi_NN_G, 3))
})

test_that("G edge-correction spellings are normalised and bad ones rejected", {
  expect_identical(match_g_correction("rs"), "rs")
  expect_identical(match_g_correction("km"), "km")
  expect_identical(match_g_correction("han"), "han")
  # "hans" was propagated through bi_NN_G_sample() and the compute_metrics() docs
  # before 2.0.0 and was handed straight to spatstat, which only knows "han".
  expect_identical(match_g_correction("hans"), "han")
  expect_identical(match_g_correction("Hanisch"), "han")
  expect_identical(match_g_correction("none"), "none")
  expect_error(match_g_correction("bogus"), "Unsupported")
  expect_error(match_g_correction(c("rs", "km")), "single string")
})

test_that("bi_NN_G rejects a single marker and NN_G accepts one", {
  f <- g_fixture()
  expect_error(bi_NN_G(f$mif, mnames = f$mnames[1], r_range = f$r), "univariate")
  expect_s3_class(
    NN_G(f$mif, mnames = f$mnames[1], r_range = f$r, num_permutations = 1,
         workers = 1, overwrite = TRUE),
    "mif"
  )
})

test_that("sparse markers yield NA rather than NaN", {
  # TMA1_[3,B] has a single PDL1-positive cell. The pre-2.0.0 hand-rolled path
  # produced NaN here; the guard now returns a clean NA.
  mif <- create_mif(
    clinical_data = example_clinical %>%
      dplyr::mutate(deidentified_id = as.character(deidentified_id)),
    sample_data = example_summary %>%
      dplyr::mutate(deidentified_id = as.character(deidentified_id)),
    spatial_list = example_spatial["TMA1_[3,B].tif"],
    patient_id = "deidentified_id", sample_id = "deidentified_sample")
  out <- bi_NN_G(mif, mnames = c("CD3..Opal.570..Positive", "PDL1..Opal.540..Positive"),
                 r_range = seq(0, 60, 10), num_permutations = 2, workers = 1,
                 overwrite = TRUE)$derived$bivariate_NN
  expect_true(all(is.na(out$`Observed G`)))
  expect_false(any(is.nan(out$`Observed G`)))
})
