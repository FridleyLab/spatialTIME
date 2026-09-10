# Holds the refactor to the numbers v1.4.0 produced, for the quantities where that
# is a meaningful thing to ask.
#
# Only the columns named in the fixture's `deterministic_cols` attribute are
# compared exactly. The v1.4.0 permutation paths are NOT reproducible -- 24 nested
# mclapply() calls omitted mc.cores and forked with their own RNG streams, so
# set.seed() had no effect and two runs at the same seed differed by thousands.
# Comparing `Permuted CSR` against that baseline would be comparing against noise.
#
# Two exact-column differences ARE expected and are asserted as differences, so
# that reverting either fix would fail this file:
#
#   edge_correction = "none"  -- now uses spatstat's whist binning rather than
#                                Kest's fast-path binning, which disagree on tied
#                                distances. The example data has half-integer cell
#                                centres, so ties are common.
#   samples over `big`        -- now keep the requested edge correction instead of
#                                being silently downgraded to "none".

baseline <- function() {
  path <- test_path("fixtures", "baseline-v1.4.0.rds")
  skip_if_not(file.exists(path), "v1.4.0 baseline fixture not present")
  readRDS(path)
}

baseline_mif <- function(b) {
  create_mif(
    clinical_data = example_clinical %>%
      dplyr::mutate(deidentified_id = as.character(deidentified_id)),
    sample_data = example_summary %>%
      dplyr::mutate(deidentified_id = as.character(deidentified_id)),
    spatial_list = example_spatial[attr(b, "sample")],
    patient_id = "deidentified_id", sample_id = "deidentified_sample")
}

test_that("the baseline fixture is what we think it is", {
  b <- baseline()
  expect_identical(attr(b, "spatialTIME_version"), "1.4.0")
  expect_identical(attr(b, "sample"), "TMA3_[9,K].tif")
  expect_true(length(attr(b, "deterministic_cols")) > 0)
  expect_match(attr(b, "nondeterministic_note"), "irreproducible")
})

test_that("Observed K matches v1.4.0 for translation and isotropic", {
  b <- baseline()
  mk <- attr(b, "markers"); rr <- attr(b, "r_range")
  for (ec in c("translation", "isotropic")) {
    key <- c(translation = "ripleys_k.exact", isotropic = "ripleys_k.iso")[[ec]]
    old <- b[[key]]$univariate_Count
    old <- old[order(old$Marker, old$r), ]
    new <- ripleys_k(baseline_mif(b), mnames = mk, r_range = rr, permute = FALSE,
                     workers = 1, edge_correction = ec,
                     overwrite = TRUE)$derived$univariate_Count
    new <- new[order(new$Marker, new$r), ]
    # K here is ~1e4, so 1e-9 absolute is ~1e-13 relative: floating point only.
    expect_equal(new$`Observed K`, old$`Observed K`, tolerance = 1e-9, info = ec)
    expect_equal(new$`Theoretical CSR`, old$`Theoretical CSR`, info = ec)
    expect_equal(new$`Exact CSR`, old$`Exact CSR`, tolerance = 1e-9, info = ec)
  }
})

test_that("edge_correction = 'none' DIFFERS from v1.4.0, by design", {
  # v1.4.0 asked Kest for "none" alone, which takes a fast C path binning d <= r.
  # The engine always uses the whist path, so its "none" is bit-identical to
  # Kest(correction = c("none", "translation")) and consistent with its own
  # translation output. Only visible on tied distances -- which the half-integer
  # centres of this data produce constantly.
  b <- baseline()
  mk <- attr(b, "markers"); rr <- attr(b, "r_range")
  old <- b$ripleys_k.none$univariate_Count
  old <- old[order(old$Marker, old$r), ]
  new <- ripleys_k(baseline_mif(b), mnames = mk, r_range = rr, permute = FALSE,
                   workers = 1, edge_correction = "none",
                   overwrite = TRUE)$derived$univariate_Count
  new <- new[order(new$Marker, new$r), ]
  expect_false(isTRUE(all.equal(new$`Observed K`, old$`Observed K`)))

  # And the new value is the one spatstat gives through its whist path.
  spat <- example_spatial[[attr(b, "sample")]]
  xy <- cbind((spat$XMax + spat$XMin) / 2, (spat$YMax + spat$YMin) / 2)
  win <- spatstat.geom::convexhull.xy(xy[, 1], xy[, 2])
  pos <- spat[[mk[1]]] == 1
  ref <- as.data.frame(spatstat.explore::Kest(
    spatstat.geom::ppp(xy[pos, 1], xy[pos, 2], window = win, check = FALSE),
    r = rr, correction = c("none", "translation")))$un
  got <- new[new$Marker == mk[1], ]
  expect_equal(got$`Observed K`[order(got$r)], ref, tolerance = 1e-9)
})

test_that("bivariate NN G Observed G is unchanged from the hand-rolled version", {
  # bi_NN_G's rs and han estimators were hand-rolled on a full as.matrix(dist(...))
  # before 2.0.0; replacing them with Gcross() was expected to be numerically free,
  # and this is the assertion that says so.
  b <- baseline()
  mk <- attr(b, "markers"); rr <- attr(b, "r_range")
  for (ec in c("rs", "han")) {
    key <- c(rs = "bi_NN_G.rs", han = "bi_NN_G.han")[[ec]]
    old <- b[[key]]$bivariate_NN
    old <- old[order(old$Anchor, old$Counted, old$r), ]
    new <- bi_NN_G(baseline_mif(b), mnames = mk[1:2], r_range = rr,
                   num_permutations = 2, workers = 1, edge_correction = ec,
                   overwrite = TRUE)$derived$bivariate_NN
    new <- new[order(new$Anchor, new$Counted, new$r), ]
    expect_equal(new$`Observed G`, old$`Observed G`, info = ec)
  }
})

test_that("univariate NN G Observed G is unchanged for rs", {
  b <- baseline()
  mk <- attr(b, "markers"); rr <- attr(b, "r_range")
  old <- b$NN_G.rs$univariate_NN
  old <- old[order(old$Marker, old$r), ]
  new <- NN_G(baseline_mif(b), mnames = mk, r_range = rr, num_permutations = 2,
              workers = 1, edge_correction = "rs",
              overwrite = TRUE)$derived$univariate_NN
  new <- new[order(new$Marker, new$r), ]
  expect_equal(new$`Observed G`, old$`Observed G`)
})

test_that("NN_G(edge_correction = 'km') output was MALFORMED in v1.4.0", {
  # Not a numeric regression -- a structural one, and the baseline preserves the
  # evidence. v1.4.0 reordered columns POSITIONALLY (res[,c(7,6,4,1,2,5,3)]), which
  # assumed Gest returns 3 columns. For correction = "km" it returns five
  # (r, theo, km, hazard, theohaz), so the indices pointed at the wrong columns: the
  # captured baseline has NO sample-id column, NO Marker column, a leaked
  # `theohaz.x`, and 95 rows where there should be 21.
  b <- baseline()
  broken <- b$NN_G.km$univariate_NN
  expect_false("Marker" %in% names(broken))
  expect_false("deidentified_sample" %in% names(broken))
  expect_true("theohaz.x" %in% names(broken))

  # Selecting by name fixes it: correct shape, and equal to spatstat's km estimate.
  mk <- attr(b, "markers"); rr <- attr(b, "r_range")
  new <- NN_G(baseline_mif(b), mnames = mk, r_range = rr, num_permutations = 2,
              workers = 1, edge_correction = "km",
              overwrite = TRUE)$derived$univariate_NN
  expect_true(all(c("Marker", "deidentified_sample") %in% names(new)))
  expect_false(any(grepl("theohaz|hazard", names(new))))
  expect_equal(nrow(new), length(mk) * length(rr))

  spat <- example_spatial[[attr(b, "sample")]]
  xy <- cbind((spat$XMax + spat$XMin) / 2, (spat$YMax + spat$YMin) / 2)
  win <- spatstat.geom::convexhull.xy(xy[, 1], xy[, 2])
  pos <- spat[[mk[1]]] == 1
  ref <- as.data.frame(spatstat.explore::Gest(
    spatstat.geom::ppp(xy[pos, 1], xy[pos, 2], window = win, check = FALSE),
    r = rr, correction = "km"))$km
  got <- new[new$Marker == mk[1], ]
  expect_equal(got$`Observed G`[order(got$r)], ref)
})

test_that("permutations are now reproducible, which v1.4.0 could not manage", {
  # The converse assertion. This is a capability the baseline cannot satisfy, so it
  # is stated forwards rather than as a comparison.
  b <- baseline()
  mk <- attr(b, "markers")[1:2]; rr <- attr(b, "r_range")
  run <- function(w) {
    set.seed(42)
    ripleys_k(baseline_mif(b), mnames = mk, r_range = rr, permute = TRUE,
              num_permutations = 3, keep_permutation_distribution = TRUE,
              workers = w, overwrite = TRUE)$derived$univariate_Count$`Permuted CSR`
  }
  expect_identical(run(1), run(1))
  # And independent of how many workers were used, which is the part that was
  # structurally impossible before.
  expect_identical(run(1), run(3))
})
