# The observation window for a sample must be the convex hull of EVERY cell in
# that sample's spatial data frame -- never of a marker-positive subset, an
# anchor/counted subset, or a dual-positive-filtered subset -- and it must stay
# fixed across all markers, all marker pairs and all permutations within a sample.
#
# This is easy to break. The natural-looking "build a ppp per marker" refactor gets
# it wrong, and the failure is silent: you still get plausible K values, just
# inflated by the ratio of the full window's area to the subset's. The fixture
# below makes that ratio ~8x so the failure is unmistakable.
#
# Why it matters statistically: K is normalised by the observation area. Shrinking
# the window to the marker's own bounding hull both shrinks the area and removes
# the edge correction appropriate to the real tissue boundary, so markers with
# different spatial extents stop being comparable to each other and to the
# Exact/Permuted CSR references computed from all cells.

skip_if_not_installed("spatstat.geom")
skip_if_not_installed("spatstat.explore")

# A sample whose cells fill a large square, but whose markers are each confined to
# one corner of it.
corner_marker_mif <- function(seed = 6, n = 3000) {
  set.seed(seed)
  spat <- data.frame(
    deidentified_sample = "S1",
    XMin = runif(n, 0, 1000),
    YMin = runif(n, 0, 1000)
  )
  spat$XMax <- spat$XMin
  spat$YMax <- spat$YMin
  spat$cornerA <- as.integer(spat$XMin < 350 & spat$YMin < 350)
  spat$cornerB <- as.integer(spat$XMin > 650 & spat$YMin > 650)
  # A third marker spread over the whole sample, to contrast with the corner ones.
  spat$spread <- as.integer(seq_len(n) %% 7 == 0)
  list(
    spat = spat,
    mif = create_mif(
      clinical_data = data.frame(deidentified_id = "p1"),
      sample_data   = data.frame(deidentified_id = "p1", deidentified_sample = "S1"),
      spatial_list  = list(S1 = spat),
      patient_id = "deidentified_id",
      sample_id  = "deidentified_sample"
    )
  )
}

test_that("the corner fixture really does have a big area discrepancy", {
  # Guard the guard: if this ever stops holding, the tests below stop testing
  # anything, because full-window and subset-window answers would coincide.
  f <- corner_marker_mif()
  W_all <- spatstat.geom::convexhull.xy(f$spat$XMin, f$spat$YMin)
  pos <- f$spat$cornerA == 1
  W_sub <- spatstat.geom::convexhull.xy(f$spat$XMin[pos], f$spat$YMin[pos])
  expect_gt(spatstat.geom::area(W_all) / spatstat.geom::area(W_sub), 4)
  expect_gt(sum(pos), 50)
})

test_that("ripleys_k uses the full-sample window, not the marker subset's", {
  f <- corner_marker_mif()
  r <- seq(0, 80, 10)
  got <- ripleys_k(f$mif, mnames = "cornerA", r_range = r, permute = FALSE,
                   workers = 1, edge_correction = "translation",
                   overwrite = TRUE)$derived$univariate_Count
  x <- f$spat$XMin; y <- f$spat$YMin; pos <- f$spat$cornerA == 1
  ref_full <- as.data.frame(spatstat.explore::Kest(
    spatstat.geom::ppp(x[pos], y[pos],
                       window = spatstat.geom::convexhull.xy(x, y), check = FALSE),
    r = r, correction = "translation"))$trans
  ref_sub <- as.data.frame(spatstat.explore::Kest(
    spatstat.geom::ppp(x[pos], y[pos],
                       window = spatstat.geom::convexhull.xy(x[pos], y[pos]), check = FALSE),
    r = r, correction = "translation"))$trans

  expect_equal(got$`Observed K`, ref_full, tolerance = 1e-12)
  # And is emphatically NOT the subset-window answer.
  expect_gt(max(abs(got$`Observed K` - ref_sub), na.rm = TRUE), 1)
})

test_that("bi_ripleys_k uses the full-sample window for disjoint corner markers", {
  f <- corner_marker_mif()
  r <- seq(0, 80, 10)
  got <- bi_ripleys_k(f$mif, mnames = c("cornerA", "cornerB"), r_range = r,
                      permute = FALSE, workers = 1, edge_correction = "translation",
                      overwrite = TRUE)$derived$bivariate_Count
  x <- f$spat$XMin; y <- f$spat$YMin
  pa <- f$spat$cornerA == 1; pb <- f$spat$cornerB == 1
  keep <- pa | pb
  Y <- spatstat.geom::ppp(x[keep], y[keep],
                          window = spatstat.geom::convexhull.xy(x, y), check = FALSE)
  spatstat.geom::marks(Y) <- factor(ifelse(pa[keep], "a", "b"), levels = c("a", "b"))
  ref_full <- as.data.frame(spatstat.explore::Kcross(Y, "a", "b", r = r,
                                                    correction = "translation"))$trans
  mine <- got[got$Anchor == "cornerA" & got$Counted == "cornerB", ]
  mine <- mine[order(mine$r), ]
  expect_equal(mine$`Observed K`, ref_full, tolerance = 1e-12)
})

test_that("NN_G uses the full-sample window, not the marker subset's", {
  f <- corner_marker_mif()
  r <- seq(0, 80, 10)
  got <- NN_G(f$mif, mnames = "cornerA", r_range = r, num_permutations = 2,
              workers = 1, edge_correction = "rs",
              overwrite = TRUE)$derived$univariate_NN
  x <- f$spat$XMin; y <- f$spat$YMin; pos <- f$spat$cornerA == 1
  gk <- function(W) as.data.frame(spatstat.explore::Gest(
    spatstat.geom::ppp(x[pos], y[pos], window = W, check = FALSE),
    r = r, correction = "rs"))$rs
  expect_equal(got$`Observed G`, gk(spatstat.geom::convexhull.xy(x, y)), tolerance = 1e-12)
  # G is bounded in [0,1] so the discrepancy cannot be huge, but the edge
  # correction differs enough to be well outside tolerance.
  expect_gt(max(abs(got$`Observed G` -
                      gk(spatstat.geom::convexhull.xy(x[pos], y[pos]))), na.rm = TRUE), 1e-6)
})

test_that("bi_NN_G uses the full-sample window for disjoint corner markers", {
  f <- corner_marker_mif()
  r <- seq(0, 200, 25)   # corners are ~450 units apart, so use a wider range
  got <- bi_NN_G(f$mif, mnames = c("cornerA", "cornerB"), r_range = r,
                 num_permutations = 2, workers = 1, edge_correction = "rs",
                 overwrite = TRUE)$derived$bivariate_NN
  x <- f$spat$XMin; y <- f$spat$YMin
  pa <- f$spat$cornerA == 1; pb <- f$spat$cornerB == 1; keep <- pa | pb
  Y <- spatstat.geom::ppp(x[keep], y[keep],
                          window = spatstat.geom::convexhull.xy(x, y), check = FALSE)
  spatstat.geom::marks(Y) <- factor(ifelse(pa[keep], "i", "j"), levels = c("i", "j"))
  ref <- as.data.frame(spatstat.explore::Gcross(Y, "i", "j", r = r,
                                               correction = "rs"))$rs
  mine <- got[got$Anchor == "cornerA" & got$Counted == "cornerB", ]
  mine <- mine[order(mine$r), ]
  expect_equal(mine$`Observed G`, ref, tolerance = 1e-12)
})

test_that("a marker's result does not depend on which other markers were requested", {
  # If the window or area were derived from the requested markers rather than from
  # all cells, asking for more markers would move every marker's numbers.
  f <- corner_marker_mif()
  r <- seq(0, 60, 10)
  one <- ripleys_k(f$mif, mnames = "cornerA", r_range = r, permute = FALSE,
                   workers = 1, overwrite = TRUE)$derived$univariate_Count
  many <- ripleys_k(f$mif, mnames = c("cornerA", "cornerB", "spread"), r_range = r,
                    permute = FALSE, workers = 1, overwrite = TRUE)$derived$univariate_Count
  many_a <- many[many$Marker == "cornerA", ]
  many_a <- many_a[order(many_a$r), ]
  expect_equal(one$`Observed K`, many_a$`Observed K`)
  # Exact CSR is a property of the sample, so it must be identical too.
  expect_equal(one$`Exact CSR`, many_a$`Exact CSR`)
})

test_that("Observed, Exact and Permuted CSR all share one window within a sample", {
  # Exact CSR is the K of all cells; the mean of the permuted distribution
  # estimates the same quantity. They can only agree if both use the same window.
  f <- corner_marker_mif()
  r <- seq(0, 60, 10)
  set.seed(31)
  perm <- ripleys_k(f$mif, mnames = "spread", r_range = r, permute = TRUE,
                    num_permutations = 200, keep_permutation_distribution = FALSE,
                    workers = 1, overwrite = TRUE)$derived$univariate_Count
  exact <- ripleys_k(f$mif, mnames = "spread", r_range = r, permute = FALSE,
                     workers = 1, overwrite = TRUE)$derived$univariate_Count
  ok <- exact$`Exact CSR` > 0
  rel <- abs(perm$`Permuted CSR`[ok] - exact$`Exact CSR`[ok]) / exact$`Exact CSR`[ok]
  # Monte Carlo error over 200 permutations, not a systematic window difference.
  expect_lt(max(rel), 0.10)
})
