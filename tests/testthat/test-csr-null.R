# Is the CSR reference right?
#
# `Degree of Clustering X = Observed - X CSR`, so trusting that number needs two
# separate things to hold:
#
#   (i)  Observed is exactly what spatstat computes  -- test-k-engine.R
#   (ii) the CSR reference is unbiased for the null  -- this file
#
# (ii) had no coverage at all. test-k-engine.R pins the estimator, and the
# regression fixture pins the numbers to v1.4.0, but nothing asked whether the
# thing being subtracted is the right thing to subtract. A CSR reference that was
# systematically 10% low would pass every other test in the suite while making
# every degree of clustering in every published table wrong.
#
# The load-bearing fact is that the identity below is EXACT, not asymptotic, so it
# can be tested by enumeration with no RNG and no guessed tolerance:
#
#   E[K_permuted(r)] = K of ALL cells = `Exact CSR`
#
# Why: `sample.int(n, m)` samples without replacement, so an ordered pair (i,j) is
# retained with probability m(m-1)/(n(n-1)); the denominator is m(m-1)/area; the
# two cancel, leaving area * sum(w * 1{d<=r}) / (n(n-1)), which is precisely
# `k_from_pairs(pairs, rep(TRUE, n))`. Holds for every r, every m >= 2, and every
# edge correction whose denominator is fixed once m is -- translation, isotropic
# and none. NOT border, which is why border is excluded here and why `Exact CSR`
# is NA for it; see the test at the bottom.
#
# The same cancellation with P(i in anchor and j in counted) = n_i*n_j/(n(n-1))
# against a denominator of n_i*n_j/area is what makes `bi_ripleys_k`'s `Exact CSR`
# the UNIVARIATE K of all cells. That is the only test anywhere of that claim.

skip_if_not_installed("spatstat.geom")
skip_if_not_installed("spatstat.explore")

# Corrections for which the unbiasedness identity holds exactly. Border is absent
# deliberately -- see the file header and the final test.
EXACT_CSR_CORR <- c("translation", "isotropic", "none")

small_pattern <- function(n = 10, seed = 4) {
  set.seed(seed)
  x <- stats::runif(n, 0, 100); y <- stats::runif(n, 0, 100)
  W <- spatstat.geom::convexhull.xy(x, y)
  spatstat.geom::ppp(x, y, window = W, check = FALSE)
}

# --- (1) exact unbiasedness, by enumeration -----------------------------------

test_that("Exact CSR is the exact mean of K over EVERY random labelling", {
  # Complete enumeration of all C(10,3) = 120 labellings. No RNG, so there is no
  # Monte Carlo error to tolerate: the tolerance below is floating-point only.
  pp <- small_pattern()
  n  <- spatstat.geom::npoints(pp)
  m  <- 3L
  combs <- utils::combn(n, m)

  # r deliberately runs past where spatstat stops reporting, so the identity is
  # checked across the NA boundary too. That is the cheapest detector of the
  # k_from_pairs() early-return bug, which reported 0 instead of NA whenever a
  # marker's selected cells had no pair within rmax.
  #
  # Each correction truncates at a different radius -- on this fixture translation
  # stops at 72.7 and isotropic at 58.0 -- and "none" has rmax_valid = Inf and so
  # never truncates at all. So the NA assertion below is keyed off the correction's
  # own rmax_valid rather than assuming one radius covers all three.
  r <- c(0, 10, 20, 30, 80)

  for (ec in EXACT_CSR_CORR) {
    pairs <- k_pairs(pp, r, ec)
    mean_over_labellings <- rowMeans(apply(combs, 2, function(idx) {
      keep <- logical(n); keep[idx] <- TRUE
      k_from_pairs(pairs, keep)
    }))
    expect_equal(mean_over_labellings, k_from_pairs(pairs, rep(TRUE, n)),
                 tolerance = 1e-12, info = ec)

    truncates <- is.finite(pairs$rmax_valid) && max(r) >= pairs$rmax_valid
    if (truncates) {
      # The NA region really was exercised, rather than the test passing on a
      # vector that happened to be finite everywhere.
      expect_true(anyNA(mean_over_labellings), info = ec)
      expect_true(all(is.na(mean_over_labellings[r >= pairs$rmax_valid])), info = ec)
    } else {
      expect_false(anyNA(mean_over_labellings), info = ec)
    }
  }

  # At least one correction must have exercised the NA path, or the paragraph above
  # is decoration.
  expect_true(any(vapply(EXACT_CSR_CORR, function(ec) {
    is.finite(k_pairs(pp, r, ec)$rmax_valid)
  }, logical(1))))
})

test_that("bivariate Exact CSR is the exact mean of Kcross over every labelling", {
  # The claim under test is that bi_ripleys_k's `Exact CSR` -- the UNIVARIATE K of
  # all cells -- is E[Kcross] under random labelling. Enumerating all
  # C(8,2) * C(6,3) = 560 anchor/counted configurations settles it.
  pp <- small_pattern(n = 8, seed = 12)
  n  <- spatstat.geom::npoints(pp)
  n_i <- 2L; n_j <- 3L
  r <- c(0, 10, 20, 30, 60)

  anchor_sets <- utils::combn(n, n_i, simplify = FALSE)
  configs <- unlist(lapply(anchor_sets, function(a) {
    lapply(utils::combn(setdiff(seq_len(n), a), n_j, simplify = FALSE),
           function(b) list(a = a, b = b))
  }), recursive = FALSE)
  expect_equal(length(configs), choose(8, 2) * choose(6, 3))

  for (ec in EXACT_CSR_CORR) {
    pairs <- k_pairs(pp, r, ec)
    mean_over_configs <- rowMeans(vapply(configs, function(cf) {
      ki <- logical(n); ki[cf$a] <- TRUE
      kj <- logical(n); kj[cf$b] <- TRUE
      k_from_pairs(pairs, ki, kj, n_i, n_j, univariate = FALSE)
    }, numeric(length(r))))
    expect_equal(mean_over_configs, k_from_pairs(pairs, rep(TRUE, n)),
                 tolerance = 1e-12, info = ec)
  }
})

# --- (2) observed and CSR must be the SAME estimator --------------------------

test_that("Observed K of an all-positive marker IS Exact CSR", {
  # Both are k_from_pairs(pairs, rep(TRUE, n)). If the observed path ever forks
  # away from the CSR path -- which is exactly what v1.4.0 did, with a separate
  # hand-rolled estimator on the permute branch -- these stop being equal.
  f <- dispersion_mif()
  spat <- f$spat
  spat$all_cells <- 1L
  m <- toy_mif(list(S1 = spat))

  for (ec in EXACT_CSR_CORR) {
    out <- ripleys_k(m, mnames = "all_cells", r_range = c(0, 20, 40, 60),
                     edge_correction = ec, permute = FALSE, workers = 1,
                     overwrite = TRUE)$derived$univariate_Count
    expect_equal(out$`Observed K`, out$`Exact CSR`, tolerance = 1e-12, info = ec)
  }
})

# --- (3) the permutation distribution must not be degenerate -----------------

test_that("every permutation is a distinct draw", {
  # A null where all permutations return the same value has the CORRECT mean and
  # would pass tests 1, 2 and 5. This is the only thing here that catches a
  # set.seed() moved inside the loop, or a reused permutation.
  f <- dispersion_mif()
  B <- 25L
  set.seed(3)
  out <- ripleys_k(f$mif, mnames = "labelled", r_range = c(0, 20, 40),
                   permute = TRUE, num_permutations = B,
                   keep_permutation_distribution = TRUE, workers = 1,
                   overwrite = TRUE)$derived$univariate_Count

  # Precondition, asserted rather than assumed: at small radii many relabellings
  # genuinely produce the identical statistic (at r = 0 every one gives 0), so
  # distinctness is only meaningful at a radius where the statistic can vary.
  # Measured: 20 of 300 distinct at r = 5, all 300 distinct at r >= 20.
  at_zero <- out$`Permuted CSR`[out$r == 0]
  expect_length(unique(at_zero), 1L)

  for (rr in c(20, 40)) {
    v <- out$`Permuted CSR`[out$r == rr]
    expect_length(v, B)
    expect_length(unique(v), B)
  }
})

# --- (4) two-sided: clustering AND inhibition --------------------------------

test_that("degree of clustering has the right sign for clustered and inhibited markers", {
  # Every other fixture in the suite is random-labelled or clustered, so a measure
  # that reported the right magnitude with a flipped sign would pass all of them.
  f <- dispersion_mif()
  expect_gt(f$n_kept, 50)   # the hard core must still leave a usable marker

  set.seed(5)
  out <- ripleys_k(f$mif, mnames = c("clustered", "inhibited", "labelled"),
                   r_range = c(0, 20, 40, 60), permute = TRUE,
                   num_permutations = 200, workers = 1,
                   overwrite = TRUE)$derived$univariate_Count

  # Inequalities, not == 0 / == B: the exact count depends on the fixture's seed.
  # Restricted to r >= 20, because at r = 0 every K is 0 by construction.
  big_r <- out$r >= 20
  for (mk in c("clustered", "inhibited", "labelled")) {
    s <- out[big_r & out$Marker == mk, ]
    doc <- s$`Degree of Clustering Exact`
    p   <- s$`Permutation p-value`
    if (mk == "clustered") {
      expect_true(all(doc > 0), info = mk)
      expect_true(all(p < 0.05), info = mk)
    } else if (mk == "inhibited") {
      expect_true(all(doc < 0), info = mk)
      expect_true(all(p > 0.95), info = mk)
    } else {
      # Random labelling: the null is true by construction, so the degree of
      # clustering must be small RELATIVE to the reference rather than zero.
      expect_true(all(abs(doc) / s$`Exact CSR` < 0.25), info = mk)
    }
  }
})

# --- (5) Monte Carlo convergence, with a derived threshold -------------------

test_that("the permuted mean converges on Exact CSR at the Monte Carlo rate", {
  # Because the identity in test 1 is exact, the only discrepancy left is Monte
  # Carlo error -- so the tolerance can be DERIVED from the permutations instead of
  # guessed. z = sqrt(B) * (mean_B - exact) / sd_B is N(0,1) under the CLT.
  #
  # |z| < 5 is a per-radius two-sided failure probability of 5.7e-7. The radii are
  # strongly positively correlated (K is cumulative, and every radius shares the
  # same permutations), so the effective number of independent tests is ~1-2 and
  # the family-wise flake stays below 1e-6.
  f <- dispersion_mif()
  n <- nrow(f$spat)
  keep <- f$spat$labelled == 1L
  pp <- spatstat.geom::ppp(f$spat$XMin, f$spat$YMin,
                           window = spatstat.geom::convexhull.xy(f$spat$XMin, f$spat$YMin),
                           check = FALSE)
  r <- c(0, 20, 40, 60, 80)
  B <- 200L

  for (ec in EXACT_CSR_CORR) {
    pairs <- k_pairs(pp, r, ec)
    exact <- k_from_pairs(pairs, rep(TRUE, n))
    set.seed(41)
    P <- vapply(seq_len(B), function(b) {
      kp <- logical(n); kp[sample.int(n, sum(keep))] <- TRUE
      k_from_pairs(pairs, kp)
    }, numeric(length(r)))

    mu <- rowMeans(P); sd_b <- apply(P, 1, stats::sd)
    # Drop radii where z is undefined: r = 0 has zero spread, and anything past
    # rmax_valid is NA. Without this the test would be comparing against NaN.
    usable <- sd_b > 0 & !is.na(exact) & !is.na(mu)
    expect_gt(sum(usable), 2L)

    z <- sqrt(B) * (mu[usable] - exact[usable]) / sd_b[usable]
    expect_lt(max(abs(z)), 5, label = sprintf("max|z| (%s)", ec))
  }
})

# --- (6) p-values must be calibrated under a true null ----------------------

test_that("permutation p-values are uniform when the null is true", {
  # Builds the marker BY random labelling, so the null holds by construction and
  # the p-values must be ~Uniform. This is the only end-to-end check that the
  # `Permutations Larger than Observed` / `Permutation p-value` bookkeeping is
  # right rather than merely self-consistent.
  #
  # A binned chi-square, not a rejection rate at alpha = 0.05 (which looks at one
  # tail and has almost no power here) and not KS (p is discrete on {0, 1/B, ..., 1},
  # which KS is not valid for).
  skip_on_cran()
  set.seed(77)
  n <- 500L; m <- 80L; B <- 100L; R <- 200L
  x <- stats::runif(n, 0, 500); y <- stats::runif(n, 0, 500)
  W  <- spatstat.geom::convexhull.xy(x, y)
  pp <- spatstat.geom::ppp(x, y, window = W, check = FALSE)
  r0 <- 50
  pairs <- k_pairs(pp, c(0, r0), "translation")

  draw <- function() { k <- logical(n); k[sample.int(n, m)] <- TRUE; k }
  ps <- vapply(seq_len(R), function(i) {
    obs <- k_from_pairs(pairs, draw())[2]
    perm <- vapply(seq_len(B), function(b) k_from_pairs(pairs, draw())[2], numeric(1))
    (1 + sum(perm >= obs)) / (B + 1)
  }, numeric(1))

  # Precondition: a heavily tied permutation distribution compresses p toward the
  # extremes and would make this test vacuous. Assert the radius is non-degenerate.
  perm_check <- vapply(seq_len(B), function(b) k_from_pairs(pairs, draw())[2], numeric(1))
  expect_gt(length(unique(perm_check)) / B, 0.95)

  # Measured at this design: mean 0.497, 5-bin counts 48 37 31 42 42, chi-sq 4.05.
  expect_lt(abs(mean(ps) - 0.5), 0.085)          # 4 sigma; sd = sqrt(1/12)/sqrt(R)
  counts <- as.integer(table(cut(ps, breaks = seq(0, 1, 0.2), include.lowest = TRUE)))
  chisq <- sum((counts - R / 5)^2) / (R / 5)
  expect_lt(chisq, stats::qchisq(0.9999, df = 4))
})

test_that("the wrapper's p-value bookkeeping matches a hand count", {
  # Test 6 works at engine level. This one goes through ripleys_k() itself, so the
  # column actually shipped to users is the one being checked.
  f <- dispersion_mif()
  B <- 50L
  set.seed(8)
  out <- ripleys_k(f$mif, mnames = "labelled", r_range = c(0, 20, 40),
                   permute = TRUE, num_permutations = B,
                   keep_permutation_distribution = TRUE, workers = 1,
                   overwrite = TRUE)$derived$univariate_Count

  for (rr in c(0, 20, 40)) {
    s <- out[out$r == rr, ]
    obs <- unique(s$`Observed K`)
    expect_length(obs, 1L)
    hand_count <- sum(s$`Permuted CSR` >= obs)
    expect_equal(unique(s$`Permutations Larger than Observed`), hand_count)
    expect_equal(unique(s$`Permutation p-value`), (1 + hand_count) / (B + 1))
  }
})

test_that("a p-value is never exactly zero, and is NA where nothing was estimated", {
  f <- dispersion_mif()
  B <- 20L
  set.seed(9)
  # Isotropic truncates at the window's bounding radius, 341.7 on this fixture, so
  # r = 400 is past it and Observed K is NA there. Translation would need r > 493
  # and "none" never truncates, so the correction has to be named explicitly.
  out <- ripleys_k(f$mif, mnames = "clustered", r_range = c(0, 20, 400),
                   edge_correction = "isotropic", permute = TRUE,
                   num_permutations = B, workers = 1,
                   overwrite = TRUE)$derived$univariate_Count

  p <- out$`Permutation p-value`
  expect_true(all(p > 0, na.rm = TRUE))
  expect_true(all(p >= 1 / (B + 1), na.rm = TRUE))

  # The NA case: a count of 0 here would read as "no permutation exceeded the
  # observation", i.e. maximal clustering, at a radius where nothing was computed.
  na_row <- is.na(out$`Observed K`)
  expect_true(any(na_row))
  expect_true(all(is.na(out$`Permutations Larger than Observed`[na_row])))
  expect_true(all(is.na(out$`Permutation p-value`[na_row])))
})

# --- border: no exact CSR, and that is deliberate ---------------------------

test_that("border reports no Exact CSR, because it has none", {
  # Border's denominator counts the SELECTED cells still further than r from the
  # window edge, so it varies from one relabelling to the next and
  # E[numerator/denominator] != E[numerator]/E[denominator]. The K of all cells is
  # therefore NOT E[K_permuted] for border.
  #
  # Measured on 600 cells, m = 120, B = 3000: the K of all cells sits below the
  # mean permuted K by 0.3% at r = 20 rising to 1.0% at r = 80, z = -0.9 to -7.9.
  # Real bias, not noise, growing with r. A column named "Exact CSR" must not be
  # 1% wrong, so it is NA -- same as G and the pair correlation function.
  f <- dispersion_mif()
  out <- ripleys_k(f$mif, mnames = "labelled", r_range = c(0, 20, 40, 60),
                   edge_correction = "border", permute = FALSE, workers = 1,
                   overwrite = TRUE)$derived$univariate_Count
  expect_true(all(is.na(out$`Exact CSR`)))
  expect_true(all(is.na(out$`Degree of Clustering Exact`)))
  # Observed is still populated -- it is the estimator that is fine, not the reference.
  expect_false(anyNA(out$`Observed K`))

  # And demonstrate the bias the NA is avoiding, so a future change that
  # "helpfully" fills the column in gets a failure that explains itself.
  n <- nrow(f$spat)
  pp <- spatstat.geom::ppp(f$spat$XMin, f$spat$YMin,
                           window = spatstat.geom::convexhull.xy(f$spat$XMin, f$spat$YMin),
                           check = FALSE)
  r <- c(0, 40, 80)
  pairs <- k_pairs(pp, r, "border")
  set.seed(23)
  B <- 1500L
  P <- vapply(seq_len(B), function(b) {
    kp <- logical(n); kp[sample.int(n, 120)] <- TRUE
    k_from_pairs(pairs, kp)
  }, numeric(length(r)))
  mu <- rowMeans(P, na.rm = TRUE)
  k_all <- k_from_pairs(pairs, rep(TRUE, n))
  se <- apply(P, 1, stats::sd, na.rm = TRUE) / sqrt(B)

  # Both forms, because they fail for different reasons. The relative gap is the
  # stable quantity (0.0067-0.0084 measured across B = 400, 800, 1500); the z is
  # what says the gap is not Monte Carlo noise, but its denominator is itself
  # noisy, so it needs the larger B to sit clear of the threshold.
  expect_gt(abs(mu[3] - k_all[3]) / k_all[3], 0.004)
  expect_gt(abs(mu[3] - k_all[3]) / se[3], 2.5)
})
