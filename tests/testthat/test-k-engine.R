# The whole justification for the k_pairs()/k_from_pairs() engine is that it is
# EXACT -- equal to spatstat's own Kest/Kcross to floating-point precision -- while
# using memory proportional to the number of close pairs rather than n^2. If these
# tests ever fail, the engine has stopped being a drop-in for spatstat and the
# refactor's premise is broken. Do not relax a tolerance here without
# understanding which of the two documented spatstat quirks you are hitting
# (see the header comment of R/utils-k-engine.R).

skip_if_not_installed("spatstat.geom")
skip_if_not_installed("spatstat.explore")

# spatstat names its output column differently from its `correction` argument.
#
# EVERY correction the engine accepts must appear here. Before 2.0.0 "border" was
# missing from this vector, and that omission is the entire reason a border
# estimator that silently returned the UNCORRECTED numbers shipped: k_pairs() had
# no border branch, so "border" fell through with weight 1 and was bit-identical
# to "none", while match_edge_correction() accepted it and the docs advertised it.
# Nothing failed because nothing compared it to Kest. Do not shorten this vector.
CORR <- c(translation = "trans", isotropic = "iso", none = "un", border = "border")

# "none" and "border" are both in Kest's `fastcorrections`, so requesting either
# ALONE takes a fast C path that bins right-closed (d <= r). Pairing it with a
# non-fast correction routes it through whist() instead, which bins left-closed and
# is the convention this engine uses throughout. The two disagree whenever a pair
# distance lands exactly on a bin edge -- never on continuous coordinates, but
# constantly on the integer coordinates HALO and Vectra emit.
FAST_CORR <- c("none", "border")

ref_corrections <- function(ec) {
  if (ec %in% FAST_CORR) c(ec, "translation") else ec
}

kest_ref <- function(X, r, ec) {
  as.data.frame(spatstat.explore::Kest(
    X, r = r, correction = ref_corrections(ec)))[[CORR[[ec]]]]
}

kcross_ref <- function(Y, r, ec) {
  as.data.frame(spatstat.explore::Kcross(
    Y, "a", "b", r = r, correction = ref_corrections(ec)))[[CORR[[ec]]]]
}

# --- fixtures -----------------------------------------------------------------

continuous_case <- function() {
  set.seed(7)
  X0 <- spatstat.random::rpoispp(900)
  W  <- spatstat.geom::convexhull.xy(X0$x, X0$y)
  X  <- spatstat.geom::ppp(X0$x, X0$y, window = W)
  set.seed(8)
  lab <- sample(c("a", "b", "other"), spatstat.geom::npoints(X),
                replace = TRUE, prob = c(.25, .35, .40))
  list(X = X, lab = lab, r = seq(0, 0.09, length.out = 19))
}

# Integer coordinates are what HALO and Vectra actually emit, and they are the
# only place binning conventions are observable, so they get their own fixture.
integer_case <- function() {
  g <- expand.grid(x = seq(0, 400, 10), y = seq(0, 400, 10))
  W <- spatstat.geom::convexhull.xy(c(-5, 405, 405, -5), c(-5, -5, 405, 405))
  X <- spatstat.geom::ppp(g$x, g$y, window = W)
  set.seed(3)
  lab <- sample(c("a", "b", "other"), spatstat.geom::npoints(X),
                replace = TRUE, prob = c(.3, .3, .4))
  list(X = X, lab = lab, r = c(0, 5, 10, 15, 20, 25, 30))
}

marked <- function(X, lab) {
  keep <- lab %in% c("a", "b")
  Y <- X[keep]
  spatstat.geom::marks(Y) <- factor(lab[keep], levels = c("a", "b"))
  Y
}

# --- univariate ---------------------------------------------------------------

for (case_name in c("continuous", "integer")) {
  cs <- if (case_name == "continuous") continuous_case() else integer_case()
  n <- spatstat.geom::npoints(cs$X)

  for (ec in names(CORR)) {
    test_that(sprintf("univariate K equals Kest (%s, %s coords)", ec, case_name), {
      mine <- k_from_pairs(k_pairs(cs$X, cs$r, ec), rep(TRUE, n))
      expect_equal(mine, kest_ref(cs$X, cs$r, ec), tolerance = 1e-12)
    })

    test_that(sprintf("marker-subset K equals Kest on that subset (%s, %s coords)", ec, case_name), {
      # This is the window invariant in miniature: the pair list is built from ALL
      # cells, then masked. The reference subsets the pattern but keeps the same
      # window, so agreement proves the mask does not shrink the window.
      keep <- cs$lab == "a"
      mine <- k_from_pairs(k_pairs(cs$X, cs$r, ec), keep)
      expect_equal(mine, kest_ref(cs$X[keep], cs$r, ec), tolerance = 1e-12)
    })
  }
}

# --- bivariate ----------------------------------------------------------------

test_that("bivariate K equals Kcross for translation, none and border", {
  for (case_name in c("continuous", "integer")) {
    cs <- if (case_name == "continuous") continuous_case() else integer_case()
    Y <- marked(cs$X, cs$lab)
    # Isotropic is excluded here and tested separately below -- it has a genuine
    # closepairs/crosspairs degeneracy. Border has no such problem: its denominator
    # comes from boundary distances, which are identical either way.
    for (ec in c("translation", "none", "border")) {
      mine <- k_from_pairs(k_pairs(cs$X, cs$r, ec),
                           cs$lab == "a", cs$lab == "b", univariate = FALSE)
      expect_equal(mine, kcross_ref(Y, cs$r, ec), tolerance = 1e-12,
                   info = sprintf("%s / %s coords", ec, case_name))
    }
  }
})

test_that("bivariate isotropic K matches Kcross except at circle-through-vertex degeneracies", {
  # Kest uses closepairs(); Kmulti uses crosspairs(). Those compute the same pair
  # distance with up to 1 ulp of difference. Normally irrelevant -- but the
  # isotropic weight is 1/(arc fraction inside W), which jumps discontinuously when
  # the circle of radius d passes exactly through a vertex of W. Because W is the
  # convex hull of the cells, its vertices ARE cells, so such pairs do occur. When
  # one does, a 1-ulp distance difference lands on opposite sides of the jump.
  # Verified against a 2e6-point Monte Carlo integration of the true arc fraction:
  # the discrepancy is confined to the single degenerate pair.
  cs <- continuous_case()
  Y <- marked(cs$X, cs$lab)
  mine <- k_from_pairs(k_pairs(cs$X, cs$r, "isotropic"),
                       cs$lab == "a", cs$lab == "b", univariate = FALSE)
  ref  <- kcross_ref(Y, cs$r, "isotropic")
  expect_equal(mine, ref, tolerance = 1e-4)
  # All but a handful of radii must still agree to full precision.
  close_enough <- abs(mine - ref) < 1e-12
  expect_gte(sum(close_enough), length(ref) - 2L)

  # The integer-grid fixture has no such degeneracy, so it must be exact.
  ci <- integer_case()
  Yi <- marked(ci$X, ci$lab)
  expect_equal(
    k_from_pairs(k_pairs(ci$X, ci$r, "isotropic"),
                 ci$lab == "a", ci$lab == "b", univariate = FALSE),
    kcross_ref(Yi, ci$r, "isotropic"), tolerance = 1e-12
  )
})

# --- border is the one correction with an r-dependent denominator -------------

test_that("border is a reduced-sample estimator, not a reweighting of 'none'", {
  # The bug this guards: border with no implementation falls through k_pairs() with
  # weight 1, which makes it bit-identical to "none". Agreement with Kest$border
  # (above) already catches that, but assert the inequality directly so the failure
  # names the actual confusion rather than showing two long numeric vectors.
  cs <- continuous_case()
  n <- spatstat.geom::npoints(cs$X)
  bord <- k_from_pairs(k_pairs(cs$X, cs$r, "border"), rep(TRUE, n))
  none <- k_from_pairs(k_pairs(cs$X, cs$r, "none"), rep(TRUE, n))
  expect_gt(max(abs(bord - none), na.rm = TRUE), 1e-6)
})

test_that("border eligibility shrinks with r and empties past the inradius", {
  # Unlike every other correction, border's denominator counts only the anchors
  # whose distance to the window edge still exceeds r. Past the inradius that count
  # is zero and the estimator is genuinely 0/0 -- spatstat returns NaN there and so
  # must we, rather than silently reporting 0 or carrying the last finite value.
  set.seed(21)
  x <- runif(300, 0, 100); y <- runif(300, 0, 100)
  W <- spatstat.geom::convexhull.xy(x, y)
  X <- spatstat.geom::ppp(x, y, window = W, check = FALSE)
  n <- spatstat.geom::npoints(X)
  inradius <- max(spatstat.geom::bdist.points(X))
  r <- c(0, 10, 20, 40, inradius + 10, inradius + 30)

  mine <- k_from_pairs(k_pairs(X, r, "border"), rep(TRUE, n))
  ref  <- kest_ref(X, r, "border")
  expect_equal(mine, ref, tolerance = 1e-12)
  # Same NaN placement as spatstat, not merely the same finite values.
  expect_identical(is.nan(mine), is.nan(ref))
  expect_true(any(is.nan(mine)))
  expect_false(anyNA(mine[r < inradius]))
})

# --- memory bounding must not change the answer -------------------------------

test_that("chunked edge-weight computation is exact", {
  cs <- continuous_case()
  n <- spatstat.geom::npoints(cs$X)
  for (ec in names(CORR)) {
    expect_identical(
      k_from_pairs(k_pairs(cs$X, cs$r, ec, block = 500), rep(TRUE, n)),
      k_from_pairs(k_pairs(cs$X, cs$r, ec, block = Inf), rep(TRUE, n)),
      info = ec
    )
  }
})

test_that("ripleys_k's `big` changes memory only, never the numbers", {
  # `big` is documented as a memory knob: above that cell count the per-pair edge
  # weights are computed in chunks. It must not be observable in the output. This
  # matters because `big` USED to mean something else entirely -- in v1.4.0,
  # exceeding it silently replaced the requested edge correction with "none" -- so
  # an invariance test here is what keeps the new meaning honest.
  # bi_ripleys_k has had such a test since 2.0.0; ripleys_k did not.
  skip_if_not_installed("spatstat.explore")
  m <- example_mif_full()
  rr <- c(0, 10, 20, 30, 40, 50)
  run <- function(ec, b) {
    ripleys_k(m, mnames = mnames_good()[1], r_range = rr, permute = FALSE,
              edge_correction = ec, big = b, workers = 1,
              overwrite = TRUE)$derived$univariate_Count
  }
  for (ec in names(CORR)) {
    # 1803 cells in this fixture, so big = 500 genuinely chunks and big = 1e9 does not.
    expect_identical(run(ec, 500), run(ec, 1e9), info = ec)
  }
})

# --- weight reuse across permutations ----------------------------------------

test_that("reusing one pair list across permutations equals recomputing Kcross each time", {
  # This is the optimisation that removes the need for the old tiling code: edge
  # weights depend only on a pair's displacement and the window, never on which
  # other points are present, so they are computed once per sample and reused.
  cs <- continuous_case()
  X <- cs$X; n <- spatstat.geom::npoints(X); r <- cs$r
  pairs <- k_pairs(X, r, "translation")
  set.seed(11)
  for (i in 1:5) {
    s  <- sample(n, 400)
    ia <- logical(n); ia[s[1:150]]   <- TRUE
    ib <- logical(n); ib[s[151:400]] <- TRUE
    Y <- X[c(which(ia), which(ib))]
    spatstat.geom::marks(Y) <- factor(rep(c("a", "b"), c(150, 250)), levels = c("a", "b"))
    expect_equal(k_from_pairs(pairs, ia, ib, univariate = FALSE),
                 kcross_ref(Y, r, "translation"), tolerance = 1e-12)
  }
})

# --- guard rails --------------------------------------------------------------

test_that("edge_correction spellings are normalised and bad ones rejected loudly", {
  expect_identical(match_edge_correction("trans"), "translation")
  expect_identical(match_edge_correction("translation"), "translation")
  expect_identical(match_edge_correction("iso"), "isotropic")
  expect_identical(match_edge_correction("Ripley"), "isotropic")
  expect_identical(match_edge_correction("none"), "none")
  expect_identical(match_edge_correction("border"), "border")
  # v1.4.0 left `edge` undefined for an unrecognised string and failed later with
  # "object 'edge' not found"; fail at the boundary with a usable message instead.
  expect_error(match_edge_correction("periodic"), "Unsupported")
  expect_error(match_edge_correction("rigid"), "Unsupported")
  expect_error(match_edge_correction(c("trans", "iso")), "single string")
})

test_that("too-few-cell cases return NA rather than dividing by zero", {
  cs <- continuous_case()
  n <- spatstat.geom::npoints(cs$X)
  pairs <- k_pairs(cs$X, cs$r, "translation")
  one <- logical(n); one[1] <- TRUE
  expect_true(all(is.na(k_from_pairs(pairs, one))))            # n*(n-1) == 0
  none <- logical(n)
  expect_true(all(is.na(k_from_pairs(pairs, none))))           # no cells at all
})
