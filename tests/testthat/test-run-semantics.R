# Every function that writes a derived slot must agree on what `overwrite` means:
#   overwrite = TRUE  -> replace, Run = 1
#   overwrite = FALSE -> append,  Run = previous max + 1
#
# Two functions failed this in two different ways before 2.0.0, and both failed on
# their DEFAULT arguments, which is why nobody noticed:
#
#   dixons_s(overwrite = FALSE)        errored: max(Run) on a column that does not
#                                      exist yet ("object 'Run' not found")
#   marker_freq_diff(overwrite = FALSE) had its branches inverted -- it DESTROYED
#                                      the previous run on a populated mif, and on
#                                      a fresh one shipped Run = -Inf
#
# Both now share write_derived() with the seven metric functions. This file is the
# contract; add new derived-slot writers here rather than trusting them.

# Uses toy data with guaranteed disjoint marker counts rather than the example
# cores: this file is about bookkeeping, not statistics, and a fixture whose
# marker pair collapses under dual-positive removal would fail here for reasons
# that have nothing to do with `overwrite`. See mnames_bivariate() for why the
# obvious CD3/CD8 choice does not work.
run_mif <- function() {
  toy_mif(list(S1 = toy_spatial("S1", n = 250, markers = c(A = 50, B = 50)),
               S2 = toy_spatial("S2", n = 250, markers = c(A = 45, B = 55), seed = 12)))
}

# name -> (function of (mif, overwrite), derived slot it writes)
run_writers <- function() {
  m  <- c("A", "B")
  r  <- seq(0, 40, 10)
  list(
    ripleys_k = list(slot = "univariate_Count", f = function(mif, ow)
      ripleys_k(mif, mnames = m, r_range = r, permute = FALSE, workers = 1, overwrite = ow)),
    bi_ripleys_k = list(slot = "bivariate_Count", f = function(mif, ow)
      bi_ripleys_k(mif, mnames = m, r_range = r, permute = FALSE, workers = 1, overwrite = ow)),
    NN_G = list(slot = "univariate_NN", f = function(mif, ow)
      NN_G(mif, mnames = m, r_range = r, num_permutations = 2, workers = 1, overwrite = ow)),
    bi_NN_G = list(slot = "bivariate_NN", f = function(mif, ow)
      bi_NN_G(mif, mnames = m, r_range = r, num_permutations = 2, workers = 1, overwrite = ow)),
    pair_correlation = list(slot = "univariate_pair_correlation", f = function(mif, ow)
      pair_correlation(mif, mnames = m, r_range = r, num_permutations = 2, workers = 1, overwrite = ow)),
    bi_pair_correlation = list(slot = "bivariate_pair_correlation", f = function(mif, ow)
      bi_pair_correlation(mif, mnames = m, r_range = r, num_permutations = 2, workers = 1, overwrite = ow)),
    interaction_variable = list(slot = "interaction_variable", f = function(mif, ow)
      interaction_variable(mif, mnames = m, r_range = r, num_permutations = 2, workers = 1, overwrite = ow)),
    marker_freq_diff = list(slot = "frequency_difference", f = function(mif, ow)
      marker_freq_diff(mif, classifier = "Classifier.Label", ref_level = "Tumor",
                       diff_level = "Stroma", mnames = m, overwrite = ow)),
    dixons_s_Z = list(slot = "Dixon_Z", f = function(mif, ow)
      dixons_s(mif, mnames = m, num_permutations = 5, workers = 1, overwrite = ow)),
    dixons_s_C = list(slot = "Dixon_C", f = function(mif, ow)
      dixons_s(mif, mnames = m, num_permutations = 5, workers = 1, overwrite = ow))
  )
}

test_that("overwrite = FALSE works on a FRESH mif and yields Run = 1", {
  # This is the case dixons_s() and marker_freq_diff() got wrong, and it is the
  # default, so it is the first thing to assert.
  for (nm in names(run_writers())) {
    w <- run_writers()[[nm]]
    out <- suppressWarnings(suppressMessages(w$f(run_mif(), FALSE)))
    expect_s3_class(out, "mif")
    tbl <- out$derived[[w$slot]]
    expect_true(!is.null(tbl) && nrow(tbl) > 0, info = nm)
    expect_identical(unique(tbl$Run), 1, info = nm)
    expect_false(any(is.infinite(tbl$Run)), info = nm)
  }
})

test_that("overwrite = TRUE replaces and resets Run to 1", {
  for (nm in names(run_writers())) {
    w <- run_writers()[[nm]]
    once  <- suppressWarnings(suppressMessages(w$f(run_mif(), TRUE)))
    twice <- suppressWarnings(suppressMessages(w$f(once, TRUE)))
    expect_identical(unique(twice$derived[[w$slot]]$Run), 1, info = nm)
    expect_equal(nrow(twice$derived[[w$slot]]), nrow(once$derived[[w$slot]]), info = nm)
  }
})

test_that("overwrite = FALSE appends and increments Run", {
  for (nm in names(run_writers())) {
    w <- run_writers()[[nm]]
    first  <- suppressWarnings(suppressMessages(w$f(run_mif(), TRUE)))
    second <- suppressWarnings(suppressMessages(w$f(first, FALSE)))
    third  <- suppressWarnings(suppressMessages(w$f(second, FALSE)))
    n1 <- nrow(first$derived[[w$slot]])
    expect_identical(sort(unique(third$derived[[w$slot]]$Run)), c(1, 2, 3), info = nm)
    expect_equal(nrow(third$derived[[w$slot]]), 3 * n1, info = nm)
    # Appending must not invent or drop slots.
    expect_identical(sort(names(third$derived)), sort(names(first$derived)), info = nm)
  }
})

test_that("appending emits no warning about Run bookkeeping", {
  # `Run = -Inf` arrived with a max()/min() warning that was easy to miss. Scoped to
  # that class of warning rather than expect_no_warning(), because spatstat.explore
  # >= 3.8 emits its own one-time "default settings have changed" notice through the
  # pcf functions, which is not ours to swallow.
  bookkeeping <- "min|max|Inf|Run"
  for (nm in names(run_writers())) {
    w <- run_writers()[[nm]]
    first <- suppressWarnings(suppressMessages(w$f(run_mif(), TRUE)))
    got <- character(0)
    withCallingHandlers(
      suppressMessages(w$f(first, FALSE)),
      warning = function(cnd) {
        got <<- c(got, conditionMessage(cnd)); invokeRestart("muffleWarning")
      })
    expect_length(grep(bookkeeping, got), 0)
  }
})
