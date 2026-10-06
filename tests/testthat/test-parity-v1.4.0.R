# Holds 2.0.0 to the numbers 1.4.0 produced, for the quantities where that is a
# meaningful thing to ask -- and holds it to DIFFERING where 1.4.0 was wrong.
#
# This is contract-driven rather than hand-written. `tests/testthat/fixtures/
# parity-v1.4.0.rds` carries the 1.4.0 results, the call table that produced them,
# and a contract saying for each (entry, column) whether the two versions must
# agree, must differ, or cannot be compared. The loops below replay the calls and
# check each verdict. Adding a comparison means adding a contract row in
# tools/parity/compare.R and regenerating -- not editing this file.
#
# Why MUST_DIFFER rows are tested at all: a column where 1.4.0 was wrong would
# otherwise just be excluded, and nothing would notice the fix being reverted. Each
# MUST_DIFFER row also carries an `independent` note naming what pins the NEW value
# (usually a spatstat call), which is checked separately in test-k-engine.R and
# test-csr-null.R rather than re-derived here.
#
# The fixture and tools/parity/ are both in .Rbuildignore, so `R CMD check` on the
# tarball skips this entire file and neither ships to CRAN. It runs from a git
# checkout under devtools::test(). Regenerate with:
#
#     bash tools/parity/run.sh
#
# Three things to know before changing a verdict:
#
#   * 1.4.0's permutation paths were IRREPRODUCIBLE -- 24 nested mclapply() calls
#     omitted mc.cores, so set.seed() had no effect. Permuted CSR and anything
#     derived from it is NOT_COMPARABLE, not merely noisy.
#   * The "none" binning delta exists only for EVENLY SPACED r, because that is
#     when Kest takes its fast C path, and only for quantities Kest computes --
#     bi_ripleys_k's Observed K (from Kcross) agrees while its Exact CSR (from
#     Kest) differs, in the same call.
#   * `allow_na_mismatch` on a contract row means the finite values agree but the
#     two versions disagree about where the estimate stops. It relaxes MUST_MATCH
#     only; for a MUST_DIFFER column the NA pattern may BE the difference.

skip_if_not_installed("spatstat.geom")
skip_if_not_installed("spatstat.explore")

fixture <- function() {
  path <- test_path("fixtures", "parity-v1.4.0.rds")
  skip_if_not(file.exists(path), "v1.4.0 parity fixture not present (see tools/parity/)")
  readRDS(path)
}

# Rebuilt exactly as tools/parity/capture.R builds it.
parity_mif <- function(design) {
  create_mif(
    clinical_data = dplyr::mutate(example_clinical,
                                  deidentified_id = as.character(deidentified_id)),
    sample_data   = dplyr::mutate(example_summary,
                                  deidentified_id = as.character(deidentified_id)),
    spatial_list  = example_spatial[design$samples],
    patient_id    = "deidentified_id",
    sample_id     = "deidentified_sample")
}

#' Replay one captured call against the current code
replay <- function(entry, design) {
  mif <- parity_mif(design)
  set.seed(design$seed)
  out <- suppressWarnings(suppressMessages(
    do.call(get(entry$fn, mode = "function"), c(list(mif), entry$args))))
  got <- out$derived[[entry$slot]]
  if (is.null(got)) NULL else as.data.frame(got)
}

# Same alignment rule as tools/parity/compare.R: `iter` is NOT a key, because its
# value is one of the documented differences ("Estimater" vs "Estimate"), and every
# column under comparison is constant across iterations.
align_on <- function(o, n, keys) {
  k <- intersect(intersect(names(o), names(n)), keys)
  if (!length(k)) return(NULL)
  ok <- do.call(paste, c(o[k], sep = "\r")); nk <- do.call(paste, c(n[k], sep = "\r"))
  o <- o[!duplicated(ok), , drop = FALSE]; ok <- ok[!duplicated(ok)]
  n <- n[!duplicated(nk), , drop = FALSE]; nk <- nk[!duplicated(nk)]
  common <- intersect(ok, nk)
  if (!length(common)) return(NULL)
  list(o = o[match(common, ok), , drop = FALSE],
       n = n[match(common, nk), , drop = FALSE])
}

# Non-finite values must agree as categories: the Hanisch G estimator returns Inf at
# radii where its eroded area reaches zero, on both versions, and Inf - Inf is NaN.
finite_agree <- function(a, b, tol, allow_na = FALSE) {
  cls <- function(x) ifelse(is.na(x), NA_integer_,
                            ifelse(is.infinite(x), as.integer(sign(x)), 0L))
  if (!allow_na && any(is.na(a) != is.na(b))) return(FALSE)
  if (sum(cls(a) != cls(b), na.rm = TRUE) > 0L) return(FALSE)
  fin <- is.finite(a) & is.finite(b)
  !any(fin) || max(abs(a[fin] - b[fin])) <= tol
}


test_that("the parity fixture is what we think it is", {
  # A guard on the guard. If the fixture were regenerated from the wrong tree, or
  # the contract emptied, every test below would pass vacuously.
  fx <- fixture()
  expect_identical(attr(fx, "old_provenance")$spatialTIME_version, "1.4.0")
  expect_identical(attr(fx, "new_provenance")$spatialTIME_version, "2.0.0")
  expect_gt(length(attr(fx, "contract")), 20L)
  expect_gt(length(attr(fx, "calls")), 20L)
  expect_identical(attr(fx, "design")$samples,
                   c("TMA3_[9,K].tif", "TMA1_[3,B].tif"))
  # Both verdicts must actually be represented, or the contract has collapsed to
  # one kind of claim.
  ct <- attr(fx, "contract")
  expect_gt(sum(vapply(ct, function(r) length(r$match) > 0L, logical(1))), 10L)
  expect_gt(sum(vapply(ct, function(r) length(r$differ) > 0L, logical(1))), 5L)
})

test_that("every captured entry is covered by the contract or explicitly structural", {
  # Stops a new capture entry from being silently uncompared.
  fx <- fixture()
  ct <- attr(fx, "contract")
  structural <- c("subset_mif", "plot_immunoflo")
  uncovered <- setdiff(names(fx), c(names(ct), structural))
  expect_identical(uncovered, character(0))
})

test_that("MUST_MATCH columns still match 1.4.0", {
  fx <- fixture()
  ct <- attr(fx, "contract")
  design <- attr(fx, "design")
  calls <- attr(fx, "calls")
  names(calls) <- vapply(calls, `[[`, character(1), "id")
  tol <- attr(fx, "tolerance")
  keys <- attr(fx, "keys")

  checked <- 0L
  for (id in names(ct)) {
    rule <- ct[[id]]
    plain <- setdiff(rule$match, c("__counts__", "__pvalues__", "__structure__"))
    special <- intersect(rule$match, c("__counts__", "__pvalues__"))
    if (!length(plain) && !length(special)) next
    if (is.null(calls[[id]])) next

    old <- fx[[id]]
    if (!is.data.frame(old)) next
    new <- replay(calls[[id]], design)
    expect_true(is.data.frame(new), info = id)

    al <- align_on(old, new, keys)
    expect_false(is.null(al), info = paste(id, "- no rows aligned"))

    want <- plain
    if ("__counts__" %in% special) {
      num <- names(new)[vapply(new, is.numeric, logical(1))]
      want <- c(want, intersect(setdiff(grep("_p\\.value$", num, value = TRUE,
                                             invert = TRUE), want), names(old)))
    }
    for (cn in want) {
      expect_true(cn %in% names(old) && cn %in% names(new),
                  info = paste(id, cn, "present on both sides"))
      a <- al$o[[cn]]; b <- al$n[[cn]]
      if (is.numeric(a) && is.numeric(b)) {
        expect_true(finite_agree(a, b, tol, isTRUE(rule$allow_na_mismatch)),
                    info = paste0(id, " / ", cn, ": max|diff| = ",
                                  format(suppressWarnings(
                                    max(abs(a - b), na.rm = TRUE)), digits = 4)))
      } else {
        expect_identical(as.character(a), as.character(b), info = paste(id, cn))
      }
      checked <- checked + 1L
    }
  }
  # The loop must have done work; a contract that silently stopped matching
  # anything would otherwise pass.
  expect_gt(checked, 40L)
})

test_that("MUST_DIFFER columns still differ from 1.4.0", {
  # These are the fixes. A column here agreeing again means a correctness fix was
  # reverted -- the border estimator falling back to "none", the Fisher table
  # regaining its double-counted margin, interaction_variable scrambling its marks.
  fx <- fixture()
  ct <- attr(fx, "contract")
  design <- attr(fx, "design")
  calls <- attr(fx, "calls")
  names(calls) <- vapply(calls, `[[`, character(1), "id")
  tol <- attr(fx, "tolerance")
  keys <- attr(fx, "keys")

  checked <- 0L
  for (id in names(ct)) {
    rule <- ct[[id]]
    plain <- setdiff(rule$differ, c("__counts__", "__pvalues__", "__structure__"))
    wants_pvalues <- "__pvalues__" %in% rule$differ
    structural <- "__structure__" %in% rule$differ
    if (!length(plain) && !wants_pvalues && !structural) next
    if (is.null(calls[[id]])) next

    old <- fx[[id]]
    new <- tryCatch(replay(calls[[id]], design), error = function(e) e)

    if (structural) {
      # Either 1.4.0 could not produce a value at all, or the shapes differ.
      old_failed <- !is.data.frame(old)
      new_failed <- inherits(new, "error") || is.null(new)
      if (old_failed || new_failed) {
        expect_true(old_failed || new_failed, info = id)
      } else {
        differs <- nrow(old) != nrow(new) ||
          !all(names(old) %in% names(new))
        expect_true(differs, info = paste(id, "- shapes must differ"))
      }
      checked <- checked + 1L
      next
    }

    expect_true(is.data.frame(old), info = id)
    expect_false(inherits(new, "error"), info = id)
    al <- align_on(old, new, keys)
    expect_false(is.null(al), info = id)

    want <- plain
    if (wants_pvalues) {
      want <- c(want, intersect(grep("_p\\.value$", names(old), value = TRUE),
                                names(new)))
    }
    for (cn in want) {
      a <- al$o[[cn]]; b <- al$n[[cn]]
      expect_false(finite_agree(a, b, tol),
                   info = paste0(id, " / ", cn,
                                 ": declared MUST_DIFFER but agrees to ", tol))
      checked <- checked + 1L
    }
  }
  expect_gt(checked, 10L)
})

test_that("the column renames from 1.4.0 are all accounted for", {
  # Every old name in the map must be gone from the current schema, and every new
  # name present. This is what makes the rename table in the fixture trustworthy.
  fx <- fixture()
  map <- attr(fx, "rename_map")
  design <- attr(fx, "design")
  out <- ripleys_k(parity_mif(design), mnames = design$markers[1],
                   r_range = c(0, 10, 20), permute = TRUE, num_permutations = 2,
                   workers = 1, overwrite = TRUE)$derived$univariate_Count
  g <- bi_NN_G(parity_mif(design), mnames = design$pairs[[1]],
               r_range = c(0, 10, 20), num_permutations = 2, workers = 1,
               overwrite = TRUE)$derived$bivariate_NN
  current <- union(names(out), names(g))
  # Old names must not reappear. `From`/`To` are excluded: dixons_s legitimately
  # keeps them as Dixon's own from-type/to-type columns.
  gone <- setdiff(names(map), c("From", "To"))
  expect_identical(intersect(gone, current), character(0))
  expect_true(all(c("Theoretical CSR", "Permuted CSR", "Exact CSR", "Anchor",
                    "Counted", "Permutations Larger than Observed") %in% current))
})
