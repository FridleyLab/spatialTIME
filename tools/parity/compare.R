#!/usr/bin/env Rscript
#
# Diff two capture.R outputs against a declared contract, print a verdict table, and
# write the committed fixture.
#
#   Rscript tools/parity/compare.R <old.rds> <new.rds> --fixture=<path.rds>
#
# The contract says, for each (entry, column), whether the two versions MUST agree,
# MUST differ, or cannot be compared at all. This script exits non-zero if any
# measured result contradicts its declared verdict -- so the generator is itself a
# check, and "the fixture was regenerated" cannot quietly mean "the deltas changed".
#
# Why three verdicts rather than two
# ----------------------------------
# MUST_MATCH is the parity claim: the refactor was supposed to preserve this number.
# MUST_DIFFER is a deliberate fix: 1.4.0's value was wrong, and a test that merely
# ignored the column would not notice the fix being reverted. NOT_COMPARABLE is for
# quantities 1.4.0 could not produce reproducibly at all -- its permutation paths ran
# 24 nested mclapply() calls without mc.cores, so set.seed() had no effect and two
# runs at the same seed could differ by thousands. Comparing those is comparing
# against noise, and a tolerance loose enough to pass would be loose enough to hide
# anything.

args_cli <- commandArgs(trailingOnly = TRUE)
pos <- args_cli[!grepl("^--", args_cli)]
if (length(pos) < 2L) stop("Usage: compare.R <old.rds> <new.rds> [--fixture=<path>]",
                           call. = FALSE)
FIXTURE <- {
  hit <- grep("^--fixture=", args_cli, value = TRUE)
  if (length(hit)) sub("^--fixture=", "", hit[[1]]) else NA_character_
}

old <- readRDS(pos[[1]])
new <- readRDS(pos[[2]])
stopifnot(identical(attr(old, "side"), "old"), identical(attr(new, "side"), "new"))

# ---- old -> new column names ----------------------------------------------
RENAME <- c(
  "Theoretical G" = "Theoretical CSR", "Theoretical g" = "Theoretical CSR",
  "Permuted G" = "Permuted CSR", "Permuted g" = "Permuted CSR",
  "Permuted Interaction" = "Permuted CSR",
  "From" = "Anchor", "To" = "Counted",
  "Degree of Correlation Theoretical" = "Degree of Clustering Theoretical",
  "Degree of Correlation Permuted" = "Degree of Clustering Permutation",
  "Degree of Interaction Permuted" = "Degree of Clustering Permutation",
  "Permuted_larger_than_Observed" = "Permutations Larger than Observed"
)
# From/To -> Anchor/Counted is right for interaction_variable and
# bi_pair_correlation, where those columns ARE the pair identifiers. It is wrong for
# dixons_s: there `From`/`To` are Dixon's own output columns (the from-type and
# to-type of a count), present under those names on BOTH sides, and 2.0.0 ADDED
# Anchor/Counted alongside them as the pair identifiers. Renaming old's From to
# Anchor therefore aliases two different quantities and the join matches almost
# nothing -- it cut dixons_s.Z from 24 rows to 6 "common" ones and made every column
# look like a huge regression.
NO_FROM_RENAME <- c("dixons_s.Z", "dixons_s.C")

rename_old <- function(d, id = NULL) {
  if (!is.data.frame(d)) return(d)
  map <- RENAME
  if (!is.null(id) && id %in% NO_FROM_RENAME) map <- map[setdiff(names(map), c("From", "To"))]
  hit <- names(d) %in% names(map)
  names(d)[hit] <- map[names(d)[hit]]
  # 1.4.0 de-spaced only tablaZ's column names, so Dixon_C arrived with padded ones
  # ("  df ", "  P.rand") and, where a success row was bound to an early-return row,
  # BOTH "  P.rand" and "P.rand". Trimming is the honest mapping of that 2.0.0 fix;
  # the duplicate it creates is collapsed to the column that actually holds values.
  if (!is.null(id) && id %in% NO_FROM_RENAME) {
    names(d) <- trimws(names(d))
    if (anyDuplicated(names(d))) {
      keep <- !duplicated(names(d)) |
        vapply(seq_along(d), function(i) !all(is.na(d[[i]])), logical(1))
      dup <- names(d)[duplicated(names(d))]
      for (nm in unique(dup)) {
        at <- which(names(d) == nm)
        best <- at[which.max(vapply(at, function(i) sum(!is.na(d[[i]])), numeric(1)))]
        keep[setdiff(at, best)] <- FALSE
      }
      d <- d[, keep, drop = FALSE]
    }
  }
  d
}

# ---- the contract ---------------------------------------------------------
# Columns 1.4.0 could not produce reproducibly, in any entry.
NOT_COMPARABLE_COLS <- c("Permuted CSR", "Degree of Clustering Permutation", "iter")

rule <- function(entry, match = character(0), differ = character(0),
                 note = "", independent = NA_character_,
                 allow_na_mismatch = FALSE) {
  list(entry = entry, match = match, differ = differ, note = note,
       independent = independent, allow_na_mismatch = allow_na_mismatch)
}

K_MATCH <- c("Observed K", "Theoretical CSR", "Exact CSR",
             "Degree of Clustering Theoretical", "Degree of Clustering Exact")

CONTRACT <- list(
  rule("ripleys_k.translation", match = K_MATCH,
       note = "the parity claim: translation K is unchanged"),
  rule("ripleys_k.isotropic", match = K_MATCH,
       note = "the parity claim: isotropic K is unchanged"),
  rule("ripleys_k.none", match = c("Observed K", "Theoretical CSR", "Exact CSR"),
       note = paste("Agrees here, and that is informative: Kest only takes its fast",
                    "C path when r is EVENLY SPACED, and this entry's r_range is not",
                    "(it carries the truncation radii). So 1.4.0 routed through whist",
                    "too and the binning difference does not arise. See",
                    "ripleys_k.none.even for the case where it does.")),
  rule("ripleys_k.none.even", differ = c("Observed K", "Exact CSR"),
       match = c("Theoretical CSR"),
       note = paste("Evenly spaced r, so 1.4.0's Kest took its fast C path, which",
                    "bins d <= r (right-closed) where whist bins left-closed. Those",
                    "disagree exactly when a pair distance lands on a bin edge --",
                    "never on continuous coordinates, constantly on the half-integer",
                    "cell centres HALO and Vectra emit. 2.0.0 always uses whist, so",
                    "all four corrections share one convention."),
       independent = "Kest(correction = c('none','translation'))$un"),
  rule("ripleys_k.border", match = c("Observed K", "Theoretical CSR"),
       differ = c("Exact CSR", "Degree of Clustering Exact"),
       note = paste("Observed K must MATCH -- 1.4.0 passed border through to Kest,",
                    "and 2.0.0's reduced-sample engine reproduces it. Exact CSR must",
                    "DIFFER: 1.4.0 reported the K of all cells, which is not E[K] of",
                    "a subset for border (measured 1% low), so 2.0.0 reports NA."),
       independent = "Kest(correction = c('border','translation'))$border"),
  rule("ripleys_k.big500", differ = c("Observed K"),
       match = c("Theoretical CSR"),
       note = paste("1.4.0 silently replaced the requested correction with 'none'",
                    "above `big` AND took the tiled getTile() branch. 2.0.0 honours",
                    "the request; `big` is memory-only.")),
  rule("ripleys_k.perm", match = c("Theoretical CSR"),
       note = "permutation path: numerically incomparable, pins shape only"),
  rule("ripleys_k.sparse",
       match = c("Observed K", "Theoretical CSR",
                 "Degree of Clustering Theoretical", "Degree of Clustering Exact"),
       differ = c("Exact CSR"), allow_na_mismatch = TRUE,
       note = paste("a marker with 1 positive. Both versions keep the shape and emit",
                    "NA rows, and 1.4.0 labelled iter 'Estimator' here against",
                    "'Estimater' in its own main branch. Exact CSR DIFFERS: 1.4.0",
                    "still reported the K of all cells (which does not depend on the",
                    "marker), 2.0.0 NAs the whole unestimable row. Defensible either",
                    "way; recorded so the change is not silent.")),
  rule("bi_ripleys_k.translation", match = K_MATCH, allow_na_mismatch = TRUE,
       note = paste("the parity claim for the cross-K. NA mismatch allowed: 2.0.0",
                    "truncates at rmax_valid (1276 on the dense core) where 1.4.0's",
                    "Kcross path reported a number, so the two disagree about WHERE",
                    "the estimate stops. Every finite value still agrees to 7e-10.")),
  rule("bi_ripleys_k.none", match = c("Observed K", "Theoretical CSR", "Exact CSR"),
       note = "as ripleys_k.none: uneven r, so both sides used whist"),
  rule("bi_ripleys_k.none.even", match = c("Observed K", "Theoretical CSR"),
       differ = c("Exact CSR"),
       note = paste("This entry locates the binning difference exactly. On an even",
                    "grid, Observed K AGREES while Exact CSR DIFFERS -- within the",
                    "same call. The fast C path is Kest's alone: Kcross/Kmulti have",
                    "none and always route through whist, so 1.4.0's cross-K already",
                    "used 2.0.0's convention. But 1.4.0 computed the bivariate Exact",
                    "CSR with Kest() over all cells, which on an even grid DID take",
                    "the fast path. So the delta follows the function, not the",
                    "statistic."),
       independent = "Kest(correction = c('none','translation'))$un"),
  rule("bi_ripleys_k.border", match = c("Observed K", "Theoretical CSR"),
       differ = c("Exact CSR", "Degree of Clustering Exact"),
       allow_na_mismatch = TRUE, note = "as ripleys_k.border"),
  rule("NN_G.rs",   match = c("Observed G", "Theoretical CSR"),
       note = "delegating to Gest was numerically free"),
  rule("NN_G.han",  match = c("Observed G", "Theoretical CSR"), note = "as NN_G.rs"),
  rule("NN_G.none", match = c("Observed G", "Theoretical CSR"), note = "as NN_G.rs"),
  rule("NN_G.km",   differ = c("__structure__"),
       note = paste("1.4.0's output was MALFORMED: it reordered columns by position",
                    "assuming Gest returns three, but 'km' returns five, so the table",
                    "had no sample-id column, no Marker column and a leaked",
                    "hazard/theohaz. There is no valid 1.4.0 value to compare."),
       independent = "Gest(correction = 'km')$km"),
  rule("bi_NN_G.rs",  match = c("Observed G", "Theoretical CSR"),
       note = "replacing the hand-rolled estimators with Gcross was numerically free"),
  rule("bi_NN_G.han", match = c("Observed G", "Theoretical CSR"), note = "as bi_NN_G.rs"),
  rule("bi_NN_G.km",  differ = c("__structure__"),
       note = "1.4.0 had no km branch at all; G_cross_df2 was undefined and it errored"),
  rule("pair_correlation", match = c("Observed g", "Theoretical CSR"),
       note = "pcf estimator unchanged"),
  rule("pair_correlation.sparse", differ = c("__structure__"),
       note = paste("markers with < 3 positives: 2.0.0 returns NULL and DROPS the",
                    "row, 1.4.0 had no guard and emitted rows. Row counts differ by",
                    "design, so the shape is the comparison.")),
  rule("bi_pair_correlation", match = c("Observed g", "Theoretical CSR"),
       note = "pcfcross estimator unchanged"),
  rule("interaction_variable", differ = c("Observed Interaction"),
       note = paste("1.4.0's values are WRONG, and this was not previously known.",
                    "get_bi_rows() returns rows marker-major (every anchor cell, then",
                    "every counted cell), so cells$cell is not globally ascending --",
                    "but subset(ppp, cells$cell) returns points in SORTED order, and",
                    "marks(ps) <- cells$Marker then attached the marker labels in",
                    "marker-major order to points in index order. 1.4.0 therefore",
                    "measured distances between SCRAMBLED marker sets. Verified by",
                    "brute force on example_spatial[['TMA3_[9,K].tif']] with",
                    "FOXP3/CD8 at r=20: exactly 4 of 109 anchors are within 20 units,",
                    "so 3.669725 is correct; 1.4.0 reported 4.587156 (5/109), and",
                    "re-running its exact code path reproduces that."),
       independent = "brute-force nearest-neighbour count over the two marker sets"),
  rule("interaction_variable.nested", differ = c("Observed Interaction"),
       note = "as interaction_variable; near-total marker nesting"),
  rule("interaction_variable.sparse", differ = c("__structure__"),
       note = paste("one core has zero counted cells: 2.0.0 drops that pair,",
                    "1.4.0 emitted a stub")),
  rule("dixons_s.Z", match = c("Obs.Count", "Exp.Count", "S", "Z", "p-val.Z"),
       note = "dixon::dixon() is untouched; only the surrounding bookkeeping changed"),
  rule("dixons_s.C", match = c("df", "Chi-sq", "P.asymp"),
       note = "asymptotic columns only; P.rand is a permutation result"),
  rule("marker_freq_diff", differ = c("__pvalues__"),
       match = c("__counts__"),
       note = paste("every p-value 1.4.0 produced was wrong: the Fisher table used",
                    "the compartment TOTAL as its second row, double-counting the",
                    "positives, and marker selection used a substring match so",
                    "CD3..CD8. absorbed CD3..CD8..FOXP3."),
       independent = "fisher.test on positives-vs-negatives by level")
)
names(CONTRACT) <- vapply(CONTRACT, `[[`, character(1), "entry")

# dixon's permutation columns, like every other permutation result.
NOT_COMPARABLE_COLS <- c(NOT_COMPARABLE_COLS, "P.rand", "p-val.Nobs", "Simulations")

# ---- comparison -----------------------------------------------------------
# `iter` is deliberately NOT a key. Its VALUE is one of the documented differences
# -- 1.4.0 wrote "Estimater" (sic) on the exact path and "Estimator" in its sparse
# branch, against "Estimate" in 2.0.0, and bivariate K used bare numbers -- so keying
# on it matches nothing at all. Every column declared MUST_MATCH is constant across
# iterations (observed and theoretical values are recycled over the permutations), so
# collapsing to the first row per key is exact for them, and everything that does
# vary by iteration is NOT_COMPARABLE anyway.
KEYS <- c("deidentified_sample", "Marker", "Anchor", "Counted", "r",
          "From", "To", "Direction")

#' Compare two numeric columns, tolerating legitimately non-finite values
#'
#' The Hanisch G estimator returns `Inf` at radii where its eroded area reaches
#' zero, on both versions and at the same radii. `Inf - Inf` is `NaN`, which
#' propagates through `max()` and made an exactly-agreeing column report as
#' differing. So finite positions are compared numerically and non-finite ones are
#' required to agree exactly, as categories.
num_diff <- function(a, b) {
  na_mismatch <- sum(is.na(a) != is.na(b))
  # Classify every non-NA position: -Inf, finite, or +Inf. These must agree.
  cls <- function(x) ifelse(is.na(x), NA_integer_,
                            ifelse(is.infinite(x), as.integer(sign(x)), 0L))
  inf_mismatch <- sum(cls(a) != cls(b), na.rm = TRUE)
  fin <- is.finite(a) & is.finite(b)
  list(n_compared = sum(fin),
       max_abs = if (any(fin)) max(abs(a[fin] - b[fin])) else NA_real_,
       max_rel = if (any(fin)) max(abs(a[fin] - b[fin]) /
                                     pmax(abs(b[fin]), .Machine$double.eps)) else NA_real_,
       na_mismatch = na_mismatch,
       inf_mismatch = inf_mismatch)
}

#' Does a numeric comparison count as agreement?
#'
#' Agreement means the NA pattern matches, the infinity pattern matches, and the
#' finite values agree to `TOL`. A column that is entirely non-finite on both sides
#' and matches categorically agrees even though there is nothing to subtract.
num_agree <- function(d, tol, allow_na = FALSE) {
  (allow_na || d$na_mismatch == 0L) && d$inf_mismatch == 0L &&
    (d$n_compared == 0L || (!is.na(d$max_abs) && d$max_abs <= tol))
}

#' Line two tables up on whatever key columns they share
align <- function(o, n) {
  k <- intersect(intersect(names(o), names(n)), KEYS)
  if (!length(k)) return(NULL)
  ok <- do.call(paste, c(o[k], sep = "\r"))
  nk <- do.call(paste, c(n[k], sep = "\r"))
  # First row per key: see the note on KEYS above.
  o1 <- !duplicated(ok); n1 <- !duplicated(nk)
  o <- o[o1, , drop = FALSE]; ok <- ok[o1]
  n <- n[n1, , drop = FALSE]; nk <- nk[n1]
  common <- intersect(ok, nk)
  if (!length(common)) return(NULL)
  list(o = o[match(common, ok), , drop = FALSE],
       n = n[match(common, nk), , drop = FALSE],
       keys = k, n_common = length(common),
       n_old_only = sum(!ok %in% nk), n_new_only = sum(!nk %in% ok))
}

TOL <- 1e-9
report <- list(); failures <- character(0)

note_failure <- function(...) failures <<- c(failures, paste0(...))

for (id in names(new)) {
  o <- rename_old(old[[id]], id); n <- new[[id]]
  ct <- CONTRACT[[id]]
  oerr <- inherits(o, "parity_error") || inherits(o, "parity_absent")
  nerr <- inherits(n, "parity_error") || inherits(n, "parity_absent")

  if (is.null(ct)) {
    # subset_mif / plot_immunoflo are structural; reported, not contract-checked.
    report[[id]] <- list(kind = "structural",
                         old = if (oerr) o$reason else "captured",
                         new = if (nerr) n$reason else "captured")
    next
  }

  if (oerr || nerr) {
    structural_expected <- "__structure__" %in% ct$differ
    report[[id]] <- list(kind = "error",
                         old = if (oerr) substr(o$reason, 1, 90) else "ok",
                         new = if (nerr) substr(n$reason, 1, 90) else "ok",
                         verdict = if (structural_expected) "MUST_DIFFER (as declared)"
                                   else "UNEXPECTED")
    if (!structural_expected) {
      note_failure(id, ": one side errored but the contract did not declare a ",
                   "structural difference")
    }
    next
  }

  al <- align(o, n)
  rows <- list(shape_old = dim(o), shape_new = dim(n),
               common = if (is.null(al)) 0L else al$n_common,
               old_only = if (is.null(al)) nrow(o) else al$n_old_only,
               new_only = if (is.null(al)) nrow(n) else al$n_new_only)

  # Declared structural difference: assert the shape really does differ.
  if ("__structure__" %in% ct$differ) {
    same_shape <- identical(dim(o)[[1]], dim(n)[[1]]) &&
      identical(sort(intersect(names(o), names(n))), sort(names(o)))
    if (same_shape) {
      note_failure(id, ": declared a structural difference, but old and new have ",
                   "the same rows and old's columns are a subset of new's")
    }
    report[[id]] <- c(list(kind = "structure", verdict = "MUST_DIFFER"), rows)
    next
  }

  cols <- list()
  check_cols <- function(which, verdict) {
    for (cn in which) {
      if (cn == "__counts__") {
        cn_list <- grep("_p\\.value$", names(n), value = TRUE, invert = TRUE)
        cn_list <- intersect(cn_list, names(o))
        cn_list <- cn_list[vapply(o[cn_list], is.numeric, logical(1))]
      } else if (cn == "__pvalues__") {
        cn_list <- intersect(grep("_p\\.value$", names(o), value = TRUE), names(n))
      } else {
        cn_list <- cn
      }
      for (c1 in cn_list) {
        if (!c1 %in% names(o) || !c1 %in% names(n)) {
          cols[[c1]] <<- list(verdict = verdict, status = "ABSENT",
                              where = if (!c1 %in% names(o)) "old" else "new")
          if (verdict == "MUST_MATCH") {
            note_failure(id, ": column \"", c1, "\" declared MUST_MATCH but is ",
                         "absent from the ", if (!c1 %in% names(o)) "old" else "new",
                         " side")
          }
          next
        }
        if (is.null(al)) {
          cols[[c1]] <<- list(verdict = verdict, status = "NO_COMMON_ROWS")
          if (verdict == "MUST_MATCH") {
            note_failure(id, ": column \"", c1, "\" declared MUST_MATCH but the two ",
                         "sides share no key rows")
          }
          next
        }
        a <- al$o[[c1]]; b <- al$n[[c1]]
        if (!is.numeric(a) || !is.numeric(b)) {
          agree <- identical(as.character(a), as.character(b))
          cols[[c1]] <<- list(verdict = verdict, status = if (agree) "IDENTICAL" else "DIFFERS")
          if (verdict == "MUST_MATCH" && !agree) {
            note_failure(id, ": column \"", c1, "\" declared MUST_MATCH but differs")
          }
          next
        }
        d <- num_diff(a, b)
        # allow_na_mismatch relaxes MUST_MATCH only. For a MUST_DIFFER column the
        # NA pattern IS the difference (border's Exact CSR is NA throughout in
        # 2.0.0), so relaxing it there would report the column as agreeing.
        agree <- num_agree(d, TOL,
                           allow_na = isTRUE(ct$allow_na_mismatch) &&
                                      verdict == "MUST_MATCH")
        cols[[c1]] <<- c(list(verdict = verdict,
                              status = if (agree) "MATCH" else "DIFFER"), d)
        if (verdict == "MUST_MATCH" && !agree) {
          note_failure(id, ": column \"", c1, "\" declared MUST_MATCH but max|diff| = ",
                       format(d$max_abs, digits = 4), " over ", d$n_compared,
                       " finite rows (na_mismatch=", d$na_mismatch,
                       ", inf_mismatch=", d$inf_mismatch, ")")
        }
        if (verdict == "MUST_DIFFER" && agree) {
          note_failure(id, ": column \"", c1, "\" declared MUST_DIFFER but the two ",
                       "sides agree to ", TOL)
        }
      }
    }
  }
  check_cols(ct$match, "MUST_MATCH")
  check_cols(ct$differ, "MUST_DIFFER")

  report[[id]] <- c(list(kind = "columns", note = ct$note), rows, list(cols = cols))
}

# ---- print ----------------------------------------------------------------
cat("\n", strrep("=", 78), "\n", sep = "")
cat("PARITY: spatialTIME ", attr(old, "provenance")$spatialTIME_version,
    "  ->  ", attr(new, "provenance")$spatialTIME_version, "\n", sep = "")
cat(strrep("=", 78), "\n\n", sep = "")

for (id in names(report)) {
  r <- report[[id]]
  cat(id, "\n", sep = "")
  if (r$kind == "error") {
    cat("    ", r$verdict, "   old: ", r$old, "\n", sep = "")
    cat("                       new: ", r$new, "\n", sep = "")
  } else if (r$kind == "structural") {
    cat("    structural capture (not contract-checked)\n")
  } else if (r$kind == "structure") {
    cat(sprintf("    MUST_DIFFER (shape)  old %dx%d  new %dx%d  common=%d old_only=%d new_only=%d\n",
                r$shape_old[1], r$shape_old[2], r$shape_new[1], r$shape_new[2],
                r$common, r$old_only, r$new_only))
  } else {
    cat(sprintf("    rows: old %d, new %d, common %d\n",
                r$shape_old[1], r$shape_new[1], r$common))
    for (cn in names(r$cols)) {
      cc <- r$cols[[cn]]
      extra <- if (!is.null(cc$n_compared))
        sprintf("max|d|=%-10.4g rel=%-10.4g n=%-4d na_mm=%d inf_mm=%d",
                cc$max_abs, cc$max_rel, cc$n_compared, cc$na_mismatch,
                cc$inf_mismatch) else ""
      cat(sprintf("      %-36s %-12s %-9s %s\n", cn, cc$verdict, cc$status, extra))
    }
  }
  cat("\n")
}

if (length(failures)) {
  cat(strrep("!", 78), "\n", sep = "")
  cat("CONTRACT VIOLATIONS (", length(failures), ")\n\n", sep = "")
  for (f in failures) cat("  * ", f, "\n", sep = "")
  cat(strrep("!", 78), "\n", sep = "")
} else {
  cat("All ", length(report), " entries agree with the declared contract.\n", sep = "")
}

# ---- fixture --------------------------------------------------------------
if (!is.na(FIXTURE)) {
  fx <- stats::setNames(lapply(names(old), function(i) rename_old(old[[i]], i)),
                        names(old))
  attr(fx, "contract") <- CONTRACT
  attr(fx, "rename_map") <- RENAME
  attr(fx, "not_comparable_cols") <- NOT_COMPARABLE_COLS
  attr(fx, "tolerance") <- TOL
  attr(fx, "keys") <- KEYS
  attr(fx, "design") <- attr(old, "design")
  attr(fx, "calls") <- attr(old, "calls")
  attr(fx, "old_provenance") <- attr(old, "provenance")
  attr(fx, "new_provenance") <- attr(new, "provenance")
  dir.create(dirname(FIXTURE), recursive = TRUE, showWarnings = FALSE)
  saveRDS(fx, FIXTURE)
  cat("\nwrote fixture: ", FIXTURE, " (",
      format(file.size(FIXTURE), big.mark = ","), " bytes)\n", sep = "")
}

if (length(failures)) quit(status = 1L)
