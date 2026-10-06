#!/usr/bin/env Rscript
#
# Capture one side of the cross-branch parity comparison.
#
#   Rscript tools/parity/capture.R --pkg=<dir> --side=old|new --out=<file.rds>
#
# Loads ONE package tree with pkgload and runs a fixed table of calls against the
# shipped example data, writing the results plus full provenance to an RDS.
# tools/parity/compare.R then diffs two of those and writes the committed fixture.
#
# Why a separate process per side
# -------------------------------
# Both trees declare `Package: spatialTIME`, so one R session cannot hold both
# namespaces regardless of `lib.loc`. Installing them into separate libraries would
# not help -- it still needs two sessions -- and would add a third installed copy
# alongside the CRAN one already in the system library, making "which did I just
# measure?" a real hazard. Two processes, no installs.
#
# The interlock below exists because of that hazard: a stray `spatialTIME::` or a
# failed load would silently measure the installed CRAN build instead of the tree
# under test, and the comparison would look fine. DO NOT write `spatialTIME::`
# anywhere in this file -- use bare names, which resolve into the loaded namespace.

suppressWarnings(suppressMessages({
  args_cli <- commandArgs(trailingOnly = TRUE)
}))

arg_of <- function(flag, default = NULL) {
  hit <- grep(paste0("^--", flag, "="), args_cli, value = TRUE)
  if (!length(hit)) {
    if (is.null(default)) stop("Missing required argument --", flag, call. = FALSE)
    return(default)
  }
  sub(paste0("^--", flag, "="), "", hit[[1]])
}

PKG  <- normalizePath(arg_of("pkg"), mustWork = TRUE)
SIDE <- arg_of("side")
OUT  <- arg_of("out")
stopifnot(SIDE %in% c("old", "new"))

EXPECT_VERSION <- c(old = "1.4.0", new = "2.0.0")

message("capture: side=", SIDE, " pkg=", PKG)
suppressMessages(pkgload::load_all(PKG, quiet = TRUE, export_all = TRUE,
                                   helpers = FALSE, attach_testthat = FALSE))

# ---- interlock ------------------------------------------------------------
# `helpers = FALSE` matters: both trees have tests/testthat/helper-mif.R, and
# loading one side's helpers would shadow functions used below.
got_version <- as.character(utils::packageVersion("spatialTIME"))
if (!identical(got_version, EXPECT_VERSION[[SIDE]])) {
  stop("Side \"", SIDE, "\" should be spatialTIME ", EXPECT_VERSION[[SIDE]],
       " but the loaded namespace reports ", got_version,
       ".\n  Did pkgload fall back to the installed package?", call. = FALSE)
}
loaded_from <- normalizePath(pkgload::pkg_path(PKG))
if (!identical(loaded_from, PKG)) {
  stop("Loaded ", loaded_from, " but was asked for ", PKG, call. = FALSE)
}
message("  spatialTIME ", got_version, " from ", loaded_from)

# ---- fixture --------------------------------------------------------------
# Two cores, unthinned.
#   TMA3_[9,K].tif  1803 cells, 536 CD3+ / 83 CD8+ / 34 FOXP3+ -- the only core with
#                   usable positives, so it carries the numeric comparisons.
#   TMA1_[3,B].tif  3803 cells, 17 / 7 / 3 -- the degenerate case, and where the
#                   row-count differences actually live: pair_correlation and
#                   interaction_variable DROP sparse markers in 2.0.0 where 1.4.0
#                   emitted stub rows.
SAMPLES <- c("TMA3_[9,K].tif", "TMA1_[3,B].tif")
MARKERS <- c("CD3..Opal.570..Positive", "CD8..Opal.520..Positive",
             "FOXP3..Opal.620..Positive")
# CD8/FOXP3 is the only workable real pair (79 vs 30 disjoint on the dense core);
# CD3/CD8 is near-total nesting, which pins the one-counted-cell degeneracy.
PAIRS   <- list(c("CD8..Opal.520..Positive", "FOXP3..Opal.620..Positive"),
                c("CD3..Opal.570..Positive", "CD8..Opal.520..Positive"))
# 0-60 keeps comparability with the existing narrow fixture. 650/700 straddle the
# dense core's isotropic rmax_valid (689.8) and 1300 crosses translation's (1276.0);
# those three radii are the only thing that pins 2.0.0's NA truncation.
R_RANGE <- c(0, 10, 20, 30, 40, 50, 60, 650, 700, 1300)
# spatstat's pcf requires evenly spaced r ("r values must be evenly spaced"), so the
# pair correlation entries cannot take the truncation radii above. They get the even
# grid instead -- which costs nothing, because pcf has no rmax_valid truncation to
# pin and its Exact CSR is NA on both sides anyway.
PCF_R   <- seq(0, 60, 10)
# Markers sparse enough to trip the degenerate-case guards, which differ between the
# two versions in SHAPE and not only in value. On TMA3_[9,K] / TMA1_[3,B]:
#   CD3..PD.L1.  1 and 0 positives -> under pair_correlation's `n_pos < 3`
#   CD8 / PD1    PD1 is 5 and 0 -> under interaction_variable's `n_j < 1` on the
#                second core only, so one core survives and one is dropped
# 2.0.0 returns NULL for these and drops the row entirely; 1.4.0 emitted stub rows.
# Nothing in MARKERS triggers it -- the rarest there is FOXP3 at 3, and the guard is
# `< 3`, so the row-count delta would otherwise go unexercised.
#
# The bivariate pair has to avoid nesting as well as being sparse. Every bivariate
# measure here drops cells positive for both markers, and on this data PD1, FOXP3 and
# PD-L1 are all entirely inside CD3 -- so CD3/PD1 leaves ZERO counted cells on both
# cores and the pair fails everywhere instead of on one core. CD8 is the only marker
# with positives outside CD3 (15 of 83 on the dense core), and no PD1+ cell is CD8+,
# so CD8/PD1 gives 83 against 5 on the dense core and 7 against 0 on the sparse one.
SPARSE_UNI  <- c("CD3..Opal.570..Positive", "CD3..PD.L1.")
SPARSE_PAIR <- c("CD8..Opal.520..Positive", "PD1..Opal.650..Positive")
SEED    <- 8675309

build_mif <- function() {
  # as.character() on the ids because 1.4.0's create_mif() computed a discarded
  # full_join() that failed on mismatched integer/character id types -- which is why
  # every pre-2.0.0 example carried this mutate.
  create_mif(
    clinical_data = dplyr::mutate(example_clinical,
                                  deidentified_id = as.character(deidentified_id)),
    sample_data   = dplyr::mutate(example_summary,
                                  deidentified_id = as.character(deidentified_id)),
    spatial_list  = example_spatial[SAMPLES],
    patient_id    = "deidentified_id",
    sample_id     = "deidentified_sample")
}

# ---- the call table -------------------------------------------------------
# One table for both sides. Anything a side cannot accept is dropped by adapt().
uni <- function(id, slot, ...) {
  list(id = id, fn = "ripleys_k", slot = slot,
       args = c(list(mnames = MARKERS, r_range = R_RANGE, workers = 1,
                     overwrite = TRUE), list(...)))
}
bi <- function(id, fn, slot, ...) {
  list(id = id, fn = fn, slot = slot,
       args = c(list(mnames = PAIRS[[1]], r_range = R_RANGE, workers = 1,
                     overwrite = TRUE), list(...)))
}

CALLS <- c(
  # Ripley's K, every correction, exact path. These carry the MUST_MATCH claims.
  lapply(c("translation", "isotropic", "none", "border"), function(ec)
    uni(paste0("ripleys_k.", ec), "univariate_Count",
        edge_correction = ec, permute = FALSE)),
  list(
    # 1.4.0: exceeding `big` silently replaced the correction with "none" AND took
    # the tiled getTile() branch, so one entry exercises both legacy behaviours.
    # Both cores are over 500.
    uni("ripleys_k.big500", "univariate_Count",
        edge_correction = "translation", permute = FALSE, big = 500),
    # Permutation entries are NOT_COMPARABLE numerically (1.4.0's RNG was
    # irreproducible) but pin iter values, column names/order, class and row counts.
    uni("ripleys_k.perm", "univariate_Count",
        edge_correction = "translation", permute = TRUE, num_permutations = 5,
        keep_permutation_distribution = TRUE),
    # The "none" binning difference only exists when r is EVENLY SPACED. Kest picks
    # its fast C path (which bins d <= r) via will.do.fast, and that requires evenly
    # spaced r; with R_RANGE's truncation radii it falls back to whist and agrees
    # with 2.0.0 exactly. So the delta needs its own entry on an even grid --
    # without it the documented difference goes unexercised and the contract would
    # record "none agrees", which is true only for uneven r.
    list(id = "ripleys_k.none.even", fn = "ripleys_k", slot = "univariate_Count",
         args = list(mnames = MARKERS, r_range = PCF_R, edge_correction = "none",
                     permute = FALSE, workers = 1, overwrite = TRUE)),
    list(id = "bi_ripleys_k.none.even", fn = "bi_ripleys_k",
         slot = "bivariate_Count",
         args = list(mnames = PAIRS[[1]], r_range = PCF_R,
                     edge_correction = "none", permute = FALSE, workers = 1,
                     overwrite = TRUE))
  ),
  lapply(c("translation", "none", "border"), function(ec)
    bi(paste0("bi_ripleys_k.", ec), "bi_ripleys_k", "bivariate_Count",
       edge_correction = ec, permute = FALSE)),
  lapply(c("rs", "han", "none", "km"), function(ec)
    list(id = paste0("NN_G.", ec), fn = "NN_G", slot = "univariate_NN",
         args = list(mnames = MARKERS, r_range = R_RANGE, edge_correction = ec,
                     num_permutations = 5, keep_permutation_distribution = TRUE,
                     workers = 1, overwrite = TRUE))),
  lapply(c("rs", "han", "km"), function(ec)
    list(id = paste0("bi_NN_G.", ec), fn = "bi_NN_G", slot = "bivariate_NN",
         args = list(mnames = PAIRS[[1]], r_range = R_RANGE, edge_correction = ec,
                     num_permutations = 5, keep_permutation_distribution = TRUE,
                     workers = 1, overwrite = TRUE))),
  list(
    list(id = "pair_correlation", fn = "pair_correlation",
         slot = "univariate_pair_correlation",
         args = list(mnames = MARKERS, r_range = PCF_R, num_permutations = 5,
                     workers = 1, overwrite = TRUE)),
    list(id = "bi_pair_correlation", fn = "bi_pair_correlation",
         slot = "bivariate_pair_correlation",
         args = list(mnames = PAIRS[[1]], r_range = PCF_R, num_permutations = 5,
                     workers = 1, overwrite = TRUE)),
    list(id = "interaction_variable", fn = "interaction_variable",
         slot = "interaction_variable",
         args = list(mnames = PAIRS[[1]], r_range = R_RANGE, num_permutations = 5,
                     workers = 1, overwrite = TRUE)),
    # The nested pair is where 1.4.0's stub rows differ from 2.0.0's drops.
    list(id = "interaction_variable.nested", fn = "interaction_variable",
         slot = "interaction_variable",
         args = list(mnames = PAIRS[[2]], r_range = R_RANGE, num_permutations = 5,
                     workers = 1, overwrite = TRUE)),
    # Degenerate-marker shape. 2.0.0 drops these rows; 1.4.0 emitted stubs.
    list(id = "pair_correlation.sparse", fn = "pair_correlation",
         slot = "univariate_pair_correlation",
         args = list(mnames = SPARSE_UNI, r_range = PCF_R, num_permutations = 5,
                     workers = 1, overwrite = TRUE)),
    list(id = "interaction_variable.sparse", fn = "interaction_variable",
         slot = "interaction_variable",
         args = list(mnames = SPARSE_PAIR, r_range = R_RANGE,
                     num_permutations = 5, workers = 1, overwrite = TRUE)),
    # Ripley's K keeps the shape and emits NA rows instead, but 1.4.0 labelled the
    # sparse branch "Estimator" where its own main branch said "Estimater".
    list(id = "ripleys_k.sparse", fn = "ripleys_k", slot = "univariate_Count",
         args = list(mnames = SPARSE_UNI, r_range = R_RANGE,
                     edge_correction = "translation", permute = FALSE,
                     workers = 1, overwrite = TRUE)),
    list(id = "dixons_s.Z", fn = "dixons_s", slot = "Dixon_Z",
         args = list(mnames = MARKERS, num_permutations = 20, workers = 1,
                     overwrite = TRUE)),
    list(id = "dixons_s.C", fn = "dixons_s", slot = "Dixon_C",
         args = list(mnames = MARKERS, num_permutations = 20, workers = 1,
                     overwrite = TRUE)),
    list(id = "marker_freq_diff", fn = "marker_freq_diff",
         slot = "frequency_difference",
         args = list(classifier = "Classifier.Label", ref_level = "Tumor",
                     diff_level = "Stroma", mnames = MARKERS, overwrite = TRUE))
  )
)

#' Translate one call's arguments for the side being captured
#'
#' Driven off `formals()` rather than a hand-maintained list of what each version
#' accepts, so an argument that exists on one side and not the other is dropped
#' automatically instead of erroring. The one genuine rename is spelled out.
adapt <- function(fn_name, args, side) {
  if (identical(side, "old") && fn_name %in% c("NN_G", "bi_NN_G")) {
    names(args)[names(args) == "keep_permutation_distribution"] <- "keep_perm_dis"
  }
  f <- get(fn_name, mode = "function")
  accepted <- names(formals(f))
  if (!"..." %in% accepted) args <- args[names(args) %in% accepted]
  args
}

# ---- run ------------------------------------------------------------------
results <- list()
for (entry in CALLS) {
  mif <- build_mif()
  args <- adapt(entry$fn, entry$args, SIDE)
  set.seed(SEED)
  # Several old entries are EXPECTED to fail (NN_G km was malformed, dixons_s with
  # overwrite = FALSE always errored). Record the condition instead of aborting --
  # that a side could not produce a value is part of the record.
  val <- tryCatch({
    out <- suppressWarnings(suppressMessages(
      do.call(get(entry$fn, mode = "function"), c(list(mif), args))))
    got <- out$derived[[entry$slot]]
    if (is.null(got)) {
      structure(list(reason = paste0("slot \"", entry$slot, "\" absent; present: ",
                                     paste(names(out$derived), collapse = ", "))),
                class = "parity_absent")
    } else {
      as.data.frame(got)
    }
  }, error = function(e) {
    structure(list(reason = conditionMessage(e)), class = "parity_error")
  })
  results[[entry$id]] <- val
  message("  ", format(entry$id, width = 30),
          if (inherits(val, "parity_error")) paste0("ERROR: ", substr(val$reason, 1, 60))
          else if (inherits(val, "parity_absent")) "ABSENT"
          else paste0(nrow(val), " x ", ncol(val)))
}

# subset_mif and plot_immunoflo are not metric slots, so they are captured apart.
results[["subset_mif"]] <- tryCatch({
  a <- adapt("subset_mif", list(classifier = "Classifier.Label", level = "Tumor",
                                markers = MARKERS), SIDE)
  s <- suppressWarnings(suppressMessages(do.call(subset_mif, c(list(build_mif()), a))))
  list(sample = as.data.frame(s$sample),
       spatial_rows = vapply(s$spatial, nrow, numeric(1)))
}, error = function(e) structure(list(reason = conditionMessage(e)),
                                 class = "parity_error"))

# A serialised ggplot carries its whole environment, so only structure is kept.
results[["plot_immunoflo"]] <- tryCatch({
  p <- suppressWarnings(suppressMessages(plot_immunoflo(
    build_mif(), plot_title = "deidentified_sample", mnames = MARKERS,
    cell_type = "Classifier.Label")))
  g <- p$derived$spatial_plots[[1]]
  list(n_plots = length(p$derived$spatial_plots),
       names = names(p$derived$spatial_plots),
       nrow_data = nrow(g$data),
       mapping = sort(vapply(g$mapping, rlang::as_label, character(1))),
       n_layers = length(g$layers),
       shape_levels = tryCatch(
         length(unique(ggplot2::layer_data(g, 1)$shape)), error = function(e) NA_integer_))
}, error = function(e) structure(list(reason = conditionMessage(e)),
                                 class = "parity_error"))

# ---- provenance -----------------------------------------------------------
# spatstat versions are not optional: "agrees with spatstat" IS the claim, so a
# future failure after a spatstat upgrade must be attributable.
pkg_ver <- function(p) tryCatch(as.character(utils::packageVersion(p)),
                                error = function(e) NA_character_)
git <- function(...) tryCatch(
  system2("git", c("-C", PKG, ...), stdout = TRUE, stderr = FALSE),
  error = function(e) NA_character_)

attr(results, "side") <- SIDE
attr(results, "provenance") <- list(
  spatialTIME_version = got_version,
  pkg_path   = loaded_from,
  git_sha    = git("rev-parse", "HEAD"),
  git_dirty  = git("status", "--porcelain"),
  r_version  = R.version.string,
  platform   = R.version$platform,
  packages   = vapply(c("spatstat.explore", "spatstat.geom", "spatstat.univar",
                        "dixon", "arrow", "dplyr", "spatstat.random"),
                      pkg_ver, character(1)),
  r_md5      = tools::md5sum(list.files(file.path(PKG, "R"), "[.]R$",
                                        full.names = TRUE)),
  captured   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
)
attr(results, "design") <- list(
  samples = SAMPLES, markers = MARKERS, pairs = PAIRS,
  sparse_uni = SPARSE_UNI, sparse_pair = SPARSE_PAIR,
  r_range = R_RANGE, pcf_r = PCF_R, seed = SEED, workers = 1)
# The call table travels with the results so the committed test can REPLAY exactly
# what was captured instead of keeping its own copy of the table, which would drift.
# Everything in it is a string or a plain list, so it serialises cleanly.
attr(results, "calls") <- CALLS

dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
saveRDS(results, OUT)
message("capture: wrote ", OUT, " (", length(results), " entries)")
