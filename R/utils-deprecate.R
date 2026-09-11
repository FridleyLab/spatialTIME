# Argument-name compatibility for the 2.0.0 unification.
#
# Before 2.0.0 the metric functions disagreed about what to call the same thing:
# keep_perm_dis vs keep_permutation_distribution, big vs nlarge, a `force` flag
# that only existed to let you past a size check, and a documented `method`
# argument in ripleys_k() that the body never read. 2.0.0 settles on one name per
# concept, and every renamed argument keeps working for the 2.x line via `...`,
# with a warning naming its replacement.
#
# Adding `...` to a user-facing function normally means typos get silently
# swallowed. apply_deprecated_args() closes that hole: anything in `...` that is
# not a known former argument of that function is an error, so `workerss = 4`
# fails loudly instead of being ignored.


#' Former argument names, per function
#'
#' `NA_character_` marks an argument that is gone with no replacement -- it is
#' accepted and ignored so old scripts still run, but it no longer does anything.
#'
#' @keywords internal
#' @noRd
deprecated_arg_map <- function(fn) {
  common <- c(keep_perm_dis = "keep_permutation_distribution")
  switch(fn,
    ripleys_k = c(common,
      # Documented as "not used currently" and never referenced in the body.
      method = NA_character_),
    bi_ripleys_k = c(common,
      nlarge = "big",
      # `force` existed only to let you past a hard stop at 10,000 cells that
      # pushed you towards bi_ripleys_k_WSI(). There is no size limit any more.
      force = NA_character_),
    bi_ripleys_k_WSI = c(common, nlarge = "big", force = NA_character_),
    NN_G = common,
    bi_NN_G = common,
    pair_correlation = common,
    bi_pair_correlation = common,
    interaction_variable = common,
    # split_tissue()/plot_tissue_split() are new in 2.1.0 and never had a
    # `keep_perm_dis` argument to be deprecated -- falling through to `common`
    # would silently accept it instead of erroring on the typo/unknown argument.
    split_tissue = character(0),
    plot_tissue_split = character(0),
    common
  )
}


#' Map former argument names onto their replacements
#'
#' Called once at the top of each metric function. Assigns each recognised former
#' argument's value to its new name in the caller's frame, so the rest of the
#' function body only ever sees canonical names.
#'
#' @param dots `list(...)` from the caller.
#' @param fn name of the calling function, used to pick the mapping and to write
#'   useful messages.
#' @param env frame to assign into; defaults to the caller's.
#' @return invisibly `NULL`; called for its side effects.
#' @keywords internal
#' @noRd
apply_deprecated_args <- function(dots, fn, env = parent.frame()) {
  if (!length(dots)) return(invisible(NULL))
  mapping <- deprecated_arg_map(fn)

  nms <- names(dots)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop(sprintf("`%s()` received unnamed arguments in `...`. ", fn),
         "All arguments must be named.", call. = FALSE)
  }

  unknown <- setdiff(nms, names(mapping))
  if (length(unknown)) {
    stop(sprintf("Unknown argument%s passed to `%s()`: %s.\n",
                 if (length(unknown) > 1) "s" else "", fn,
                 paste0("`", unknown, "`", collapse = ", ")),
         "Check for a typo, or see ?", fn, " for the current arguments.",
         call. = FALSE)
  }

  for (old in nms) {
    new <- mapping[[old]]
    if (is.na(new)) {
      warning(sprintf(
        "`%s` is defunct as of spatialTIME 2.0.0 and is being ignored.\n  %s",
        old, defunct_reason(old)), call. = FALSE)
    } else {
      warning(sprintf(
        "`%s` is deprecated as of spatialTIME 2.0.0; use `%s` instead.\n  Forwarding the value for now; `%s` will be removed in 3.0.0.",
        old, new, old), call. = FALSE)
      assign(new, dots[[old]], envir = env)
    }
  }
  invisible(NULL)
}


#' Why a defunct argument no longer does anything
#' @keywords internal
#' @noRd
defunct_reason <- function(old) {
  switch(old,
    method = paste("`ripleys_k()` has always ignored it -- it was documented as",
                   "\"not used currently\" and never read by the function body."),
    force = paste("There is no longer a cell-count limit to force past:",
                  "`bi_ripleys_k()` handles whole slide images directly."),
    sprintf("See ?spatialTIME-deprecated for details on `%s`.", old)
  )
}
