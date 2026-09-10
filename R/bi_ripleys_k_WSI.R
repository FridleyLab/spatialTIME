#' Bivariate Ripley's K for Whole Slide Images (removed)
#'
#' @description
#' `r lifecycle_note_wsi()`
#'
#' @details
#' `bi_ripleys_k_WSI()` existed because [bi_ripleys_k()] used to build an n-by-n
#' distance matrix, which is impossible at whole-slide scale, and so a separate
#' tiled implementation was maintained alongside it. The two functions were
#' otherwise near-duplicates, and the tiled one traded accuracy for memory: above
#' a cell-count threshold it replaced the requested edge correction with `"none"`,
#' its tiled branch hard-coded translation regardless of the `edge_correction`
#' argument, and it silently returned garbage for `edge_correction = "isotropic"`
#' because that branch of its inner loop returned `NULL`.
#'
#' As of 2.0.0 that trade-off is gone. [bi_ripleys_k()] only ever materialises the
#' cell pairs closer than `max(r_range)`, so whole-slide images are handled
#' directly -- at 200,000 cells roughly 75 MB rather than 319 GB -- with the
#' requested edge correction always honoured and results equal to
#' [spatstat.explore::Kcross()] to floating-point precision.
#'
#' Migration is a rename. `big` and `nlarge` have no equivalent because they no
#' longer affect the statistic; [bi_ripleys_k()] takes a `big` argument that bounds
#' peak memory only.
#'
#' \preformatted{
#' # before
#' bi_ripleys_k_WSI(mif, mnames, r_range = 0:100, big = 1000, nlarge = 1000)
#'
#' # after
#' bi_ripleys_k(mif, mnames, r_range = 0:100)
#' }
#'
#' @param mif,mnames,r_range,edge_correction,num_permutations,permute Ignored.
#' @param keep_permutation_distribution,overwrite,workers,big,nlarge,xloc,yloc Ignored.
#' @param ... Ignored.
#'
#' @return Nothing. Always throws an error directing you to [bi_ripleys_k()].
#' @seealso [bi_ripleys_k()]
#' @export
bi_ripleys_k_WSI = function(mif, mnames, r_range = 0:100,
                            edge_correction = "translation",
                            num_permutations = 50, permute = FALSE,
                            keep_permutation_distribution = FALSE,
                            overwrite = FALSE, workers = 1,
                            big = 1000, nlarge = 1000,
                            xloc = NULL, yloc = NULL, ...){
  stop("`bi_ripleys_k_WSI()` was removed in spatialTIME 2.0.0.\n",
       "  Use `bi_ripleys_k()`, which now handles whole slide images natively:\n",
       "  it only materialises cell pairs within `max(r_range)`, so memory scales\n",
       "  with nearby pairs rather than with the square of the cell count, and the\n",
       "  edge correction you ask for is always the one applied.\n\n",
       "  Drop `big`/`nlarge` -- they no longer change the statistic. `bi_ripleys_k()`\n",
       "  has its own `big` argument that bounds peak memory only.\n",
       "  See ?bi_ripleys_k_WSI for the full rationale.",
       call. = FALSE)
}

#' Text for the WSI removal notice
#' @keywords internal
#' @noRd
lifecycle_note_wsi <- function() {
  paste("**Removed in spatialTIME 2.0.0.** Use [bi_ripleys_k()] instead, which now",
        "handles whole slide images directly. Calling this function raises an error.")
}
