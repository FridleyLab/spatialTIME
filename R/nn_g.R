#' Univariate Nearest Neighbor G(r)
#'
#' @param mif object of class `mif` created by function `create_mif()`
#' @param mnames character vector of column names within the spatial files
#'   indicating whether a cell is positive for a phenotype
#' @param r_range numeric vector of radii at which to evaluate G(r)
#' @param num_permutations integer number of permutations used to estimate
#'   sample-specific complete spatial randomness (CSR)
#' @param edge_correction edge correction method: one of "rs", "km", "han" or "none"
#' @param keep_permutation_distribution boolean; keep each permutation's result or
#'   average them to one row per marker and radius
#' @param workers integer number of CPU cores used to process samples in parallel
#' @param overwrite boolean; replace an existing `univariate_NN` slot rather than
#'   appending it as a new `Run`
#' @param xloc,yloc the x and y columns giving cell centres. If left `NULL`,
#'   `XMin`, `XMax`, `YMin` and `YMax` must be present.
#' @param ... support for deprecated argument names. `keep_perm_dis` is accepted as
#'   an alias for `keep_permutation_distribution`. Anything else is an error.
#'
#' @description
#' `NN_G()` computes the nearest-neighbour distance distribution function G(r) for
#' each marker: the proportion of marker-positive cells whose nearest
#' marker-positive neighbour lies within r. CSR is estimated by permutation --
#' relabelling the same cell locations -- so the null accounts for the sample's own
#' geometry rather than assuming a homogeneous Poisson process.
#'
#' @section Accuracy:
#' Estimates come from [spatstat.explore::Gest()], so they agree with spatstat
#' exactly. The observation window is the convex hull of **every** cell in the
#' sample and is held fixed across all markers and permutations, which is what
#' makes G comparable between markers within a sample.
#'
#' @section Why there is no exact CSR for G:
#' [ripleys_k()] can skip permutations because the K of all cells is the closed
#' form for the expected K of a random subset of them. No such shortcut exists for
#' G, because G depends on the *intensity* of the point set, not just its geometry:
#' a random subset of the cells has a lower intensity and therefore larger
#' nearest-neighbour distances. Concretely, on one of the example samples
#' `Gest()` over all cells gives 0.52 at a radius where the mean permuted G is
#' 0.13. The `Exact CSR` column is therefore present for schema consistency with
#' the other metrics but is always `NA`; use `Permuted CSR` or `Theoretical CSR`.
#'
#' @return object of class `mif` with a `univariate_NN` table in the `derived` slot
#' @export
#'
#' @examples
#' x <- spatialTIME::create_mif(clinical_data = spatialTIME::example_clinical %>%
#'   dplyr::mutate(deidentified_id = as.character(deidentified_id)),
#'   sample_data = spatialTIME::example_summary %>%
#'   dplyr::mutate(deidentified_id = as.character(deidentified_id)),
#'   spatial_list = spatialTIME::example_spatial[1],
#'   patient_id = "deidentified_id",
#'   sample_id = "deidentified_sample")
#'
#' x2 = NN_G(mif = x, mnames = "CD3..Opal.570..Positive",
#'   r_range = seq(0, 100, 10), num_permutations = 10,
#'   edge_correction = "rs", workers = 1, overwrite = TRUE)
NN_G = function(mif,
                mnames,
                r_range = 0:100,
                num_permutations = 50,
                edge_correction = "rs",
                keep_permutation_distribution = FALSE,
                workers = 1,
                overwrite = FALSE,
                xloc = NULL,
                yloc = NULL,
                ...){
  apply_deprecated_args(list(...), "NN_G")
  if(!inherits(mif, "mif")){
    stop("Please submit a mif created with `create_mif()`.")
  }
  if(!(0 %in% r_range)){
    r_range = sort(c(0, r_range))
  }
  edge_correction = match_g_correction(edge_correction)

  #Seeds drawn in the parent so results depend only on the user's set.seed() and
  #not on `workers`.
  seeds = sample.int(.Machine$integer.max, length(mif$spatial))

  out = parallel::mclapply(seq_along(mif$spatial), function(sample_i){
    set.seed(seeds[[sample_i]])
    spat = mif$spatial[[sample_i]]
    spat = add_cell_centres(spat, xloc, yloc)
    label = as.character(spat[[mif$sample_id]][1])

    #Window from EVERY cell in the sample, fixed across markers and permutations.
    win = spatstat.geom::convexhull.xy(spat$xloc, spat$yloc)
    pp  = spatstat.geom::ppp(spat$xloc, spat$yloc, window = win, check = FALSE)
    n   = spatstat.geom::npoints(pp)

    res = lapply(mnames, function(marker){
      keep  = !is.na(spat[[marker]]) & spat[[marker]] == 1
      n_pos = sum(keep)

      if(n_pos < 3){
        return(g_result_frame(label, marker, r_range,
                              theo = NA_real_, observed = NA_real_, permuted = NA_real_,
                              iter = if(keep_permutation_distribution)
                                       as.character(seq_len(num_permutations)) else "Permuted",
                              sample_id = mif$sample_id, larger = NA_integer_))
      }

      obs = g_univariate(pp, keep, r_range, edge_correction)

      permuted = vapply(seq_len(num_permutations), function(p){
        keep_p = logical(n)
        keep_p[sample.int(n, n_pos)] = TRUE
        g_univariate(pp, keep_p, r_range, edge_correction)$est
      }, numeric(length(r_range)))

      larger = rowSums(permuted > obs$est, na.rm = TRUE)

      if(keep_permutation_distribution){
        g_result_frame(label, marker, r_range, obs$theo, obs$est,
                       permuted = as.vector(permuted),
                       iter = as.character(seq_len(num_permutations)),
                       sample_id = mif$sample_id, larger = larger)
      } else {
        g_result_frame(label, marker, r_range, obs$theo, obs$est,
                       permuted = rowMeans(permuted, na.rm = TRUE),
                       iter = "Permuted", sample_id = mif$sample_id, larger = larger)
      }
    })

    dplyr::bind_rows(res)
  }, mc.cores = workers, mc.preschedule = FALSE) %>%
    do.call(dplyr::bind_rows, .) %>%
    add_degrees_of_clustering("Observed G") %>%
    as_standard_metric(mif$sample_id, "Observed G", bivariate = FALSE)

  write_derived(mif, "univariate_NN", out, overwrite)
}


#' Assemble one marker's univariate G results
#'
#' `Exact CSR` is present and always NA -- see the "Why there is no exact CSR for
#' G" section of [NN_G()]. It exists so that every metric in the package shares one
#' output schema.
#'
#' @keywords internal
#' @noRd
g_result_frame <- function(label, marker, r_range, theo, observed, permuted,
                           iter, sample_id, larger) {
  d <- data.frame(
    Label             = label,
    Marker            = marker,
    iter              = rep(iter, each = length(r_range)),
    r                 = r_range,
    `Theoretical CSR` = theo,
    `Permuted CSR`    = permuted,
    `Exact CSR`       = NA_real_,
    `Observed G`      = observed,
    check.names = FALSE
  )
  d[["Permutations Larger than Observed"]] <- larger
  names(d)[names(d) == "Label"] <- sample_id
  d
}
