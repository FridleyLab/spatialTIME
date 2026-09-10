#' Bivariate Nearest Neighbor G(r)
#'
#' @param mif object of class `mif` created by function `create_mif()`
#' @param mnames character vector of column names within the spatial files
#'   indicating whether a cell is positive for a phenotype, or a two-column data
#'   frame of specific anchor/counted combinations
#' @param r_range numeric vector of radii at which to evaluate G(r)
#' @param num_permutations integer number of permutations used to estimate
#'   sample-specific complete spatial randomness (CSR)
#' @param edge_correction edge correction method: one of "rs", "km", "han" or "none"
#' @param keep_permutation_distribution boolean; keep each permutation's result or
#'   average them to one row per marker pair and radius
#' @param workers integer number of CPU cores used to process samples in parallel
#' @param overwrite boolean; replace an existing `bivariate_NN` slot rather than
#'   appending it as a new `Run`
#' @param xloc,yloc the x and y columns giving cell centres. If left `NULL`,
#'   `XMin`, `XMax`, `YMin` and `YMax` must be present.
#' @param ... support for deprecated argument names. `keep_perm_dis` is accepted as
#'   an alias for `keep_permutation_distribution`. Anything else is an error.
#'
#' @description
#' `bi_NN_G()` computes the cross-type nearest-neighbour distribution function:
#' for each ordered pair of markers, the proportion of anchor-positive cells whose
#' nearest *counted*-positive neighbour lies within r. Cells positive for both
#' markers of a pair are excluded from that pair, so the anchor and counted sets
#' are always disjoint.
#'
#' @section Accuracy:
#' Estimates come from [spatstat.explore::Gcross()]. Before 2.0.0 this function
#' hand-rolled the `rs` and `han` estimators on a full `as.matrix(dist(...))`,
#' needing memory proportional to the square of the cell count. Those hand-rolled
#' results agreed with `Gcross()` exactly on well-populated markers, but returned
#' `NaN` where `Gcross()` correctly returns `0` for sparse markers, and the `rs`
#' branch referenced an undefined variable. Delegating fixes the sparse-marker
#' case and drops memory to O(n).
#'
#' The observation window is the convex hull of **every** cell in the sample and is
#' held fixed across all marker pairs and permutations.
#'
#' @section Why there is no exact CSR for G:
#' See the corresponding section of [NN_G()]. The `Exact CSR` column exists for
#' schema consistency with [ripleys_k()] and is always `NA`.
#'
#' @return object of class `mif` with a `bivariate_NN` table in the `derived` slot
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
#' x2 = bi_NN_G(mif = x,
#'       mnames = c("CD3..Opal.570..Positive", "CD8..Opal.520..Positive"),
#'       r_range = seq(0, 100, 10), num_permutations = 10,
#'       edge_correction = "rs", workers = 1, overwrite = TRUE)
bi_NN_G = function(mif,
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
  apply_deprecated_args(list(...), "bi_NN_G")
  if(!inherits(mif, "mif")){
    stop("Please submit a mif created with `create_mif()`.")
  }
  if(!inherits(mnames, "character") && !inherits(mnames, "data.frame")){
    stop("Please use either a character vector or a data frame of marker combinations for `mnames`.")
  }
  if(inherits(mnames, "character") && length(mnames) < 2){
    stop("Please use the univariate `NN_G()` for a single marker.")
  }
  if(!(0 %in% r_range)){
    r_range = sort(c(0, r_range))
  }
  edge_correction = match_g_correction(edge_correction)

  m_combos = marker_combinations(mnames)

  #Seeds drawn in the parent so results depend only on the user's set.seed() and
  #not on `workers`.
  seeds = sample.int(.Machine$integer.max, length(mif$spatial))

  out = parallel::mclapply(seq_along(mif$spatial), function(sample_i){
    set.seed(seeds[[sample_i]])
    spat = mif$spatial[[sample_i]]
    spat = add_cell_centres(spat, xloc, yloc)
    label = as.character(spat[[mif$sample_id]][1])

    #Window from EVERY cell in the sample, fixed across marker pairs and permutations.
    win = spatstat.geom::convexhull.xy(spat$xloc, spat$yloc)
    pp  = spatstat.geom::ppp(spat$xloc, spat$yloc, window = win, check = FALSE)
    n   = spatstat.geom::npoints(pp)

    res = lapply(seq_len(nrow(m_combos)), function(combo){
      anchor  = as.character(m_combos$anchor[combo])
      counted = as.character(m_combos$counted[combo])

      pos_a = !is.na(spat[[anchor]])  & spat[[anchor]]  == 1
      pos_c = !is.na(spat[[counted]]) & spat[[counted]] == 1
      #Cells positive for both markers belong to neither set.
      keep_i = pos_a & !pos_c
      keep_j = pos_c & !pos_a
      n_i = sum(keep_i); n_j = sum(keep_j)

      if(n_i < 3 || n_j < 3){
        return(bi_g_result_frame(label, anchor, counted, r_range,
                                 theo = NA_real_, observed = NA_real_, permuted = NA_real_,
                                 iter = if(keep_permutation_distribution)
                                          as.character(seq_len(num_permutations)) else "Permuted",
                                 sample_id = mif$sample_id, larger = NA_integer_))
      }

      obs = g_bivariate(pp, keep_i, keep_j, r_range, edge_correction)

      permuted = vapply(seq_len(num_permutations), function(p){
        s  = sample.int(n, n_i + n_j)
        ki = logical(n); ki[s[seq_len(n_i)]] = TRUE
        kj = logical(n); kj[s[n_i + seq_len(n_j)]] = TRUE
        g_bivariate(pp, ki, kj, r_range, edge_correction)$est
      }, numeric(length(r_range)))

      larger = rowSums(permuted > obs$est, na.rm = TRUE)

      if(keep_permutation_distribution){
        bi_g_result_frame(label, anchor, counted, r_range, obs$theo, obs$est,
                          permuted = as.vector(permuted),
                          iter = as.character(seq_len(num_permutations)),
                          sample_id = mif$sample_id, larger = larger)
      } else {
        bi_g_result_frame(label, anchor, counted, r_range, obs$theo, obs$est,
                          permuted = rowMeans(permuted, na.rm = TRUE),
                          iter = "Permuted", sample_id = mif$sample_id, larger = larger)
      }
    })

    dplyr::bind_rows(res)
  }, mc.cores = workers, mc.preschedule = FALSE) %>%
    do.call(dplyr::bind_rows, .) %>%
    add_degrees_of_clustering("Observed G") %>%
    as_standard_metric(mif$sample_id, "Observed G", bivariate = TRUE)

  write_derived(mif, "bivariate_NN", out, overwrite)
}


#' Assemble one marker pair's bivariate G results
#' @keywords internal
#' @noRd
bi_g_result_frame <- function(label, anchor, counted, r_range, theo, observed,
                              permuted, iter, sample_id, larger) {
  d <- data.frame(
    Label             = label,
    Anchor            = anchor,
    Counted           = counted,
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
