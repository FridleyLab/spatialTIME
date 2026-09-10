#' Bivariate Ripley's K
#'
#' @param mif mIF object with spatial data frames, clinical, and per-sample summary information
#' @param mnames vector of column names for phenotypes, or a two-column data frame
#'   of specific anchor/counted marker combinations to run
#' @param r_range vector range of radii at which to calculate co-localization *K*
#' @param edge_correction edge correction method: one of "translation",
#'   "isotropic", "border" or "none"
#' @param num_permutations integer number of permutations used to estimate CSR.
#'   Ignored when `permute = FALSE`.
#' @param permute whether to estimate CSR by permutation (`TRUE`) or to use the
#'   exact closed-form CSR estimate (`FALSE`, the default)
#' @param keep_permutation_distribution boolean; keep each permutation's result or
#'   average them to one row per marker pair and radius
#' @param overwrite boolean; replace an existing `bivariate_Count` slot rather than
#'   appending it as a new `Run`
#' @param workers integer number of CPU workers to use
#' @param xloc,yloc the x and y columns giving cell centres. If left `NULL`,
#'   `XMin`, `XMax`, `YMin` and `YMax` must be present.
#' @param big cell count above which per-pair edge weights are computed in chunks
#'   to bound peak memory. Memory and speed only -- results are identical either
#'   way and the requested `edge_correction` is always honoured.
#'
#' @return mif object with bivariate Ripley's K calculated
#'
#' @description
#' `bi_ripleys_k()` takes a `mIF` object plus marker names and a range of radii,
#' and measures bivariate clustering (co-localization) between each ordered pair of
#' markers. Cells positive for both markers of a pair are excluded from that pair,
#' so the anchor and counted sets are always disjoint.
#'
#' Either estimate CSR by permutation (`permute = TRUE`) or use the exact CSR
#' estimate (`permute = FALSE`). The exact estimate is the univariate K of *all*
#' cells in the sample, which is the closed form for the expected cross-K under
#' random labelling -- verified against a 500-permutation Monte Carlo to within
#' Monte Carlo error.
#'
#' @section Whole slide images:
#' As of 2.0.0 this function handles whole-slide-scale data directly and
#' `bi_ripleys_k_WSI()` is gone. Only the cell pairs closer than `max(r_range)`
#' are ever materialised, so memory scales with the number of nearby pairs rather
#' than with n squared -- at 200,000 cells roughly 75 MB instead of 319 GB. The
#' `big` argument bounds peak memory further by chunking, and unlike the old
#' `big`/`nlarge` arguments it never changes the statistic. In particular the
#' requested edge correction is no longer silently replaced with `"none"` above a
#' cell-count threshold.
#'
#' @section Accuracy:
#' Values agree with [spatstat.explore::Kcross()] to floating-point precision. The
#' observation window is the convex hull of **every** cell in the sample and is
#' held fixed across all marker pairs and permutations.
#'
#' @export
#'
#' @examples
#' x <- spatialTIME::create_mif(clinical_data = spatialTIME::example_clinical %>%
#'                                dplyr::mutate(deidentified_id = as.character(deidentified_id)),
#'                              sample_data = spatialTIME::example_summary %>%
#'                                dplyr::mutate(deidentified_id = as.character(deidentified_id)),
#'                              spatial_list = spatialTIME::example_spatial[1],
#'                              patient_id = "deidentified_id",
#'                              sample_id = "deidentified_sample")
#' x2 = bi_ripleys_k(mif = x,
#'                   mnames = c("CD3..Opal.570..Positive", "CD8..Opal.520..Positive"),
#'                   r_range = seq(0, 100, 10),
#'                   edge_correction = "translation",
#'                   permute = FALSE,
#'                   workers = 1)
bi_ripleys_k = function(mif,
                        mnames,
                        r_range = 0:100,
                        edge_correction = "translation",
                        num_permutations = 50,
                        permute = FALSE,
                        keep_permutation_distribution = FALSE,
                        overwrite = FALSE,
                        workers = 1,
                        xloc = NULL,
                        yloc = NULL,
                        big = 10000){
  if(!inherits(mif, "mif")){
    stop("Please use a mIF object for `mif`, created with `create_mif()`.")
  }
  if(!inherits(mnames, "character") && !inherits(mnames, "data.frame")){
    stop("Please use either a character vector or a data frame of marker combinations for `mnames`.")
  }
  if(keep_permutation_distribution && !permute){
    stop("Conflicting `permute` and `keep_permutation_distribution` parameters.\n",
         "\tTo keep a permutation distribution, set `permute = TRUE`.")
  }
  #r must contain 0 so that the curve starts at the origin (needed for AUC)
  if(!(0 %in% r_range)){
    r_range = sort(c(0, r_range))
  }
  edge_correction = match_edge_correction(edge_correction)

  m_combos = marker_combinations(mnames)
  all_markers = as.character(unique(unlist(mnames)))

  #Seeds drawn in the parent so results depend only on the user's set.seed() and
  #not on `workers`. See the note in ripleys_k().
  seeds = sample.int(.Machine$integer.max, length(mif$spatial))

  out = parallel::mclapply(seq_along(mif$spatial), function(sample_i){
    set.seed(seeds[[sample_i]])
    spat = mif$spatial[[sample_i]]
    spat = add_cell_centres(spat, xloc, yloc)
    label = as.character(spat[[mif$sample_id]][1])

    #Window and area from EVERY cell in the sample, fixed across all marker pairs
    #and permutations.
    win = spatstat.geom::convexhull.xy(spat$xloc, spat$yloc)
    pp  = spatstat.geom::ppp(spat$xloc, spat$yloc, window = win, check = FALSE)
    n   = spatstat.geom::npoints(pp)

    pairs = k_pairs(pp, r_range, edge_correction,
                    block = if(n > big) big else Inf)

    theo  = pi * r_range^2
    #Exact CSR: the univariate K of all cells is the closed form for E[Kcross]
    #under random labelling.
    exact = if(permute) rep(NA_real_, length(r_range)) else k_from_pairs(pairs, rep(TRUE, n))

    res = lapply(seq_len(nrow(m_combos)), function(combo){
      anchor  = as.character(m_combos$anchor[combo])
      counted = as.character(m_combos$counted[combo])

      pos_a = !is.na(spat[[anchor]])  & spat[[anchor]]  == 1
      pos_c = !is.na(spat[[counted]]) & spat[[counted]] == 1
      #Cells positive for both markers belong to neither set.
      keep_i = pos_a & !pos_c
      keep_j = pos_c & !pos_a
      n_i = sum(keep_i); n_j = sum(keep_j)

      if(n_i < 2 || n_j < 2){
        return(bi_k_result_frame(label, anchor, counted, r_range, theo,
                                 observed = NA_real_, permuted = NA_real_, exact = NA_real_,
                                 iter = if(permute) as.character(seq_len(num_permutations)) else "Estimate",
                                 sample_id = mif$sample_id, larger = NA_integer_))
      }

      observed = k_from_pairs(pairs, keep_i, keep_j, n_i, n_j, univariate = FALSE)

      if(!permute){
        return(bi_k_result_frame(label, anchor, counted, r_range, theo, observed,
                                 permuted = NA_real_, exact = exact, iter = "Estimate",
                                 sample_id = mif$sample_id, larger = NA_integer_))
      }

      #Random labelling: draw n_i + n_j cells from all cells in the sample and
      #split them into anchor and counted, then re-mask the same pair list.
      permuted = vapply(seq_len(num_permutations), function(p){
        s = sample.int(n, n_i + n_j)
        ki = logical(n); ki[s[seq_len(n_i)]] = TRUE
        kj = logical(n); kj[s[n_i + seq_len(n_j)]] = TRUE
        k_from_pairs(pairs, ki, kj, n_i, n_j, univariate = FALSE)
      }, numeric(length(r_range)))

      larger = rowSums(permuted > observed, na.rm = TRUE)

      if(keep_permutation_distribution){
        bi_k_result_frame(label, anchor, counted, r_range, theo, observed,
                          permuted = as.vector(permuted), exact = NA_real_,
                          iter = as.character(seq_len(num_permutations)),
                          sample_id = mif$sample_id, larger = larger)
      } else {
        bi_k_result_frame(label, anchor, counted, r_range, theo, observed,
                          permuted = rowMeans(permuted, na.rm = TRUE), exact = NA_real_,
                          iter = "Permuted", sample_id = mif$sample_id, larger = larger)
      }
    })

    dplyr::bind_rows(res)
  }, mc.cores = workers, mc.preschedule = FALSE) %>%
    do.call(dplyr::bind_rows, .) %>%
    dplyr::mutate(`Degree of Clustering Theoretical` = `Observed K` - `Theoretical CSR`,
                  `Degree of Clustering Permutation` = `Observed K` - `Permuted CSR`,
                  `Degree of Clustering Exact`       = `Observed K` - `Exact CSR`)

  write_derived(mif, "bivariate_Count", out, overwrite)
}


#' Expand marker input into an anchor/counted table
#'
#' Accepts either a character vector (all ordered pairs of distinct markers) or a
#' two-column data frame of specific pairs. Shared by the bivariate metrics so
#' they all interpret `mnames` identically.
#'
#' @keywords internal
#' @noRd
marker_combinations <- function(mnames) {
  if (inherits(mnames, "data.frame")) {
    if (ncol(mnames) < 2L) {
      stop("A data frame passed to `mnames` needs two columns: anchor and counted.",
           call. = FALSE)
    }
    out <- data.frame(anchor  = as.character(mnames[[1]]),
                      counted = as.character(mnames[[2]]),
                      stringsAsFactors = FALSE)
  } else {
    if (length(mnames) < 2L) {
      stop("At least two markers are needed for a bivariate measure.", call. = FALSE)
    }
    out <- expand.grid(anchor = mnames, counted = mnames, stringsAsFactors = FALSE)
    out <- out[out$anchor != out$counted, , drop = FALSE]
  }
  out <- out[out$anchor != out$counted, , drop = FALSE]
  if (!nrow(out)) {
    stop("No valid anchor/counted marker pairs were found in `mnames`.", call. = FALSE)
  }
  rownames(out) <- NULL
  out
}


#' Assemble one marker pair's bivariate Ripley's K results
#' @keywords internal
#' @noRd
bi_k_result_frame <- function(label, anchor, counted, r_range, theo,
                              observed, permuted, exact, iter, sample_id, larger) {
  d <- data.frame(
    Label               = label,
    Anchor              = anchor,
    Counted             = counted,
    iter                = rep(iter, each = length(r_range)),
    r                   = r_range,
    `Theoretical CSR`   = theo,
    `Permuted CSR`      = permuted,
    `Exact CSR`         = exact,
    `Observed K`        = observed,
    check.names = FALSE
  )
  d[["Permutations Larger than Observed"]] <- larger
  names(d)[names(d) == "Label"] <- sample_id
  d
}
