#' Bivariate Pair Correlation Function
#'
#' @param mif object of class `mif`
#' @param mnames character vector, or a two-column data frame of anchor/counted
#'   marker combinations to run
#' @param r_range numeric vector of radii. If `NULL`, `spatstat` chooses the range.
#' @param num_permutations integer number of permutations used to estimate CSR
#' @param edge_correction edge correction passed to [spatstat.explore::pcfcross()]
#' @param keep_permutation_distribution boolean; keep each permutation's result or
#'   average them to one row per marker pair and radius
#' @param workers integer number of CPU cores used to process samples in parallel
#' @param overwrite boolean; replace an existing `bivariate_pair_correlation` slot
#'   rather than appending it as a new `Run`
#' @param xloc,yloc the x and y columns giving cell centres. If left `NULL`,
#'   `XMin`, `XMax`, `YMin` and `YMax` must be present.
#' @param ... other parameters passed to [spatstat.explore::pcfcross()], plus
#'   support for deprecated argument names (see Details).
#'
#' @description
#' The cross-type pair correlation function: clustering of counted cells at
#' distance r from anchor cells, as a density rather than a cumulative count.
#' Cells positive for both markers of a pair are excluded from that pair.
#'
#' @details
#' `keep_perm_dis` is accepted as a deprecated alias for
#' `keep_permutation_distribution`.
#'
#' @section Output columns:
#' As of 2.0.0 this function returns the same columns as [bi_ripleys_k()], with
#' `Observed g` in place of `Observed K`. `From`/`To` are now `Anchor`/`Counted`,
#' `Theoretical g`/`Permuted g` are now `Theoretical CSR`/`Permuted CSR`,
#' `Degree of Correlation *` is now `Degree of Clustering *`, and
#' `Permuted_larger_than_Observed` is now `Permutations Larger than Observed`.
#'
#' @return `mif` object with the `bivariate_pair_correlation` slot filled
#' @export
bi_pair_correlation = function(mif,
                               mnames,
                               r_range = NULL,
                               num_permutations = 100,
                               edge_correction = "translation",
                               keep_permutation_distribution = FALSE,
                               workers = 1,
                               overwrite = FALSE,
                               xloc = NULL,
                               yloc = NULL,
                               ...){
  dots = list(...)
  dep  = intersect(names(dots), names(deprecated_arg_map("bi_pair_correlation")))
  apply_deprecated_args(dots[dep], "bi_pair_correlation")
  pcf_args = dots[setdiff(names(dots), dep)]

  if(!inherits(mif, "mif")){
    stop("mIF should be of class `mif` created with function `create_mif()`\n",
         "\tTo check use `inherits(mif, 'mif')`")
  }
  if(length(edge_correction) != 1){
    stop("`edge_correction` must be of length 1.")
  }
  if(!inherits(mnames, "character") && !inherits(mnames, "data.frame")){
    stop("`mnames` must either be a character vector or a data.frame.")
  }
  if(!is.null(r_range) && !(0 %in% r_range)){
    r_range = sort(c(0, r_range))
  }

  m_combos = marker_combinations(mnames)
  needed   = unique(c(m_combos$anchor, m_combos$counted))

  seeds = sample.int(.Machine$integer.max, length(mif$spatial))

  out = parallel::mclapply(seq_along(mif$spatial), function(sample_i){
    set.seed(seeds[[sample_i]])
    spat = mif$spatial[[sample_i]]
    #Fixed in 2.0.0: the old guard was
    #`if(FALSE %in% unique(unlist(mnames)) %in% colnames(spat))`, which because
    #`%in%` is left-associative evaluated `(FALSE %in% mnames) %in% colnames(spat)`
    #and so never fired.
    absent = setdiff(needed, colnames(spat))
    if(length(absent)){
      stop("Marker column(s) not found in spatial data: ",
           paste(absent, collapse = ", "), call. = FALSE)
    }
    spat = add_cell_centres(spat, xloc, yloc)
    label = as.character(spat[[mif$sample_id]][1])

    #Window from EVERY cell in the sample, fixed across marker pairs and permutations.
    win = spatstat.geom::convexhull.xy(spat$xloc, spat$yloc)
    pp  = spatstat.geom::ppp(spat$xloc, spat$yloc, window = win, check = FALSE)
    n   = spatstat.geom::npoints(pp)

    pcf_of = function(keep_i, keep_j){
      keep = keep_i | keep_j
      Y = pp[keep]
      spatstat.geom::marks(Y) = factor(ifelse(keep_i[keep], "i", "j"), levels = c("i", "j"))
      as.data.frame(do.call(spatstat.explore::pcfcross,
                            c(list(Y, i = "i", j = "j", r = r_range,
                                   correction = edge_correction), pcf_args)))
    }

    res = lapply(seq_len(nrow(m_combos)), function(combo){
      anchor  = m_combos$anchor[combo]
      counted = m_combos$counted[combo]
      pos_a = !is.na(spat[[anchor]])  & spat[[anchor]]  == 1
      pos_c = !is.na(spat[[counted]]) & spat[[counted]] == 1
      keep_i = pos_a & !pos_c
      keep_j = pos_c & !pos_a
      n_i = sum(keep_i); n_j = sum(keep_j)
      if(n_i < 3 || n_j < 3) return(NULL)

      obs = pcf_of(keep_i, keep_j)
      est = obs[[ncol(obs)]]

      permuted = vapply(seq_len(num_permutations), function(p){
        s  = sample.int(n, n_i + n_j)
        ki = logical(n); ki[s[seq_len(n_i)]] = TRUE
        kj = logical(n); kj[s[n_i + seq_len(n_j)]] = TRUE
        pcf_of(ki, kj)[[ncol(obs)]]
      }, numeric(nrow(obs)))

      larger   = rowSums(permuted > est, na.rm = TRUE)
      keep_all = keep_permutation_distribution

      d = data.frame(
        Label               = label,
        Anchor              = anchor,
        Counted             = counted,
        iter                = rep(if(keep_all) as.character(seq_len(num_permutations)) else "Permuted",
                                 each = nrow(obs)),
        r                   = obs$r,
        `Theoretical CSR`   = obs$theo,
        `Permuted CSR`      = if(keep_all) as.vector(permuted) else rowMeans(permuted, na.rm = TRUE),
        `Exact CSR`         = NA_real_,
        `Observed g`        = est,
        check.names = FALSE
      )
      d[["Permutations Larger than Observed"]] = larger
      names(d)[names(d) == "Label"] = mif$sample_id
      d
    })

    dplyr::bind_rows(res)
  }, mc.cores = workers, mc.preschedule = FALSE) %>%
    do.call(dplyr::bind_rows, .)

  if(!nrow(out)){
    stop("No marker pair had at least 3 cells of each type in any sample, so no ",
         "bivariate pair correlation could be estimated.")
  }

  out = out %>%
    add_degrees_of_clustering("Observed g") %>%
    as_standard_metric(mif$sample_id, "Observed g", bivariate = TRUE)

  write_derived(mif, "bivariate_pair_correlation", out, overwrite)
}
