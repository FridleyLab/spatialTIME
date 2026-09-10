#' Univariate Pair Correlation Function
#'
#' @param mif object of class `mif`
#' @param mnames character vector of marker names
#' @param r_range numeric vector including 0. If `NULL`, `spatstat` chooses the range.
#' @param num_permutations integer number of permutations used to estimate CSR
#' @param edge_correction edge correction passed to [spatstat.explore::pcf()]
#' @param keep_permutation_distribution boolean; keep each permutation's result or
#'   average them to one row per marker and radius
#' @param workers integer number of CPU cores used to process samples in parallel
#' @param overwrite boolean; replace an existing `univariate_pair_correlation` slot
#'   rather than appending it as a new `Run`
#' @param xloc,yloc the x and y columns giving cell centres. If left `NULL`,
#'   `XMin`, `XMax`, `YMin` and `YMax` must be present.
#' @param ... other parameters passed to [spatstat.explore::pcf()], plus support
#'   for deprecated argument names (see Details).
#'
#' @description
#' The pair correlation function g(r) is the derivative of Ripley's K, so it
#' measures clustering *at* a radius rather than cumulatively up to it. It is
#' correspondingly slower to compute and noisier at small r.
#'
#' `xloc` and `yloc`, if `NULL`, are taken as the midpoints of `XMin`/`XMax` and
#' `YMin`/`YMax`.
#'
#' @details
#' `keep_perm_dis` is accepted as a deprecated alias for
#' `keep_permutation_distribution`.
#'
#' @section Output columns:
#' As of 2.0.0 this function returns the same columns as [ripleys_k()], with
#' `Observed g` in place of `Observed K`. `Theoretical g` and `Permuted g` are now
#' `Theoretical CSR` and `Permuted CSR`, `Degree of Correlation *` is now
#' `Degree of Clustering *`, and `Permuted_larger_than_Observed` is now
#' `Permutations Larger than Observed`. `Exact CSR` is present but always `NA`:
#' the closed-form CSR shortcut available to [ripleys_k()] has not been
#' implemented for the pair correlation function yet.
#'
#' @return `mif` object with the `univariate_pair_correlation` derived slot filled
#'   or appended to
#' @export
pair_correlation = function(mif,
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
  dep  = intersect(names(dots), names(deprecated_arg_map("pair_correlation")))
  apply_deprecated_args(dots[dep], "pair_correlation")
  pcf_args = dots[setdiff(names(dots), dep)]

  if(!inherits(mif, "mif")){
    stop("mIF should be of class `mif` created with function `create_mif()`\n",
         "\tTo check use `inherits(mif, 'mif')`")
  }
  if(length(edge_correction) != 1){
    stop("`edge_correction` must be of length 1.")
  }
  if(!is.null(r_range) && !(0 %in% r_range)){
    r_range = sort(c(0, r_range))
  }

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

    pcf_of = function(keep){
      as.data.frame(do.call(spatstat.explore::pcf,
                            c(list(pp[keep], r = r_range, correction = edge_correction),
                              pcf_args)))
    }

    res = lapply(mnames, function(marker){
      keep  = !is.na(spat[[marker]]) & spat[[marker]] == 1
      n_pos = sum(keep)
      if(n_pos < 3) return(NULL)

      obs = pcf_of(keep)
      est = obs[[ncol(obs)]]

      permuted = vapply(seq_len(num_permutations), function(p){
        keep_p = logical(n)
        keep_p[sample.int(n, n_pos)] = TRUE
        pcf_of(keep_p)[[ncol(obs)]]
      }, numeric(nrow(obs)))

      larger = rowSums(permuted > est, na.rm = TRUE)
      keep_all = keep_permutation_distribution

      d = data.frame(
        Label               = label,
        Marker              = marker,
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
    stop("No marker had at least 3 positive cells in any sample, so no pair ",
         "correlation could be estimated.")
  }

  out = out %>%
    add_degrees_of_clustering("Observed g") %>%
    as_standard_metric(mif$sample_id, "Observed g", bivariate = FALSE)

  write_derived(mif, "univariate_pair_correlation", out, overwrite)
}
