#' Bivariate Interaction Variable
#'
#' Single-cell spatial-protein metric introduced by Steinhart et al.,
#' \doi{10.1158/1541-7786.MCR-21-0411}.
#'
#' @param mif object of class `mif`
#' @param mnames a character vector, or a two-column data frame of anchor/counted
#'   markers to assess
#' @param r_range numeric vector of radii at which to evaluate the interaction
#'   variable
#' @param num_permutations integer number of permutations used to derive the
#'   interaction estimate under CSR
#' @param keep_permutation_distribution boolean; keep each permutation's result or
#'   average them to one row per marker pair and radius
#' @param workers integer number of CPU cores used to process samples in parallel
#' @param overwrite boolean; replace an existing `interaction_variable` slot rather
#'   than appending it as a new `Run`
#' @param xloc,yloc the x and y columns giving cell centres. If left `NULL`,
#'   `XMin`, `XMax`, `YMin` and `YMax` must be present.
#' @param ... support for deprecated argument names (see Details).
#'
#' @description
#' For each anchor cell, find the distance to its nearest counted cell; then report
#' the cumulative percentage of cells whose nearest counted neighbour falls within
#' each radius. Cells positive for both markers of a pair are excluded from that
#' pair.
#'
#' @details
#' `keep_perm_dis` is accepted as a deprecated alias for
#' `keep_permutation_distribution`.
#'
#' @section How the percentage is normalised:
#' The numerator counts anchor cells (one nearest-neighbour distance per anchor
#' cell), but the denominator is the **total** number of anchor plus counted cells
#' in the pair. `Observed Interaction` therefore approaches
#' `100 * n_anchor / (n_anchor + n_counted)` rather than 100 as r grows, and its
#' ceiling depends on the relative abundance of the two markers.
#'
#' This is preserved exactly as implemented before 2.0.0 so that existing results
#' remain comparable -- it is documented here rather than changed, because
#' redefining a published metric is not a refactoring decision. If you want a
#' quantity that reaches 100%, divide by the anchor count:
#' `Observed Interaction * (n_anchor + n_counted) / n_anchor`.
#'
#' @section Output columns:
#' As of 2.0.0 this function returns the same columns as [bi_ripleys_k()], with
#' `Observed Interaction` in place of `Observed K`. `From`/`To` are now
#' `Anchor`/`Counted`, `Permuted Interaction` is now `Permuted CSR`, and
#' `Degree of Interaction Permuted` is now `Degree of Clustering Permutation`.
#' `Theoretical CSR` and `Exact CSR` are present but always `NA`: this metric has
#' no closed-form null, which is why it is permutation-only.
#'
#' @return object of class `mif` with the `interaction_variable` derived slot filled
#' @export
interaction_variable = function(mif,
                                mnames,
                                r_range = NULL,
                                num_permutations = 100,
                                keep_permutation_distribution = FALSE,
                                workers = 1,
                                overwrite = FALSE,
                                xloc = NULL,
                                yloc = NULL,
                                ...){
  apply_deprecated_args(list(...), "interaction_variable")

  if(!inherits(mif, "mif")){
    stop("mIF should be of class `mif` created with function `create_mif()`\n",
         "\tTo check use `inherits(mif, 'mif')`")
  }
  if(!inherits(mnames, "character") && !inherits(mnames, "data.frame")){
    stop("`mnames` must either be a character vector or a data.frame.")
  }
  if(is.null(r_range)){
    stop("`r_range` must be supplied for the interaction variable.")
  }
  if(!(0 %in% r_range)){
    r_range = sort(c(0, r_range))
  }
  if(length(r_range) < 2){
    stop("`r_range` needs at least two radii.")
  }

  m_combos = marker_combinations(mnames)
  needed   = unique(c(m_combos$anchor, m_combos$counted))

  seeds = sample.int(.Machine$integer.max, length(mif$spatial))

  out = parallel::mclapply(seq_along(mif$spatial), function(sample_i){
    set.seed(seeds[[sample_i]])
    spat = mif$spatial[[sample_i]]
    #Fixed in 2.0.0: the old guard `FALSE %in% x %in% colnames(spat)` could never
    #fire because `%in%` is left-associative.
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

    #Cumulative percentage of cells whose nearest counted neighbour is within r.
    #Denominator is n_i + n_j -- see the "How the percentage is normalised" section.
    interaction_of = function(keep_i, keep_j){
      nnd = spatstat.geom::nncross(pp[keep_i], pp[keep_j], what = "dist")
      denom = sum(keep_i) + sum(keep_j)
      c(0, cumsum(table(cut(nnd, breaks = r_range, include.lowest = TRUE)))) /
        denom * 100
    }

    res = lapply(seq_len(nrow(m_combos)), function(combo){
      anchor  = m_combos$anchor[combo]
      counted = m_combos$counted[combo]
      pos_a = !is.na(spat[[anchor]])  & spat[[anchor]]  == 1
      pos_c = !is.na(spat[[counted]]) & spat[[counted]] == 1
      keep_i = pos_a & !pos_c
      keep_j = pos_c & !pos_a
      n_i = sum(keep_i); n_j = sum(keep_j)
      if(n_i < 1 || n_j < 1) return(NULL)

      obs = unname(interaction_of(keep_i, keep_j))

      permuted = vapply(seq_len(num_permutations), function(p){
        s  = sample.int(n, n_i + n_j)
        ki = logical(n); ki[s[seq_len(n_i)]] = TRUE
        kj = logical(n); kj[s[n_i + seq_len(n_j)]] = TRUE
        unname(interaction_of(ki, kj))
      }, numeric(length(r_range)))

      larger   = rowSums(permuted > obs, na.rm = TRUE)
      keep_all = keep_permutation_distribution

      d = data.frame(
        Label                  = label,
        Anchor                 = anchor,
        Counted                = counted,
        iter                   = rep(if(keep_all) as.character(seq_len(num_permutations)) else "Permuted",
                                    each = length(r_range)),
        r                      = r_range,
        `Theoretical CSR`      = NA_real_,
        `Permuted CSR`         = if(keep_all) as.vector(permuted) else rowMeans(permuted, na.rm = TRUE),
        `Exact CSR`            = NA_real_,
        `Observed Interaction` = obs,
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
    stop("No marker pair had cells of both types in any sample, so no interaction ",
         "variable could be estimated.")
  }

  out = out %>%
    add_degrees_of_clustering("Observed Interaction") %>%
    as_standard_metric(mif$sample_id, "Observed Interaction", bivariate = TRUE)

  write_derived(mif, "interaction_variable", out, overwrite)
}
