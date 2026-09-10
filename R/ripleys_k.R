#' Calculate Ripley's K
#'
#' @param mif object of class `mif` created with `create_mif`
#' @param mnames cell phenotype markers to calculate Ripley's K for
#' @param r_range radius range (including 0)
#' @param num_permutations number of permutations to use to estimate CSR. Ignored
#'   when `permute = FALSE`.
#' @param edge_correction edge correction method: one of "translation",
#'   "isotropic", "border" or "none". Unlike previous versions this is never
#'   silently downgraded for large samples.
#' @param permute whether to estimate CSR by permutation (`TRUE`) or to use the
#'   exact closed-form CSR estimate (`FALSE`, the default and much faster)
#' @param keep_permutation_distribution whether to keep each permutation's result
#'   or average them into a single row per marker and radius
#' @param workers number of cores to use for calculations
#' @param overwrite whether to overwrite the `univariate_Count` slot within `mif$derived`
#' @param xloc,yloc columns giving the cell centre. If left `NULL`, `XMin`, `XMax`,
#'   `YMin` and `YMax` must be present and the centre is their midpoint.
#' @param big cell count above which the per-pair edge-correction weights are
#'   computed in chunks to bound peak memory. This affects memory and speed only:
#'   results are identical either way, and the requested `edge_correction` is
#'   always honoured.
#'
#' @description
#' `ripleys_k()` calculates the empirical Ripley's K for the cell types given in
#' `mnames`. This is useful for exploring the spatial clustering of single cell
#' types on TMA cores, ROI spots, or whole slide images following phenotyping with
#' a program such as HALO.
#'
#' Either estimate CSR by permutation (`permute = TRUE`) or use the exact CSR
#' estimate (`permute = FALSE`). The exact estimate is the K of *all* cells in the
#' sample, which is the closed form for the expected K of a randomly chosen subset
#' of them, so it is both faster and free of Monte Carlo error. Permutations are
#' still useful if you want the full null distribution rather than its mean --
#' run 1000 and treat an observed value outside the 95th percentile as
#' significant.
#'
#' @section Accuracy:
#' Values agree with [spatstat.explore::Kest()] to floating-point precision. The
#' observation window is the convex hull of **every** cell in the sample, and it
#' is held fixed across all markers and permutations, so K values for different
#' markers within a sample are directly comparable.
#'
#' Large samples are handled by only ever materialising the cell pairs closer than
#' `max(r_range)`, rather than a full n-by-n distance matrix. At 200,000 cells
#' that is the difference between roughly 75 MB and 319 GB. Because of this,
#' `big` no longer changes the statistic -- in versions before 2.0.0 exceeding it
#' silently replaced your `edge_correction` with `"none"`.
#'
#' @return object of class `mif`
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
#' x2 = ripleys_k(mif = x,
#'   mnames = "CD3..Opal.570..Positive",
#'   r_range = seq(0, 100, 10),
#'   edge_correction = "translation",
#'   permute = FALSE,
#'   workers = 1,
#'   overwrite = TRUE)
ripleys_k = function(mif,
                     mnames,
                     r_range = seq(0, 100, 1),
                     num_permutations = 50,
                     edge_correction = "translation",
                     permute = FALSE,
                     keep_permutation_distribution = FALSE,
                     workers = 1,
                     overwrite = FALSE,
                     xloc = NULL,
                     yloc = NULL,
                     big = 10000){

  if(keep_permutation_distribution && !permute){
    stop("Conflicting `permute` and `keep_permutation_distribution` parameters.\n",
         "\tTo keep a permutation distribution, set `permute = TRUE`.")
  }
  if(!inherits(mif, "mif")){
    stop("mIF should be of class `mif` created with function `create_mif()`\n",
         "\tTo check use `inherits(mif, 'mif')`")
  }
  #r must contain 0 so that the curve starts at the origin (needed for AUC)
  if(!(0 %in% r_range)){
    r_range = sort(c(0, r_range))
  }
  edge_correction = match_edge_correction(edge_correction)

  #Draw one seed per sample in the parent process. Every random draw below is
  #derived from these, so results depend only on the user's set.seed() and not on
  #`workers`. Before 2.0.0 permutations were irreproducible because nested
  #mclapply() calls forked with their own RNG streams.
  seeds = sample.int(.Machine$integer.max, length(mif$spatial))

  out = parallel::mclapply(seq_along(mif$spatial), function(sample_i){
    set.seed(seeds[[sample_i]])
    spat = mif$spatial[[sample_i]]

    spat = add_cell_centres(spat, xloc, yloc)
    label = as.character(spat[[mif$sample_id]][1])
    spat = spat %>%
      dplyr::select(dplyr::all_of(c("xloc", "yloc")), dplyr::any_of(mnames))

    #Window and area come from EVERY cell in the sample, never from a marker
    #subset, and stay fixed for all markers and permutations.
    win = spatstat.geom::convexhull.xy(spat$xloc, spat$yloc)
    pp  = spatstat.geom::ppp(spat$xloc, spat$yloc, window = win, check = FALSE)
    n   = spatstat.geom::npoints(pp)

    #One pair list + one set of edge weights for the whole sample.
    pairs = k_pairs(pp, r_range, edge_correction,
                    block = if(n > big) big else Inf)

    theo = pi * r_range^2
    #Exact CSR is the K of all cells: the closed form for E[K] of a random subset.
    exact = if(permute) rep(NA_real_, length(r_range)) else k_from_pairs(pairs, rep(TRUE, n))

    res = lapply(mnames, function(marker){
      keep = !is.na(spat[[marker]]) & spat[[marker]] == 1
      n_pos = sum(keep)

      if(n_pos < 3){
        #Not enough cells to estimate K; emit the shape but not the numbers.
        return(k_result_frame(label, marker, r_range, theo,
                              observed = NA_real_, permuted = NA_real_, exact = NA_real_,
                              iter = if(permute) as.character(seq_len(num_permutations)) else "Estimate",
                              sample_id = mif$sample_id, larger = NA_integer_))
      }

      observed = k_from_pairs(pairs, keep)

      if(!permute){
        return(k_result_frame(label, marker, r_range, theo, observed,
                              permuted = NA_real_, exact = exact, iter = "Estimate",
                              sample_id = mif$sample_id, larger = NA_integer_))
      }

      #Random labelling: re-mask the same pair list, which is exactly equivalent
      #to recomputing Kest on a relabelled pattern but far cheaper.
      permuted = vapply(seq_len(num_permutations), function(p){
        keep_p = logical(n)
        keep_p[sample.int(n, n_pos)] = TRUE
        k_from_pairs(pairs, keep_p)
      }, numeric(length(r_range)))

      larger = rowSums(permuted > observed, na.rm = TRUE)

      if(keep_permutation_distribution){
        k_result_frame(label, marker, r_range, theo, observed,
                       permuted = as.vector(permuted), exact = NA_real_,
                       iter = as.character(seq_len(num_permutations)),
                       sample_id = mif$sample_id, larger = larger)
      } else {
        k_result_frame(label, marker, r_range, theo, observed,
                       permuted = rowMeans(permuted, na.rm = TRUE), exact = NA_real_,
                       iter = "Permuted", sample_id = mif$sample_id, larger = larger)
      }
    })

    dplyr::bind_rows(res)
  }, mc.cores = workers, mc.preschedule = FALSE) %>%
    do.call(dplyr::bind_rows, .) %>%
    dplyr::mutate(`Degree of Clustering Permutation` = `Observed K` - `Permuted CSR`,
                  `Degree of Clustering Theoretical` = `Observed K` - `Theoretical CSR`,
                  `Degree of Clustering Exact`       = `Observed K` - `Exact CSR`)

  write_derived(mif, "univariate_Count", out, overwrite)
}


#' Assemble one marker's Ripley's K results
#'
#' Recycles `observed`/`exact`/`theo` across permutations, so the caller passes
#' one value per radius for those and either one value per radius (summarised) or
#' `length(r) * num_permutations` values (full distribution) for `permuted`.
#'
#' @keywords internal
#' @noRd
k_result_frame <- function(label, marker, r_range, theo, observed, permuted, exact,
                           iter, sample_id, larger) {
  d <- data.frame(
    iter                = rep(iter, each = length(r_range)),
    Label               = label,
    Marker              = marker,
    r                   = r_range,
    `Theoretical CSR`   = theo,
    `Observed K`        = observed,
    `Permuted CSR`      = permuted,
    `Exact CSR`         = exact,
    check.names = FALSE
  )
  d[["Permutations Larger than Observed"]] <- larger
  names(d)[names(d) == "Label"] <- sample_id
  d
}


#' Resolve cell centres from either explicit columns or a bounding box
#'
#' Shared by every metric function so the fallback rule lives in one place.
#' Unlike the pre-2.0.0 `dplyr::rename('xloc' := xloc)` form, this does not
#' depend on data-masking fallback and so cannot break when a spatial file
#' already contains a column literally named `xloc`.
#'
#' @keywords internal
#' @noRd
add_cell_centres <- function(spat, xloc = NULL, yloc = NULL) {
  if (is.null(xloc) != is.null(yloc)) {
    stop("`xloc` and `yloc` must either both be NULL or both name a column.",
         call. = FALSE)
  }
  if (is.null(xloc)) {
    need <- c("XMin", "XMax", "YMin", "YMax")
    missing <- setdiff(need, colnames(spat))
    if (length(missing)) {
      stop("With `xloc`/`yloc` left NULL the spatial data must contain ",
           paste(need, collapse = ", "), ". Missing: ",
           paste(missing, collapse = ", "), ".", call. = FALSE)
    }
    spat$xloc <- (spat$XMax + spat$XMin) / 2
    spat$yloc <- (spat$YMax + spat$YMin) / 2
  } else {
    for (nm in c(xloc, yloc)) {
      if (!nm %in% colnames(spat)) {
        stop("Column \"", nm, "\" not found in the spatial data.", call. = FALSE)
      }
    }
    spat$xloc <- spat[[xloc]]
    spat$yloc <- spat[[yloc]]
  }
  spat
}


#' Write a results table into a derived slot, overwriting or appending a new Run
#'
#' Centralised because each metric function had its own copy of this logic and
#' three of them wrote to a misspelled slot name on the append path.
#'
#' @keywords internal
#' @noRd
write_derived <- function(mif, slot, out, overwrite) {
  existing <- mif$derived[[slot]]
  if (overwrite || is.null(existing) || !nrow(existing)) {
    mif$derived[[slot]] <- dplyr::mutate(out, Run = 1)
  } else {
    mif$derived[[slot]] <- dplyr::bind_rows(
      existing,
      dplyr::mutate(out, Run = max(existing$Run, na.rm = TRUE) + 1)
    )
  }
  mif
}
