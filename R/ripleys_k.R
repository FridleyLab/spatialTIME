#' Calculate Ripley's K
#'
#' @param mif object of class `mif` created with `create_mif`
#' @param mnames cell phenotype markers to calculate Ripley's K for
#' @param r_range radius range (including 0)
#' @param num_permutations number of permutations to use to estimate CSR. Ignored
#'   when `permute = FALSE`.
#' @param edge_correction edge correction method: one of "translation",
#'   "isotropic", "border" or "none". Unlike previous versions this is never
#'   silently downgraded for large samples. All four agree with
#'   [spatstat.explore::Kest()] to floating-point precision; note that `"border"`
#'   reports no `Exact CSR` (see below) and that `"none"` applies no correction at
#'   all and so is biased downward near the window edge.
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
#' @param ... support for deprecated argument names. `keep_perm_dis` is accepted as
#'   an alias for `keep_permutation_distribution`; `method` is accepted and ignored
#'   (it was never used). Anything else is an error.
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
#' `Exact CSR` is filled in either way, so with `permute = TRUE` you can read
#' `Permuted CSR` against it and see directly whether your permutation count was
#' large enough to converge. They agree exactly in expectation.
#'
#' @section The CSR columns, and which to trust:
#' Three references are reported for each radius, and `Degree of Clustering X` is
#' always `Observed K - X CSR`:
#'
#' * `Theoretical CSR` is \eqn{\pi r^2}, the K of a homogeneous Poisson process.
#'   It ignores the shape of the window and the cell density actually present, so
#'   it is the weakest of the three on real tissue.
#' * `Exact CSR` is the K of every cell in the sample. Under random labelling of a
#'   fixed set of cell locations this is *exactly* the expected K of a marker-positive
#'   subset -- not an approximation -- for translation, isotropic and none. It is
#'   the reference to prefer.
#' * `Permuted CSR` is the mean over `num_permutations` random relabellings.
#'
#' `Exact CSR` is `NA` for `edge_correction = "border"`. Border is a reduced-sample
#' estimator whose denominator counts only the cells still further than `r` from
#' the window edge, and that count changes with which cells are marker-positive. So
#' the cancellation that makes the other three exact does not apply, and the K of
#' all cells is *not* the expected K of a subset: measured on 600 cells with 120
#' positive, it sits about 1% below the mean permuted K at larger radii. Rather
#' than report a number that is quietly wrong in a column called "Exact", it is
#' omitted -- as it is for [NN_G()] and [pair_correlation()], for the same reason.
#' Use `permute = TRUE` with border.
#'
#' @section Permutation p-values:
#' `Permutations Larger than Observed` counts the permutations at least as extreme
#' as the observation (`>=`), so ties count toward the count and therefore toward a
#' larger p-value. `Permutation p-value` is the standard Monte Carlo p-value,
#' \eqn{(1 + \#\{perm \ge obs\}) / (B + 1)}, where `B` is the number of
#' permutations that actually produced a value at that radius. The `+1`s mean it
#' can never be exactly 0, which a raw count divided by `B` can be, and which would
#' read as infinite significance rather than "nothing in this sample was more
#' extreme".
#'
#' Both are `NA` at radii where `Observed K` is `NA`.
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
                     big = 10000,
                     ...){
  apply_deprecated_args(list(...), "ripleys_k")

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
  seeds = sample.int(.Machine$integer.max, n_mif_samples(mif))

  out = parallel::mclapply(seq_len(n_mif_samples(mif)), function(sample_i){
    set.seed(seeds[[sample_i]])
    #On a disk-backed mif this reads only the columns named below, for this one
    #sample, inside this worker; on an in-memory mif it is `mif$spatial[[sample_i]]`
    #unchanged. See R/utils-mif-store.R.
    spat = mif_spatial(mif, sample_i,
                       spatial_columns(mif, mnames, xloc, yloc, i = sample_i))

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
    #
    #Computed whether or not `permute` is set. Under random labelling E[K_perm] is
    #EXACTLY this quantity -- sample.int() is sampling without replacement, so a
    #pair is retained with probability m(m-1)/(n(n-1)) while the denominator is
    #m(m-1)/area, and the two cancel to the K of all cells. So with `permute = TRUE`
    #the user can read `Permuted CSR` against `Exact CSR` and see directly whether
    #their permutation count was large enough to converge. Before 2.0.0 this was
    #hard-NA whenever `permute = TRUE`, which hid exactly that comparison. It costs
    #one extra mask over a pair list that has already been built.
    #See tests/testthat/test-csr-null.R, which proves the identity by enumeration.
    #
    #BUT NOT FOR BORDER. The cancellation above needs a denominator that is fixed
    #once m is fixed. Translation, isotropic and none all have one: m(m-1)/area.
    #Border does not -- its denominator counts the SELECTED cells still further
    #than r from the window edge, which varies from one relabelling to the next, so
    #E[numerator/denominator] is not E[numerator]/E[denominator]. Measured on 600
    #cells with m = 120 and 3000 permutations, the K of all cells sits below the
    #mean permuted K by 0.3% at r = 20 rising to 1.0% at r = 80 (z = -0.9 to -7.9),
    #i.e. a real bias and not Monte Carlo noise, growing with r.
    #
    #So border gets NA here rather than a number that is 1% wrong in a column
    #called "Exact CSR". Same precedent as G and the pair correlation function --
    #see the "Why there is no exact CSR for G" section of ?NN_G. Use
    #`permute = TRUE` with border; `Permuted CSR` is still correct for it.
    exact = if(pairs$edge_correction == "border") rep(NA_real_, length(r_range))
            else k_from_pairs(pairs, rep(TRUE, n))

    res = lapply(mnames, function(marker){
      keep = !is.na(spat[[marker]]) & spat[[marker]] == 1
      n_pos = sum(keep)

      if(n_pos < 3){
        #Not enough cells to estimate K; emit the shape but not the numbers.
        return(k_result_frame(label, marker, r_range, theo,
                              observed = NA_real_, permuted = NA_real_, exact = NA_real_,
                              iter = if(permute) as.character(seq_len(num_permutations)) else "Estimate",
                              sample_id = mif$sample_id, larger = NA_integer_,
                              p_value = NA_real_))
      }

      observed = k_from_pairs(pairs, keep)

      if(!permute){
        return(k_result_frame(label, marker, r_range, theo, observed,
                              permuted = NA_real_, exact = exact, iter = "Estimate",
                              sample_id = mif$sample_id, larger = NA_integer_,
                              p_value = NA_real_))
      }

      #Random labelling: re-mask the same pair list, which is exactly equivalent
      #to recomputing Kest on a relabelled pattern but far cheaper.
      permuted = vapply(seq_len(num_permutations), function(p){
        keep_p = logical(n)
        keep_p[sample.int(n, n_pos)] = TRUE
        k_from_pairs(pairs, keep_p)
      }, numeric(length(r_range)))

      ps = permutation_summary(permuted, observed)

      if(keep_permutation_distribution){
        k_result_frame(label, marker, r_range, theo, observed,
                       permuted = as.vector(permuted), exact = exact,
                       iter = as.character(seq_len(num_permutations)),
                       sample_id = mif$sample_id, larger = ps$larger,
                       p_value = ps$p_value)
      } else {
        k_result_frame(label, marker, r_range, theo, observed,
                       permuted = rowMeans(permuted, na.rm = TRUE), exact = exact,
                       iter = "Permuted", sample_id = mif$sample_id,
                       larger = ps$larger, p_value = ps$p_value)
      }
    })

    dplyr::bind_rows(res)
  }, mc.cores = workers, mc.preschedule = FALSE) %>%
    do.call(dplyr::bind_rows, .) %>%
    add_degrees_of_clustering("Observed K") %>%
    as_standard_metric(mif$sample_id, "Observed K", bivariate = FALSE)

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
                           iter, sample_id, larger, p_value) {
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
  d[["Permutation p-value"]] <- p_value
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
