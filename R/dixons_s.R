#' Dixon's S Segregation Statistic
#'
#' @description This function processes the spatial files in the mif object,
#' requiring a column that distinguishes between different groups i.e. tumor and 
#' stroma
#' @param mif An MIF object
#' @param mnames vector of markers corresponding to spatial columns to check Dixon's S between
#' @param num_permutations Numeric value indicating the number of permutations used. 
#'  Default is 1000.
#' @param type a character string for the type that is wanted in the output which can
#' be "Z" for z-statistic results or "C" for Chi-squared statistic results
#' @param workers Integer value for the number of workers to spawn
#' @param overwrite Logical value determining if you want the results to replace the 
#' current output (TRUE) or be to be appended (FALSE).
#' @param xloc a string corresponding to the x coordinates. If null the average of 
#' XMin and XMax will be used 
#' @param yloc a string corresponding to the y coordinates. If null the average of 
#' YMin and YMax will be used 
#' @importFrom magrittr %>%
#' 
#' @return An object of class `mif`. Depending on `type`, one or both of these
#'   tables is written to the `derived` slot.
#'
#'   `Dixon_Z` (`type = "Z"`), one row per ordered marker pair per sample:
#'   \describe{
#'     \item{<sample_id>}{sample identifier}
#'     \item{From, To}{the two cell types being compared, as `dixon::dixon()` labels them}
#'     \item{Anchor, Counted}{which requested marker pair the row belongs to. The
#'       per-marker count columns are pair-specific, because dual positives are
#'       dropped per pair, so without this the counts cannot be interpreted}
#'     \item{Obs.Count, Exp.Count, S, Z, p-val.Z, p-val.Nobs}{Dixon's statistics}
#'     \item{<marker> columns}{cells of each marker retained for that pair}
#'     \item{Simulations, Run}{permutations requested, and which call produced the row}
#'   }
#'
#'   `Dixon_C` (`type = "C"`), three rows per marker pair per sample:
#'   \describe{
#'     \item{<sample_id>}{sample identifier}
#'     \item{Direction}{"Overall segregation", or "From <marker>"}
#'     \item{Anchor, Counted}{as above}
#'     \item{df, Chi-sq, P.asymp, P.rand}{chi-squared segregation test}
#'     \item{Simulations, Run}{as above}
#'   }
#'
#'   Marker pairs with fewer than 3 cells of either type are returned with `NA`
#'   statistics rather than dropped, so the table shape does not depend on cell
#'   counts.
#' @examples 
#' #' #Create mif object
#' library(dplyr)
#' x <- create_mif(clinical_data = example_clinical %>% 
#' mutate(deidentified_id = as.character(deidentified_id)),
#' sample_data = example_summary %>% 
#' mutate(deidentified_id = as.character(deidentified_id)),
#' spatial_list = example_spatial,
#' patient_id = "deidentified_id", 
#' sample_id = "deidentified_sample")
#' 
#' @export

dixons_s = function(mif, mnames, num_permutations = 1000, type = c("Z", "C"),
                    workers = 1, overwrite = FALSE, xloc = NULL, yloc = NULL){
  
  if(!inherits(mif, "mif")){
    stop("mIF should be of class `mif` created with function `create_mif()`")
  }
  if(!inherits(mnames, "character")){
    stop("Provide a vector of marker names in the spatial files")
  }
  if(length(mnames) < 2){
    stop("Dixon's S compares two cell types, so at least 2 markers are needed in `mnames`.")
  }
  #`type` was never validated, so a typo such as type = "z" ran the whole
  #permutation loop and then returned the mif unchanged, with no message.
  type = unique(type)
  if(!length(type) || !all(type %in% c("Z", "C"))){
    stop("`type` must be \"Z\", \"C\" or both.")
  }
  data = mif$spatial
  #filter names of markers because order doesn't matter, still computes both ways
  mnames = mnames %>%
    expand.grid(., .) %>%
    dplyr::filter(Var1 != Var2) %>% 
    dplyr::rowwise() %>%
    dplyr::mutate(Var3 = paste0(sort(c(Var1, Var2)), collapse = ",")) %>%
    dplyr::distinct(Var3, .keep_all = TRUE) %>%
    dplyr::select(1, 2) %>%
    dplyr::ungroup()
  
  #per spatial file
  out = parallel::mclapply(data, function(spat){
    #if locations not provided
    if(is.null(xloc)){
      spat$xloc = (spat$XMin+spat$XMax)/2
    } else {
      spat$xloc = spat[[xloc]]
    }
    if(is.null(yloc)){
      spat$yloc = (spat$YMin+spat$YMax)/2
    } else {
      spat$yloc = spat[[yloc]]
    }
    
    #Serial: `workers` parallelises over samples in the enclosing mclapply. This
    #inner loop used to be an mclapply with no mc.cores, so it forked onto
    #getOption("mc.cores", 2L) regardless of `workers`, and any error inside it
    #surfaced as "$ operator is invalid for atomic vectors" instead of the real
    #message.
    res = lapply(seq_len(nrow(mnames)), function(r){
      markers = as.character(unlist(mnames[r,]))
      missing_markers = setdiff(markers, colnames(spat))
      if(length(missing_markers)){
        stop("Marker column(s) not found in spatial data: ",
             paste(missing_markers, collapse = ", "), call. = FALSE)
      }
      df = spat %>%
        dplyr::select(xloc, yloc, dplyr::all_of(markers)) %>%
        #drop dual positives
        dplyr::filter(!(.data[[markers[1]]] == 1 & .data[[markers[2]]] == 1)) %>%
        tidyr::gather("Marker", "Positive", -xloc, -yloc) %>%
        dplyr::mutate(Marker = factor(Marker, levels = markers)) %>%
        dplyr::filter(Positive == 1)

      #`Marker` is a factor carrying both levels, so table() always has 2 rows --
      #the old `nrow(df_tab) == 0 | nrow(df_tab) == 1` clauses were dead code.
      #Minimum 3 cells per group to estimate anything.
      df_tab = table(df$Marker)
      if(any(df_tab < 3)){
        #Shape must match the success branch exactly, otherwise bind_rows across
        #marker pairs produces a ragged table. Previously this branch invented its
        #own columns, referenced a non-existent `Image Location` (silently a no-op,
        #since mutate(x = NULL) removes rather than adds), and renamed the `df`
        #column to the sample id -- which the caller then overwrote, destroying the
        #degrees of freedom rather than reporting it.
        final_df_z = expand.grid(From = markers, To = markers,
                                 stringsAsFactors = FALSE) %>%
          dplyr::mutate(Obs.Count = NA_real_, Exp.Count = NA_real_,
                        S = NA_real_, Z = NA_real_,
                        `p-val.Z` = NA_real_, `p-val.Nobs` = NA_real_)
        final_df_c = data.frame(df = NA_real_, `Chi-sq` = NA_real_,
                                P.asymp = NA_real_, P.rand = NA_real_,
                                check.names = FALSE)[rep(1, 3), , drop = FALSE]
        rownames(final_df_c) = c("Overall segregation", paste("From", markers))
      } else {
        #dixon::dixon() prints its permutation counter to stdout. Capture it: this
        #package removed its own cat() progress spam in 2.0.0 and should not leak a
        #dependency's either.
        invisible(utils::capture.output(
          dixon_val <- dixon::dixon(df, nsim = num_permutations)
        ))
        final_df_z = dixon_val$tablaZ %>%
          dplyr::rename_with(~ .x %>% gsub(" ", "", .))
        #tablaC's own column names carry padding spaces ("  df ", "  P.rand").
        #Only tablaZ used to be de-spaced, so binding a success row to an
        #early-return row yielded BOTH "  P.rand" and "P.rand".
        final_df_c = dixon_val$tablaC %>%
          dplyr::rename_with(~ .x %>% gsub(" ", "", .))
      }
      #Which pair a row belongs to, on both tables. The per-marker count columns
      #are pair-specific -- the dual-positive filter runs per pair -- so without
      #this the counts are not interpretable.
      final_df_z = final_df_z %>%
        dplyr::mutate(Anchor = markers[1], Counted = markers[2],
                      !!markers[1] := as.numeric(df_tab[1]),
                      !!markers[2] := as.numeric(df_tab[2]))
      final_df_c$Anchor = markers[1]
      final_df_c$Counted = markers[2]
      return(list(tablaZ = final_df_z, tablaC = final_df_c))
    })
    #bring together
    Dixon_Z = lapply(res, function(marks){
      marks$tablaZ %>%
        dplyr::mutate(!!mif$sample_id := spat[[mif$sample_id]][1], .before = 1)
    }) %>%
      do.call(dplyr::bind_rows, .) %>%
      dplyr::mutate(Simulations = num_permutations)
    Dixon_C = lapply(res, function(marks){
      marks$tablaC%>%
        tibble::rownames_to_column("Direction") %>%
        dplyr::mutate(Direction = stringr::str_squish(Direction)) %>%
        dplyr::mutate(!!mif$sample_id := spat[[mif$sample_id]][1], .before = 1)
    }) %>%
      do.call(dplyr::bind_rows, .) %>%
      dplyr::mutate(Simulations = num_permutations)
    return(list(Dixon_Z=Dixon_Z, 
                Dixon_C=Dixon_C))
  }, mc.cores = workers, mc.allow.recursive = T, mc.preschedule = F)
  
  #The append path used to do `mutate(Run = max(Run) + 1)` where `Run` was expected
  #to be a column of the freshly computed table -- but `Run` was only ever created
  #in the overwrite branch, so `max(Run)` resolved lexically and threw
  #"object 'Run' not found". overwrite = FALSE is the default, so dixons_s() was
  #broken on its default arguments. write_derived() is the same helper the metric
  #functions use.
  if("Z" %in% type){
    mif = write_derived(mif, "Dixon_Z",
                        do.call(dplyr::bind_rows, lapply(out, `[[`, "Dixon_Z")),
                        overwrite)
  }
  if("C" %in% type){
    mif = write_derived(mif, "Dixon_C",
                        do.call(dplyr::bind_rows, lapply(out, `[[`, "Dixon_C")),
                        overwrite)
  }

  structure(mif, class="mif")
}

