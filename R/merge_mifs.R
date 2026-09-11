#' Merge several MIF objects together
#'
#' @description This function merges MIF objects that were run separately so they
#' can be used as a single MIF. MIF objects don't *need* but *should* have the same
#' column names in the summary file and clinical data file. The MIF objects **DO**
#' need to have the same patient_id and sample_id.
#' @param mifs A list of MIF objects to merge together
#' @param check.names whether to check names of spatial files and summary enttries
#'  
#' @return Returns a new MIF object list
#'    \item{clinical_data}{clinical information from all}
#'    \item{sample}{cell level summary data from all}
#'    \item{spatial}{contains all spatial files from all MIFs}
#'    \item{derived}{appended derived variables}
#'    \item{patient_id}{patient_id from the first MIF - this is
#'    why it is important to have the same patient_id for all MIFs}
#'    \item{sample_id}{sample_id from the first MIF - also important
#'    for all MIFs to have the same sample_id}
#' @examples 
#' #merge several MIF objects
#' library(dplyr)
#' x <- create_mif(clinical_data = example_clinical %>% 
#' mutate(deidentified_id = as.character(deidentified_id)),
#' sample_data = example_summary %>% 
#' mutate(deidentified_id = as.character(deidentified_id)),
#' spatial_list = example_spatial,
#' patient_id = "deidentified_id", 
#' sample_id = "deidentified_sample")
#' x <- merge_mifs(mifs = list(x, x), check.names = FALSE)
#' 
#'@export


merge_mifs = function(mifs = NULL, check.names = T){
  #check for proper number of MIF objects
  #`length(mifs) == 1` missed the empty case, which then reached
  #`names(sizes) = seq(length(sizes))` -> seq(0) is c(1, 0), length 2, giving
  #"'names' attribute [2] must be the same length as the vector [0]".
  if(is.null(mifs) || length(mifs) < 2){
    stop("Please enter at least 2 MIF objects to merge")
  }
  if(!all(vapply(mifs, inherits, logical(1), what = "mif"))){
    stop("Every element of `mifs` must be a mif object created with `create_mif()`")
  }
  #find patient_id values
  patient_id = sapply(mifs, function(mif){
    mif$patient_id
  }) %>% unique()
  if(length(patient_id) > 1 & check.names == T){
    stop("MIF objects have different patient_id values")
  }
  #find sample_id values
  sample_id = sapply(mifs, function(mif){
    mif$sample_id
  }) %>% unique()
  if(length(sample_id) > 1 & check.names == T){
    stop("MIF objects have different sample_id values")
  }
  
  #merge clinical variables
  clinical = lapply(mifs, function(mif){
    mif$clinical
  }) %>%
    do.call(dplyr::bind_rows, .) %>%
    dplyr::distinct()
  #merge clinical variables
  sample = lapply(mifs, function(mif){
    mif$sample
  }) %>%
    do.call(dplyr::bind_rows, .) %>%
    dplyr::distinct()
  
  #merging spatial data
  spatial = lapply(mifs, function(mif){
    mif$spatial
  }) %>% 
    do.call(c, .)
  if(TRUE %in% duplicated(names(spatial)) & check.names == TRUE){
    stop("Multiple files have the same name")
  }
  
  #merging the derived variables
  derived_names = sapply(mifs, function(mif){
    names(mif$derived)
  }) %>%
    unlist() %>%
    as.character() %>%
    unique()
  #alert no derived
  #`is.null()` never fired: as.character(unlist(list(NULL, NULL))) is character(0),
  #not NULL, so this message was unreachable for mifs with empty derived slots --
  #exactly the case it was written for.
  if(!length(derived_names)){
    message("No variables have been derived yet")
  }
  #find mif with the most to use as "base" mif for merging
  sizes = sapply(mifs, function(mif){
    length(mif$derived)
  })
  names(sizes) = seq_along(sizes)
  derived = lapply(mifs, function(mif){
    mif$derived
  })
  names(derived) = names(sizes)
  #sort derived to have order with most derived first
  derived = derived[names(sort(sizes, decreasing = T))]
  #begin merging
  #`bind_rows()` on every derived slot "succeeds" even on list-valued slots
  #(spatial_plots, and split_tissue()'s density_boundary) by silently
  #collapsing a list keyed by sample name into a nameless table -- the worst
  #possible outcome. Type-check first: bind_rows() only when every part is a
  #data frame, otherwise concatenate the lists and carry provenance forward.
  derived2 = lapply(derived_names, function(name){
    parts = lapply(derived, function(mif){
      mif[[name]]
    })
    parts = parts[!vapply(parts, is.null, logical(1))]
    if(!length(parts)){
      return(NULL)
    }
    if(all(vapply(parts, is.data.frame, logical(1)))){
      return(do.call(dplyr::bind_rows, parts))
    }
    #`parts` is named by mif index ("1", "2", ...); do.call(c, parts) treats
    #those as argument names, so c() prefixes each inner sample name with the
    #mif index ("1.S1") instead of concatenating the lists' own names.
    merged = do.call(c, unname(parts))
    if(check.names && anyDuplicated(names(merged))){
      stop("Multiple mifs' `derived$", name, "` share a sample name; merging would ",
           "silently overwrite one with another.")
    }
    call_infos = lapply(parts, attr, "call_info")
    call_infos = call_infos[!vapply(call_infos, is.null, logical(1))]
    if(length(call_infos) > 1 &&
      !all(vapply(call_infos[-1], identical, logical(1), call_infos[[1]]))){
      warning("Merged mifs' `derived$", name, "` were computed with different settings ",
             "(differing `call_info`); keeping the first mif's settings on the merged result.")
    }
    if(length(call_infos)){
      attr(merged, "call_info") = call_infos[[1]]
    }
    merged
  })
  rm(derived)
  names(derived2) = derived_names
  mif <- list(clinical = clinical,
              sample = sample,
              spatial = spatial,
              derived = derived2,
              patient_id = patient_id,
              sample_id = sample_id)
  
  structure(mif, class="mif")
}
