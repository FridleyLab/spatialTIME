#' Calculate Area under Ripley's K curve
#'
#' @param mif object of class `mif` made with `createMIF`
#'
#' @return object of class `mif` with `univariate_AUC` in derived slot
#' @export
#'
#' @examples
#' x <- spatialTIME::create_mif(clinical_data =spatialTIME::example_clinical %>% 
#'   dplyr::mutate(deidentified_id = as.character(deidentified_id)),
#'   sample_data = spatialTIME::example_summary %>% 
#'   dplyr::mutate(deidentified_id = as.character(deidentified_id)),
#'   spatial_list = spatialTIME::example_spatial,
#'   patient_id = "deidentified_id", 
#'   sample_id = "deidentified_sample")
#' mnames = "CD3..Opal.570..Positive"
#' x2 = ripleys_k2(mif = x,
#'   mnames = mnames, 
#'   r_range = seq(0, 100, 1), 
#'   num_permutations = 100,
#'   edge_correction = "translation", 
#'   method = "K", 
#'   keep_perm_dis =FALSE, 
#'   workers = 1, 
#'   overwrite =TRUE)
#' x2 = ripleys_k_auc(x2)
#' 
ripleys_k_auc = function(mif){
  if(!inherits(mif, "mif")){
    stop("Input mIF must be of class `mif`.\n\tCreate mIF object with function `createMIF()`")
  }
  
  if(!exists("univariate_Count", where = mif$derived)){
    stop("Must run `ripleys_k()` before determining area.")
  }
  #summarise univariate K in the case of permutations used
  sum_dat = mif$derived$univariate_Count %>% 
    dplyr::mutate(iter = dplyr::case_when(iter == "Estimate" ~ iter,
                                          T ~ "Permutation")) %>%
    dplyr::group_by(Run, iter, Label, Marker, r) %>% 
    dplyr::summarise_all(~mean(., na.rm=T)) %>% 
    dplyr::group_by(Run, iter, Label, Marker)
  #identify different groups
  groups = expand.grid(Run = unique(sum_dat$Run),
                       iter = unique(sum_dat$iter),
                       Label = unique(sum_dat$Label),
                       Marker = unique(sum_dat$Marker))
  
  #subset data to group and then compute area for columns
  out = lapply(data.frame(t(groups)), function(group){
    df = sum_dat %>% 
      dplyr::filter(Run == group[1],
                    iter == group[2],
                    Label == group[3],
                    Marker == group[4])
    df %>%
      dplyr::summarise(across(`Theoretical K`:`Degree of Clustering Theoretical`, 
                              list(AUC = ~ ifelse(all(is.na(.x)),
                                                   NA_real_,
                                                   flux::auc(x = r[!is.na(.x)],
                                                             y = .x[!is.na(.x)]))),
                              .names = "{.col} AUC"),
                       .groups = "keep")
    
  }) %>%
    do.call(dplyr::bind_rows, .)
  
  #add area output to derived_variables slot in mif
  mif$derived$univariate_AUC = out
  
  return(mif)
}
