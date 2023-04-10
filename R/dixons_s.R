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
#' @return Returns a data frame for Z-statistic
#'    \item{From}{}
#'    \item{To}{}
#'    \item{Obs.Count}{}
#'    \item{Exp. Count}{}
#'    \item{S}{}
#'    \item{Z}{}
#'    \item{p-val.Z}{}
#'    \item{p-val.Nobs}{}
#'    \item{Marker}{}
#'    \item{Classifier Labeled Column Counts}{}
#'    \item{Image.Tag}{}
#' @return Returns a data frame for C-statistic
#'    \item{Segregation}{}
#'    \item{df}{}
#'    \item{Chi-sq}{}
#'    \item{P.asymp}{}
#'    \item{P.rand}{}
#'    \item{Marker}{}
#'    \item{Classifier Labeled Column Counts}{}
#'    \item{Image.Tag}{}
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
  #dix_s_z
  #dix_s_c
  #make sure to assign the spatial data name to the output of dix_s_z and dix_s_c
  #check classifier label
  
  data = mif$spatial
  mnames = mnames %>%
    expand.grid(., .) %>%
    dplyr::filter(Var1 != Var2) %>% 
    dplyr::rowwise() %>%
    dplyr::mutate(Var3 = paste0(sort(c(Var1, Var2)), collapse = ",")) %>%
    distinct(Var3, .keep_all = TRUE) %>%
    select(1, 2)
  
  future::plan(future::multisession, workers = workers)
  if(overwrite == FALSE){
    if("Z" %in% type){
      mif$derived$dixons_S_Z = rbind(mif$derived$dixons_S_Z,
                                     furrr::future_map(.x = 1:length(data),
                                                       ~{
                                                         dix_s_z(data = data[[.x]], num_permutations = num_permutations,
                                                                 markers = mnames, xloc = xloc, yloc = yloc) %>%
                                                           dplyr::mutate(Image.Tag = names(data)[.x])
                                                       }, .options = furrr::furrr_options(seed = TRUE),
                                                       .progress = TRUE) %>%
                                       plyr::ldply())
    }
    if("C" %in% type){
      mif$derived$dixons_S_C = rbind(mif$derived$dixons_S_C,
                                     furrr::future_map(.x = 1:length(data),
                                                       ~{
                                                         dix_s_c(data = data[[.x]], num_permutations = num_permutations,
                                                                 markers = mnames, classifier_label = classifier_label, 
                                                                 xloc = xloc, yloc = yloc) %>%
                                                           mutate(Image.Tag = names(data)[.x])
                                                       }, .options = furrr::furrr_options(seed =TRUE),
                                                       .progress =TRUE) %>%
                                       plyr::ldply())
    }
  } else {
    if("Z" %in% type){
      mif$derived$dixons_S_Z = furrr::future_map(.x = 1:length(data),
                                                 ~{
                                                   dix_s_z(data = data[[.x]], num_permutations = num_permutations,
                                                           markers = mnames, classifier_label = classifier_label, 
                                                           xloc = xloc, yloc = yloc) %>%
                                                     mutate(Image.Tag = names(data)[.x])
                                                 }, .options = furrr::furrr_options(seed =TRUE),
                                                 .progress =TRUE) %>%
        plyr::ldply()
    }
    if("C" %in% type){
      mif$derived$dixons_S_C = furrr::future_map(.x = 1:length(data),
                                                 ~{
                                                   dix_s_c(data = data[[.x]], num_permutations = num_permutations,
                                                           markers = mnames, classifier_label = classifier_label, 
                                                           xloc = xloc, yloc = yloc) %>%
                                                     mutate(Image.Tag = names(data)[.x])
                                                 }, .options = furrr::furrr_options(seed =TRUE),
                                                 .progress =TRUE) %>%
        plyr::ldply()
    }
  }
  structure(mif, class="mif")
}