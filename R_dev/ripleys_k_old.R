

#' Calculate Count Based Measures of Spatial Clustering for IF data
#'
#' @description This function calculates count based Measures (Ripley's K, Besag 
#'   L, and Marcon's M) of IF data to characterize correlation of spatial point
#'   process.
#' @param mif An MIF object
#' @param mnames Character vector of marker names to estimate degree of 
#' spatial clustering.
#' @param r_range Numeric vector of potential r values this range must include 0. 
#' @param edge_correction A haracter bector indicating the type of edge correction 
#'  to use. Options include "translation" or "isotropic". 
#' @param method Character value indicating which measure (K, L, M) used to 
#' estimate the degree of spatial clustering. Description of the methods can be 
#' found in Details section.
#' @param num_permutations Numeric value indicating the number of permutations used. 
#'  Default is 50.   
#' @param keep_perm_dis Logical value determining whether or not to keep the full 
#'  distribution of permuted K values
#' @param workers Integer value for the number of workers to spawn
#' @param overwrite Logical value determining if you want the results to replace the 
#' current output (TRUE) or be to be appended (FALSE).
#' @param xloc a string corresponding to the x coordinates. If null the average of 
#' XMin and XMax will be used 
#' @param yloc a string corresponding to the y coordinates. If null the average of 
#' YMin and YMax will be used 
#' @importFrom magrittr %>%
#'  
#' @return Returns a data.frame
#'    \item{Theoretical CSR}{Expected value assuming complete spatial randomnessn}
#'    \item{Permuted CSR}{Average observed K, L, or M for the permuted point 
#'    process}
#'    \item{Observed}{Observed valuefor the observed point process}
#'    \item{Degree of Clustering Permuted}{Degree of spatial clustering where the
#'    reference is the permutated estimate of CSR}
#'    \item{Degree of Clustering Theoretical}{Degree of spatial clustering where the
#'    reference is the theoretical estimate of CSR}
#' @examples 
#' #Create mif object
#' library(dplyr)
#' x <- create_mif(clinical_data = example_clinical %>% 
#' mutate(deidentified_id = as.character(deidentified_id)),
#' sample_data = example_summary %>% 
#' mutate(deidentified_id = as.character(deidentified_id)),
#' spatial_list = example_spatial,
#' patient_id = "deidentified_id", 
#' sample_id = "deidentified_sample")
#' 
#' # Define the set of markers to study
#' markers <- c("CD3..Opal.570..Positive","CD8..Opal.520..Positive",
#' "FOXP3..Opal.620..Positive","CD3..CD8.","CD3..FOXP3.")
#' 
#' # Ripley's K for all markers with a neighborhood size 
#' # of  10,20,...,100 (zero must be included in the input).
#' 
#' x <- ripleys_k(mif = x, mnames = markers, num_permutations = 1,
#' edge_correction = 'translation', r = seq(0,100,10),
#' keep_perm_dis = FALSE, workers = 1)
#' 

ripleys_k_old = function(mif, mnames, r_range = seq(0, 100, 50),
                         num_permutations = 50, edge_correction = "translation",
                         method = 'K',keep_perm_dis = FALSE, workers = 1,
                         overwrite = FALSE, xloc = NULL, yloc = NULL){
  
  future::plan(future::multisession, workers = workers)
  data = mif$spatial
  id = mif$sample_id
  if(overwrite == FALSE){
    mif$derived$univariate_Count = rbind(mif$derived$univariate_Count,
                                         furrr::future_map(.x = 1:length(data), ~{
                                           uni_Rip_K(data = data[[.x]], num_iters = num_permutations, r = r_range,
                                                     markers = mnames, id  = id, correction = edge_correction, 
                                                     method = method, perm_dist = keep_perm_dis, xloc, yloc)}, 
                                           .options = furrr::furrr_options(seed=TRUE), .progress = T,
                                           xloc = xloc, yloc = yloc) %>%
                                           plyr::ldply()
    )
  }else{
    mif$derived$univariate_Count = furrr::future_map(.x = 1:length(data), ~{
      uni_Rip_K(data = data[[.x]], num_iters = num_permutations, r = r_range,
                markers = mnames, id  = id, correction = edge_correction, 
                method = method, perm_dist = keep_perm_dis, xloc, yloc)}, 
      .options = furrr::furrr_options(seed=TRUE), .progress = T) %>%
      plyr::ldply()
  }
  return(mif)
}


#' Bivariate Count Based Measures of Spatial Clustering function for IF data
#'
#' @description This function calculates count based Measures (Ripley's K, Besag 
#'   L, and Marcon's M) of IF data to characterize correlation of the observed and permyted spatial point
#'   processes for two markers.
#' @param mif An MIF object
#' @param mnames Character vector of marker names to estimate degree of 
#' spatial clustering. Spatial clustering will be computed between each 
#' combination of markers in this list.
#' @param r_range Numeric vector of potential r values this range must include 0
#' @param edge_correction Character value indicating the type of edge correction 
#'  to use. Options include "theoretical", "translation", "isotropic" or "border". 
#'  Various edges corrections are most appropriate in different settings. Default
#'  is "translation". 
#' @param method Character value indicating which measure (K, L, M) used to 
#' estimate the degree of spatial clustering. Description of the methods can be 
#' found in Details section.
#' @param num_permutations Numeric value indicating the number of permutations used. 
#'  Default is 50.   
#' @param keep_perm_dis Logical value determining whether or not to keep the full 
#'  distribution of permuted K values
#' @param exhaustive Logical. If TRUE then markers must be a vector and spatial 
#' measures will be computed all pairs of unique markers. If FALSE then markers must
#' be a data.frame with the desired combinations.
#' @param workers Integer value for the number of workers to spawn
#' @param overwrite Logical value determining if you want the results to replace the 
#' current output (TRUE) or be to be appended (FALSE).
#' @param xloc a string corresponding to the x coordinates. If null the average of 
#' XMin and XMax will be used 
#' @param yloc a string corresponding to the y coordinates. If null the average of 
#' YMin and YMax will be used 
#' @importFrom magrittr %>%
#' 
#' @return Returns a data frame 
#'    \item{anchor}{Marker for which the distances are measured from}
#'    \item{counted}{Marker for which the distances are measured to}
#'    \item{Theoretical CSR}{Expected value assuming complete spatial randomness}
#'    \item{Permuted CSR}{Average observed K, L, or M for the permuted point 
#'    process}
#'    \item{Observed}{Observed value for the observed point process}
#'    \item{Degree of Clustering Permuted}{Degree of spatial clustering where the
#'    reference is the permuted estimate of CSR}
#'    \item{Degree of Clustering Theoretical}{Degree of spatial clustering where the
#'    reference is the theoretical estimate of CSR}
#'    
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
#' #Ripley's K for the colocalization of CD3+CD8+ positive cells and
#' #CD3+FOXP3+ positive cells where CD3+FOXP3+ is the reference cell type at 
#' #neighborhood size of 10,20,...,100 (zero must be included in the input).
#' 
#' x <- bi_ripleys_k(mif = x, mnames = c("CD3..CD8.", "CD3..FOXP3."), 
#' num_permutations = 1, edge_correction = 'translation', r = seq(0,100,10),
#' keep_perm_dis = FALSE, workers = 1, exhaustive = TRUE) 
#' 


bi_ripleys_k_old <- function(mif,
                             mnames, 
                             r_range = seq(0, 100, 50),
                             num_permutations = 50,
                             edge_correction = "translation",
                             method = 'K',
                             keep_perm_dis = FALSE,
                             exhaustive = TRUE,
                             workers = 1,
                             overwrite = FALSE,
                             xloc = NULL, yloc = NULL){
  
  future::plan(future::multisession, workers = workers)
  data = mif$spatial
  id = mif$sample_id
  if(overwrite == FALSE){
    mif$derived$bivariate_Count = rbind(mif$derived$bivariate_Count,
                                        furrr::future_map(.x = 1:length(data),
                                                          ~{
                                                            bi_Rip_K(data = data[[.x]], num_iters = num_permutations, 
                                                                     markers = mnames, id  = id, r = r_range,
                                                                     correction = edge_correction, method = method, 
                                                                     perm_dist = keep_perm_dis,
                                                                     exhaustive = exhaustive, xloc, yloc) %>%
                                                              data.frame(check.names = FALSE)
                                                          }, .options = furrr::furrr_options(seed=TRUE), .progress = T) %>%
                                          plyr::ldply()
                                        
    )
  }else{
    mif$derived$bivariate_Count = furrr::future_map(.x = 1:length(data),
                                                    ~{
                                                      bi_Rip_K(data = data[[.x]], num_iters = num_permutations, 
                                                               markers = mnames, id  = id, r = r_range,
                                                               correction = edge_correction, method = method, 
                                                               perm_dist = keep_perm_dis,
                                                               exhaustive = exhaustive, xloc, yloc) %>%
                                                        data.frame(check.names = FALSE)
                                                    }, .options = furrr::furrr_options(seed=TRUE), .progress = T) %>%
      plyr::ldply()
  }
  return(mif)
}


#' #' Dixon's S Segregation Statistic
#' #'
#' #' @description This function processes the spatial files in the mif object,
#' #' requiring a column that distinguishes between different groups i.e. tumor and 
#' #' stroma
#' #' @param mif An MIF object
#' #' @param mnames vector of markers corresponding to spatial columns to check Dixon's S between
#' #' @param num_permutations Numeric value indicating the number of permutations used. 
#' #'  Default is 1000.
#' #' @param type a character string for the type that is wanted in the output which can
#' #' be "Z" for z-statistic results or "C" for Chi-squared statistic results
#' #' @param workers Integer value for the number of workers to spawn
#' #' @param overwrite Logical value determining if you want the results to replace the 
#' #' current output (TRUE) or be to be appended (FALSE).
#' #' @param xloc a string corresponding to the x coordinates. If null the average of 
#' #' XMin and XMax will be used 
#' #' @param yloc a string corresponding to the y coordinates. If null the average of 
#' #' YMin and YMax will be used 
#' #' @importFrom magrittr %>%
#' #' 
#' #' @return Returns a data frame for Z-statistic
#' #'    \item{From}{}
#' #'    \item{To}{}
#' #'    \item{Obs.Count}{}
#' #'    \item{Exp. Count}{}
#' #'    \item{S}{}
#' #'    \item{Z}{}
#' #'    \item{p-val.Z}{}
#' #'    \item{p-val.Nobs}{}
#' #'    \item{Marker}{}
#' #'    \item{Classifier Labeled Column Counts}{}
#' #'    \item{Image.Tag}{}
#' #' @return Returns a data frame for C-statistic
#' #'    \item{Segregation}{}
#' #'    \item{df}{}
#' #'    \item{Chi-sq}{}
#' #'    \item{P.asymp}{}
#' #'    \item{P.rand}{}
#' #'    \item{Marker}{}
#' #'    \item{Classifier Labeled Column Counts}{}
#' #'    \item{Image.Tag}{}
#' #' @examples 
#' #' #' #Create mif object
#' #' library(dplyr)
#' #' x <- create_mif(clinical_data = example_clinical %>% 
#' #' mutate(deidentified_id = as.character(deidentified_id)),
#' #' sample_data = example_summary %>% 
#' #' mutate(deidentified_id = as.character(deidentified_id)),
#' #' spatial_list = example_spatial,
#' #' patient_id = "deidentified_id", 
#' #' sample_id = "deidentified_sample")
#' #' 
#' #' @export
#' 
#' dixons_s = function(mif, mnames, num_permutations = 1000, type = c("Z", "C"),
#'                     workers = 1, overwrite = FALSE, xloc = NULL, yloc = NULL){
#'   #dix_s_z
#'   #dix_s_c
#'   #make sure to assign the spatial data name to the output of dix_s_z and dix_s_c
#'   #check classifier label
#'   
#'   data = mif$spatial
#'   mnames = mnames %>%
#'     expand.grid(., .) %>%
#'     dplyr::filter(Var1 != Var2) %>% 
#'     dplyr::rowwise() %>%
#'     dplyr::mutate(Var3 = paste0(sort(c(Var1, Var2)), collapse = ",")) %>%
#'     distinct(Var3, .keep_all = TRUE) %>%
#'     select(1, 2)
#'   
#'   future::plan(future::multisession, workers = workers)
#'   if(overwrite == FALSE){
#'     if("Z" %in% type){
#'       mif$derived$dixons_S_Z = rbind(mif$derived$dixons_S_Z,
#'                                      furrr::future_map(.x = 1:length(data),
#'                                                        ~{
#'                                                          dix_s_z(data = data[[.x]], num_permutations = num_permutations,
#'                                                                  markers = mnames, xloc = xloc, yloc = yloc) %>%
#'                                                            dplyr::mutate(Image.Tag = names(data)[.x])
#'                                                        }, .options = furrr::furrr_options(seed = TRUE),
#'                                                        .progress = TRUE) %>%
#'                                        plyr::ldply())
#'     }
#'     if("C" %in% type){
#'       mif$derived$dixons_S_C = rbind(mif$derived$dixons_S_C,
#'                                      furrr::future_map(.x = 1:length(data),
#'                                                        ~{
#'                                                          dix_s_c(data = data[[.x]], num_permutations = num_permutations,
#'                                                                  markers = mnames, classifier_label = classifier_label, 
#'                                                                  xloc = xloc, yloc = yloc) %>%
#'                                                            mutate(Image.Tag = names(data)[.x])
#'                                                        }, .options = furrr::furrr_options(seed =TRUE),
#'                                                        .progress =TRUE) %>%
#'                                        plyr::ldply())
#'     }
#'   } else {
#'     if("Z" %in% type){
#'       mif$derived$dixons_S_Z = furrr::future_map(.x = 1:length(data),
#'                                                  ~{
#'                                                    dix_s_z(data = data[[.x]], num_permutations = num_permutations,
#'                                                            markers = mnames, classifier_label = classifier_label, 
#'                                                            xloc = xloc, yloc = yloc) %>%
#'                                                      mutate(Image.Tag = names(data)[.x])
#'                                                  }, .options = furrr::furrr_options(seed =TRUE),
#'                                                  .progress =TRUE) %>%
#'         plyr::ldply()
#'     }
#'     if("C" %in% type){
#'       mif$derived$dixons_S_C = furrr::future_map(.x = 1:length(data),
#'                                                  ~{
#'                                                    dix_s_c(data = data[[.x]], num_permutations = num_permutations,
#'                                                            markers = mnames, classifier_label = classifier_label, 
#'                                                            xloc = xloc, yloc = yloc) %>%
#'                                                      mutate(Image.Tag = names(data)[.x])
#'                                                  }, .options = furrr::furrr_options(seed =TRUE),
#'                                                  .progress =TRUE) %>%
#'         plyr::ldply()
#'     }
#'   }
#'   structure(mif, class="mif")
#' }