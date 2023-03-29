
#' Nearest Neighbor Based Measures of Spatial Clustering for IF data
#'
#' @description For a given cell type, this function computes proportion of cells 
#' that have nearest neighbor less than r for the observed and permuted point processes.
#' @param mif An MIF object
#' @param mnames Character vector of marker names to estimate degree of 
#' nearest neighbor distribution
#' @param r_range Numeric vector of potential r values this range must include 0.
#' @param edge_correction Character value indicating the type of edge correction 
#'  to use. Options include "rs" or "hans". 
#' @param num_permutations Numeric value indicating the number of permutations used. 
#'  Default is 50.   
#' @param keep_perm_dis Logical value determining whether or not to keep the full 
#'  distribution of permuted G values
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
#'    \item{Permuted CSR}{Average observed G for the permuted point 
#'    process}
#'    \item{Observed}{Observed valuefor the observed point process}
#'    \item{Degree of Clustering Permuted}{Degree of spatial clustering where the
#'    reference is the permuted estimate of CSR}
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
#' # Nearest Neighbor distribution for all markers with a neighborhood size 
#' # of  10,20,...,100 (zero must be included in the input).
#'
#' 
#' x <- NN_G(mif = x, mnames = markers[1:2], num_permutations = 1,
#' edge_correction = 'rs', r = seq(0,100,10),
#' keep_perm_dis = FALSE, workers = 1)
#' 
#' @export


NN_G = function(mif, mnames, r_range = seq(0, 100, 50),
                num_permutations = 50, edge_correction = "rs",
                keep_perm_dis = FALSE, workers = 1,
                overwrite = FALSE, xloc = NULL, yloc = NULL){
  future::plan(future::multisession, workers = workers)
  data = mif$spatial
  id = mif$sample_id
  if(overwrite == FALSE){
    mif$derived$univariate_NN = rbind(mif$derived$univariate_NN ,
                                      furrr::future_map(.x = 1:length(data),
                                                 ~{
                                                   uni_NN_G(data = data[[.x]], num_iters = num_permutations, 
                                                            markers = mnames,  id  = id, 
                                                            correction = edge_correction, r = r_range,
                                                            perm_dist = keep_perm_dis, xloc, yloc)}, 
                                                 .options = furrr::furrr_options(seed=TRUE), .progress = T) %>%
                                        plyr::ldply()
    )
  }else{
    mif$derived$univariate_NN = furrr::future_map(.x = 1:length(data),
                                           ~{
                                             uni_NN_G(data = data[[.x]], num_iters = num_permutations, 
                                                      markers = mnames,  id  = id, 
                                                      correction = edge_correction, r = r_range,
                                                      perm_dist = keep_perm_dis, xloc, yloc)}, 
                                           .options = furrr::furrr_options(seed=TRUE), .progress = T) %>%
      plyr::ldply()
  }
  return(mif)
}


#' Bivariate Nearest Neighbor Based Measures of Spatial Clustering for IF data
#'
#' @description This function computes the nearest neighbor distribution for a 
#' particular marker relative to another marker for the observed and permuted point
#' processes.
#' @param mif An MIF object
#' @param mnames Character vector of marker names to estimate degree of 
#' nearest neighbor distribution
#' @param r_range Numeric vector of potential r values this range must include 0.
#' Note that the range selected is very different than count based measures. 
#' See details. 
#' @param edge_correction Character value indicating the type of edge correction 
#'  to use. Options include "rs" or "hans". 
#' @param num_permutations Numeric value indicating the number of permutations used. 
#'  Default is 50.   
#' @param keep_perm_dis Logical value determining whether or not to keep the full 
#'  distribution of permuted G values
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
#'    \item{Permuted CSR}{Average observed G for the permuted point 
#'    process}
#'    \item{Observed}{Observed value for the observed point process}
#'    \item{Degree of Clustering Permuted}{Degree of spatial clustering where the
#'    reference is the permuted estimate of CSR}
#'    \item{Degree of Clustering Theoretical}{Degree of spatial clustering where the
#'    reference is the theoretical estimate of CSR}
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
#' #Nearest Neighbor distribution for the colocalization of CD3+CD8+ positive 
#' #cells and CD3+FOXP3+ positive cells where CD3+FOXP3+ is the reference cell 
#' #type at neighborhood size of 10,20,...,100 (zero must be included in the 
#' #input).
#' 
#' x <- bi_NN_G(mif = x, mnames = c("CD3..CD8.", "CD3..FOXP3."), 
#' num_permutations = 1, edge_correction = 'rs', r = seq(0,100,10),
#' keep_perm_dis = FALSE, workers = 1, exhaustive = TRUE) 
#' 
#' @export


bi_NN_G = function(mif, mnames, r_range = seq(0, 100, 50),
                   num_permutations = 50, edge_correction = "rs",
                   keep_perm_dis = FALSE, exhaustive = TRUE,
                   workers = 1, overwrite = FALSE, xloc = NULL, yloc = NULL){

  future::plan(future::multisession, workers = workers)
  data = mif$spatial
  id = mif$sample_id
  if(overwrite == FALSE){
    mif$derived$bivariate_NN = rbind(mif$derived$bivariate_NN,
                                     furrr::future_map(.x = 1:length(data),
                                                ~{
                                                  bi_NN_G_sample(data = data[[.x]], num_iters = num_permutations, 
                                                                 markers = mnames, r = r_range, id  = id, 
                                                                 correction = edge_correction,
                                                                 perm_dist = keep_perm_dis,
                                                                 exhaustive, xloc, yloc) %>%
                                                    data.frame(check.names = FALSE)}, 
                                                .options = furrr::furrr_options(seed=TRUE), 
                                                .progress = T
                                     ) %>%
                                       plyr::ldply()
    )
  }else{
    mif$derived$bivariate_NN = furrr::future_map(.x = 1:length(data),
                                          ~{
                                            bi_NN_G_sample(data = data[[.x]], num_iters = num_permutations, 
                                                           markers = mnames, r = r_range, id  = id, 
                                                           correction = edge_correction,
                                                           perm_dist = keep_perm_dis,
                                                           exhaustive, xloc, yloc) %>%
                                              data.frame(check.names = FALSE)}, 
                                          .options = furrr::furrr_options(seed=TRUE), 
                                          .progress = T
    ) %>% plyr::ldply()
  }
  
  return(mif)
}


