#' Calculate Count Based Measures of Spatial Clustering for IF data
#'
#' @description This function calculates count based Measures (Ripley's K, Besag 
#'   L, and Marcon's M) of IF data to characterize correlation of spatial point
#'   process.
#' @param mif An MIF object
#' @param mnames Character vector of marker names to estimate degree of 
#' spatial clustering.
#' @param r_range Numeric vector of potential r values this range must include 0. 
#' @param edge_correction Character value indicating the type of edge correction 
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
#' @export
#'

ripleys_k = function(mif, mnames, r_range = seq(0, 100, 50),
                        num_permutations = 50, edge_correction = "translation",
                        method = 'K',keep_perm_dis = FALSE, workers = 1,
                        overwrite = FALSE){
  require(furrr)
  require(magrittr)
  plan(multisession, workers = workers)
  data = mif$spatial
  id = mif$sample_id
  if(overwrite == FALSE){
    mif$derived$univariate_Count = rbind(mif$derived$univariate_Count,
                                         future_map(.x = 1:length(data), ~{
                                           uni_Rip_K(data = data[[.x]], num_iters = num_permutations, r = r_range,
                                                     markers = mnames, id  = id, correction = edge_correction, 
                                                     method = method, perm_dist = keep_perm_dis)}, 
                                           .options = furrr_options(seed=TRUE), .progress = T) %>%
                                           plyr::ldply()
    )
  }else{
    mif$derived$univariate_Count = future_map(.x = 1:length(data), ~{
      uni_Rip_K(data = data[[.x]], num_iters = num_permutations, r = r_range,
                markers = mnames, id  = id, correction = edge_correction, 
                method = method, perm_dist = keep_perm_dis)}, 
      .options = furrr_options(seed=TRUE), .progress = T) %>%
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
#' @export
#' 
bi_ripleys_k <- function(mif,
                            mnames, 
                            r_range = seq(0, 100, 50),
                            num_permutations = 50,
                            edge_correction = "translation",
                            method = 'K',
                            keep_perm_dis = FALSE,
                            exhaustive = TRUE,
                            workers = 1,
                            overwrite = FALSE){
  require(furrr)
  require(magrittr)
  plan(multisession, workers = workers)
  data = mif$spatial
  id = mif$sample_id
  if(overwrite == FALSE){
    mif$derived$bivariate_Count = rbind( mif$derived$bivariate_Count,
                                         furrr::future_map(.x = 1:length(data),
                                                           ~{
                                                             bi_Rip_K(data = data[[.x]], num_iters = num_permutations, 
                                                                      markers = mnames, id  = id, r = r_range,
                                                                      correction = edge_correction, method = method, 
                                                                      perm_dist = keep_perm_dis,
                                                                      exhaustive = exhaustive) %>%
                                                               data.frame(check.names = FALSE)
                                                           }, .options = furrr_options(seed=TRUE), .progress = T) %>%
                                           plyr::ldply()
                                         
    )
  }else{
    mif$derived$bivariate_Count = furrr::future_map(.x = 1:length(data),
                                                    ~{
                                                      bi_Rip_K(data = data[[.x]], num_iters = num_permutations, 
                                                               markers = mnames, id  = id, r = r_range,
                                                               correction = edge_correction, method = method, 
                                                               perm_dist = keep_perm_dis,
                                                               exhaustive = exhaustive) %>%
                                                        data.frame(check.names = FALSE)
                                                    }, .options = furrr_options(seed=TRUE), .progress = T) %>%
      plyr::ldply()
  }
  return(mif)
}

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
#' @return Returns a data.frame
#'    \item{Theoretical CSR}{Expected value assuming complete spatial randomnessn}
#'    \item{Permuted CSR}{Average observed G for the permuted point 
#'    process}
#'    \item{Observed}{Observed valuefor the observed point process}
#'    \item{Degree of Clustering Permuted}{Degree of spatial clustering where the
#'    reference is the permuted estimate of CSR}
#'    \item{Degree of Clustering Theoretical}{Degree of spatial clustering where the
#'    reference is the theoretical estimate of CSR}
#' @export
#'
NN_G = function(mif, mnames, r_range = seq(0, 100, 50),
                num_permutations = 50, edge_correction = "rs",
                keep_perm_dis = FALSE, workers = 1,
                overwrite = FALSE){
  require(furrr)
  require(magrittr)
  plan(multisession, workers = workers)
  data = mif$spatial
  id = mif$sample_id
  if(overwrite == FALSE){
    mif$derived$univariate_NN = rbind(mif$derived$univariate_NN ,
                                      future_map(.x = 1:length(data),
                                                 ~{
                                                   uni_NN_G(data = data[[.x]], num_iters = num_permutations, 
                                                            markers = mnames,  id  = id, 
                                                            correction = edge_correction, r = r_range,
                                                            perm_dist = keep_perm_dis)}, 
                                                 .options = furrr_options(seed=TRUE), .progress = T) %>%
                                        plyr::ldply()
    )
  }else{
    mif$derived$univariate_NN = future_map(.x = 1:length(data),
                                           ~{
                                             uni_NN_G(data = data[[.x]], num_iters = num_permutations, 
                                                      markers = mnames,  id  = id, 
                                                      correction = edge_correction, r = r_range,
                                                      perm_dist = keep_perm_dis)}, 
                                           .options = furrr_options(seed=TRUE), .progress = T) %>%
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
#'@export
#'
bi_NN_G = function(mif, mnames, r_range = seq(0, 100, 50),
                   num_permutations = 50, edge_correction = "rs",
                   keep_perm_dis = FALSE, exhaustive = TRUE,
                   workers = 1, overwrite = FALSE){
  require(furrr)
  require(magrittr)
  plan(multisession, workers = workers)
  data = mif$spatial
  id = mif$sample_id
  if(overwrite == FALSE){
    mif$derived$bivariate_NN = rbind(mif$derived$bivariate_NN,
                                     future_map(.x = 1:length(data),
                                                ~{
                                                  bi_NN_G_sample(data = data[[.x]], num_iters = num_permutations, 
                                                                 markers = mnames, r = r_range, id  = id, 
                                                                 correction = edge_correction,
                                                                 perm_dist = keep_perm_dis,
                                                                 exhaustive) %>%
                                                    data.frame(check.names = FALSE)}, 
                                                .options = furrr_options(seed=TRUE), 
                                                .progress = T
                                     ) %>%
                                       plyr::ldply()
    )
  }else{
    mif$derived$bivariate_NN = future_map(.x = 1:length(data),
                                          ~{
                                            bi_NN_G_sample(data = data[[.x]], num_iters = num_permutations, 
                                                           markers = mnames, r = r_range, id  = id, 
                                                           correction = edge_correction,
                                                           perm_dist = keep_perm_dis,
                                                           exhaustive) %>%
                                              data.frame(check.names = FALSE)}, 
                                          .options = furrr_options(seed=TRUE), 
                                          .progress = T
    ) %>% plyr::ldply()
  }
  
  return(mif)
}
