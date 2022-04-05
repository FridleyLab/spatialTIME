#' Calculate Count Based Measures and NN Measures of Spatial Clustering for IF data
#'
#' @description This function calculates count based Measures (Ripley's K, Besag 
#'   L, and Marcon's M) of IF data to characterize correlation of spatial point
#'   process. For neareast neighbor calculations of a given cell type, this function 
#'   computes proportion of cells that have nearest neighbor less than r for the 
#'   observed and permuted point processes.
#' @param mif An MIF object
#' @param mnames Character vector of marker names to estimate degree of 
#' spatial clustering.
#' @param r_range Numeric vector of potential r values this range must include 0. 
#' @param edge_correction Character vector indicating the type of edge correction 
#'  to use. Options for count based include "translation" or "isotropic" and for 
#'  nearest neighboroOptions include "rs" or "hans". 
#' @param method Character vector indicating which count based measure (K, BiK, 
#' G, BiG) used to estimate the degree of spatial clustering. Description of the 
#' methods can be found in Details section.
#' @param num_permutations Numeric value indicating the number of permutations used. 
#'  Default is 50.   
#' @param keep_perm_dis Logical value determining whether or not to keep the full 
#'  distribution of permuted K or G values
#' @param workers Integer value for the number of workers to spawn
#' @param overwrite Logical value determining if you want the results to replace the 
#' current output (TRUE) or be to be appended (FALSE).
#' @param k_trans Character value of the transformation to apply to count based 
#' metrics (none, M, or L)
#' @param xloc a string corresponding to the x coordinates. If null the average of 
#' XMin and XMax will be used 
#' @param yloc a string corresponding to the y coordinates. If null the average of 
#' YMin and YMax will be used 
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
#' # Ripley's K and nearest neighbor G for all markers with a neighborhood size 
#' # of  10,20,...,100 (zero must be included in the input).
#' 
#' x <- compute_metrics(x, mnames, r_range = seq(0, 100, 10),
#' num_permutations = 10, edge_correction = c("translation", "rs"),
#' method = c("K", "G", "BiK", "BiG"), k_trans = "none", keep_perm_dis = T,
#' workers = 1, overwrite = T, exhaustive = T)
#' 
#'@export

compute_metrics = function(mif, mnames, r_range = seq(0, 100, 50),
                           num_permutations = 50, edge_correction = c("translation"),
                           method = c("K"), k_trans = "none", keep_perm_dis = FALSE, 
                           workers = 1, overwrite = FALSE, xloc = NULL, yloc = NULL,
                           exhaustive = T){
  #create number of workers
  future::plan(future::multisession, workers = workers)
  #set working variables
  data = mif$spatial
  id = mif$sample_id
  correction = c() #making new correction variable for edge correction methods
  if(length(edge_correction) > 2){
    stop("Please specify an edge correction for the methods of interest")
  }
  if(sum(c("translation", "isotropic") %in% edge_correction) == 2){
    stop("Please specify only 1 edge correction method for K")
  }
  if(sum(c("rs", "hans") %in% edge_correction) == 2){
    stop("Please specify only 1 edge correction method for G")
  }
  if("K" %in% method | "BiK" %in% method){
    if("translation" %in% edge_correction | "isotropic" %in% edge_correction){
      correction = append(correction, c("translation", "isotropic")[which(c("translation", "isotropic") %in% edge_correction)])
    } else {
      message("K was specified but no (or incorrect) edge correction - Defaulting to translation")
      correction = append(correction, "translation")
    }
  }
  if("G" %in% method | "BiG" %in% method){
    if("rs" %in% edge_correction | "hans" %in% edge_correction){
      correction = append(correction, c("rs", "hans")[which(c("rs", "hans") %in% edge_correction)])
    } else {
      message("G was specified but no (or incorrect) edge correction - Defaulting to rs")
      correction = append(correction, "rs")
    }
  }
  if(k_trans == "none"){
    k_trans = "K"
  }
  #create table for metrics to calculate
  calc_grid = expand.grid(Method = method, `K Transformation` = k_trans, 
                          `Range` = r_range, Permutation = rep(1, num_permutations), 
                            `Spatial File` = seq(data)) %>%
    dplyr::filter(Range != 0)
  
  calc_grid$`Correction Method` = ifelse(calc_grid$Method %in% c("K", "BiK"), 
                                         c("translation", "isotropic")[which(c("translation", "isotropic") %in% correction)],
                                         c("rs", "hans")[which(c("rs", "hans") %in% correction)]) 
  
  #setting up derived tables
  univariate_Count = data.frame(matrix(ncol = 9, nrow = 0))
  colnames(univariate_Count) = c("iter", id, "Marker", "r", "Theoretical CSR", 
                                 "Permuted K", "Observed K", 
                                 "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  bivariate_Count = data.frame(matrix(ncol = 10, nrow = 0))
  colnames(bivariate_Count) = c("iter", id, "anchor", "counted", "r", "Theoretical CSR", 
                                "Permuted K", "Observed K", 
                                "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  univariate_NN = data.frame(matrix(ncol = 9, nrow = 0))
  colnames(univariate_NN ) = c("iter", id, "Marker", "r", "Theoretical CSR", 
                               "Permuted K", "Observed K", 
                               "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  bivariate_NN = data.frame(matrix(ncol = 10, nrow = 0))
  colnames(bivariate_NN) = c("iter", id, "anchor", "counted", "r", "Theoretical CSR", 
                             "Permuted K", "Observed K", 
                             "Degree of Clustering Permutation", "Degree of Clustering Theoretical")
  
  #paralleled computing of metrics
  derived = furrr::future_map(.x = 1:nrow(calc_grid), ~{
    #set info to unique combination of what to calculate
    info = calc_grid[.x,]
    #assign("info", info, envir = .GlobalEnv)
    if(info$Method == "K"){
      univariate_Count = rbind(univariate_Count,
                               uni_Rip_K(data = data[[info$`Spatial File`]], markers = mnames_good, id = id,
                                         num_iters = info$Permutation, r = c(0, info$Range),
                                         correction = info$`Correction Method`, method = info$`K Transformation`,
                                         perm_dist = keep_perm_dis, xloc, yloc))
    }
    if(info$Method == "BiK"){
      bivariate_Count = rbind(bivariate_Count,
                              bi_Rip_K(data = data[[info$`Spatial File`]], num_iters = info$Permutation,
                                       markers = mnames, id  = id, r = c(0, info$Range),
                                       correction = info$`Correction Method`, method = info$`K Transformation`,
                                       perm_dist = keep_perm_dis, exhaustive = exhaustive, xloc, yloc))
    }
    if(info$Method == "G"){
      univariate_NN = rbind(univariate_NN,
                            uni_NN_G(data = data[[info$`Spatial File`]], num_iters = info$Permutation,
                                     markers = mnames,  id  = id, correction = info$`Correction Method`,
                                     r = c(0, info$Range), perm_dist = keep_perm_dis, xloc, yloc))
    }
    if(info$Method == "BiG"){
      bivariate_NN = rbind(bivariate_NN,
                           bi_NN_G_sample(data = data[[info$`Spatial File`]], num_iters = info$Permutation,
                                          markers = mnames, r = c(0, info$Range), id  = id,
                                          correction = info$`Correction Method`,
                                          perm_dist = keep_perm_dis,
                                          exhaustive = exhaustive, xloc, yloc))
    }
    return(list(UniK = univariate_Count,
                BiK = bivariate_Count,
                UniG = univariate_NN,
                BiK = bivariate_NN))
  }, .options = furrr::furrr_options(seed=TRUE), .progress = T)
  #collapse tables from permutations and methods
  derived2 = lapply(seq(derived[[1]]), function(metric){
    lapply(seq(derived), function(run){
      tmp = derived[[run]][[metric]]
      if(nrow(tmp) == 0 ){
        return()
      } else {
        return(tmp)
      }
    }) %>%
      do.call(bind_rows,.)
  })
  #remove derived and name new tables in derived2
  rm(derived)
  key = c("univariate_Count" = "K", "bivariate_Count" = "BiK",
          "univariate_NN" = "G", "bivariate_NN" = "BiG")
  names(derived2) = names(key)[method %in% key]
  
  #if keep permutation, set permutation number in table
  if(keep_perm_dis){
    derived = lapply(derived2, function(tmp){
      tmp %>%
        dplyr::group_by(dplyr::across(dplyr::any_of(c(id, "Marker", "anchor", "counted", "r")))) %>%
        dplyr::mutate(iter = 1:n(), .before = !!id)
        
    })
  } else { #if not keep permutation, find mean of all runs
    derived = lapply(derived2, function(tmp){
      tmp %>%
        dplyr::group_by(dplyr::across(dplyr::any_of(c(id, "Marker", "anchor", "counted", "r")))) %>%
        dplyr::summarize_all(~mean(., na.rm = TRUE)) %>%
        dplyr::mutate(iter = paste0("Mean of ", num_permutations, " permutations"))
    })
  }
  rm(derived2)
  #Overwrite or Not
  if(overwrite == FALSE){
    if("K" %in% method){
      if(!exists("mif$derived$univariate_Count")){
        mif$derived$univariate_Count = derived$univariate_Count
      } else {
        mif$derived$univariate_Count = dplyr::bind_rows(mif$derived$univariate_Count,
                                                 derived$univariate_Count)
      }
    }
    if("BiK" %in% method){
      if(!exists("mif$derived$univariate_Count")){
        mif$derived$bivariate_Count = derived$bivariate_Count
      } else {
        mif$derived$bivariate_Count = dplyr::bind_rows(mif$derived$bivariate_Count,
                                                derived$bivariate_Count)
      }
    }
    if("G" %in% method){
      if(!exists("mif$derived$univariate_Count")){
        mif$derived$univariate_NN = derived$univariate_NN
      } else {
        mif$derived$univariate_NN = dplyr::bind_rows(mif$derived$univariate_NN,
                                              derived$univariate_NN)
      }
    }
    if("BiG" %in% method){
      if(!exists("mif$derived$univariate_Count")){
        mif$derived$bivariate_NN = derived$bivariate_NN
      } else {
        mif$derived$bivariate_NN = dplyr::bind_rows(mif$derived$bivariate_NN,
                                             derived$bivariate_NN)
      }
    }
  } else {
    if("K" %in% method){
      mif$derived$univariate_Count = derived$univariate_Count
    }
    if("BiK" %in% method){
      mif$derived$bivariate_Count = derived$bivariate_Count
    }
    if("G" %in% method){
      mif$derived$univariate_NN = derived$univariate_NN
    }
    if("BiG" %in% method){
      mif$derived$bivariate_NN = derived$bivariate_NN
    }
  }
  
  #return resulting mIF object
  return(mif)
}

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
#'@export

ripleys_k = function(mif, mnames, r_range = seq(0, 100, 50),
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
#' @export


bi_ripleys_k <- function(mif,
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
#' x <- NN_G(mif = x, mnames = markers, num_permutations = 1,
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
