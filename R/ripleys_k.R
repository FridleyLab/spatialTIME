#' Calculate Univariate Ripley's K function for IF data
#'
#' @description This function calculates Ripley's K function of IF data to 
#'   characterize correlation of spatial point process using tranlation and
#'   isotropic edge correction method.
#' @param mif An MIF object
#' @param id Character string of variable name for subject ID in TMA data.
#' @param mnames Character vector of marker names to calculate Ripley's K on.
#' @param wshape Character string of window shape. Potenital values are
#'  'rectangle" for rectangular window or "circle" for
#'  circular window. Default is circle.
#' @param r_range Numeric vector of potential r values to estimate K at. 
#' @param csr_calculation Character value indicating the method of calculating Ripley's K
#' @param edge_correction Character value indicating the type of edge correction 
#'  to use. Options include "theoretical", "translation", "isotropic" or "border". 
#'  Various edges corrections are most appropriate in different settings. Default
#'  is "none". 
#' @param kestimation Logical value determining the type estimation performed.
#'  TRUE estimates Ripley's reduced second moment function while FALSE 
#'  estimates Besags's transformation of Ripley's K.
#' @param num_permutations Numeric value indicating the number of permutations used. 
#'  Default is 50.   
#' @param keep_perm_dis Logical value determining whether or not to keep the full 
#'  distribution of permuted K values
#' 
#' @return Returns a list
#'    \item{r}{Subject ID in TMA data}
#'    \item{theo}{Ripley's K estimate using translation edge correction}
#'    
#' @export
#'
ripleys_k <- function(mif,
                      id,
                      mnames, 
                      wshape = "circle",
                      r_range = seq(0, 100, 50),
                      # pick permutation vs theoretical 
                      csr_calculation = "permutation",
                      # permutation number 
                      num_permutations = 50,
                      # edge correction 
                      edge_correction = "translation",
                      # k or l 
                      kestimation = TRUE,
                      keep_perm_dis = FALSE) {
  
  data <- mif[["spatial"]]
  
  # check if any/all provided marker names are not present in the data
  if (all(!mnames %in% colnames(data[[1]]))) {
    stop("No marker names are in the data")
  } else if (any(!mnames %in% colnames(data[[1]]))) {
    stop("Marker names: `", 
         paste(mnames[!mnames %in% colnames(data)], collapse = ", "),
         "` are not in the data")
  }
  
  # check if provided window shape is valid 
  if (!wshape %in% c("circle", "rectangle"))
    stop("invalid window shape name")
  
  # determine calc type
  if (!csr_calculation %in% c("permutation", "theoretical"))
    stop("invalid calculation type")
  
  # determine edge correction
  if (!edge_correction %in% c("none", "translation", "isotropic", "border"))
    stop("invalid edge correction")
  
  # check if provided window shape is valid 
  if  (csr_calculation == "theoretical" & keep_perm_dis == TRUE)
    stop("Permutation distributions not available for theoretical K/L calculations")
  
  # progress bar for k estimation
  pb <- dplyr::progress_estimated(length(data))
  
  if (csr_calculation == "theoretical") {
    
    #   # update progress bar
    #   pb$tick()$print()
      
    estimate_list <- purrr::map(data, univariate_ripleys_k, id, mnames, 
                                wshape, r_range, edge_correction, kestimation) 
    
    estimate_list <- data.frame(matrix(unlist(estimate_list), ncol=4, byrow=T))
    colnames(estimate_list) <- c("sample", "marker", "theoretical_estimate",
                                 "observed_estimate")
    
  } else {
    
    estimate_list <- lapply(data, function(data){

      # update progress bar
      pb$tick()$print()

      perms <- modelr::permute(data, n = num_permutations, mnames)

      perms_df <- lapply(perms$perm, as.data.frame)

      ripleys_estimates <- lapply(perms_df, function(perm_data){
        
        perm_k <- univariate_ripleys_k(perm_data, id, mnames, wshape, r_range,
                                       edge_correction, kestimation)

        return(perm_k)

      })

      # ripleys_estimates <- purrr::map(perms_df, univariate_ripleys_k, id, mnames,
      #                                        wshape, r_range, edge_correction, kestimation)
        
      
      results_list <- ripleys_estimates %>% 
        unlist() %>% 
        matrix(ncol = 4, byrow = TRUE) %>% 
        tibble::as_tibble() 
      
      if (keep_perm_dis == TRUE){
        results_list
      } else {
        results_list <- results_list %>% 
          dplyr::group_by(.data$V1, .data$V2) %>%
          dplyr::rename(sample = .data$V1) %>% 
          dplyr::rename(markers = .data$V2) %>%
          dplyr::summarise(avg_theoretical = mean(as.numeric(.data$V3), na.rm = TRUE),
                           avg_observed = mean(as.numeric(.data$V4), na.rm = TRUE))
      }
      
      return(results_list)
      
    })
    
    estimate_list <- dplyr::bind_rows(estimate_list, .id = "sample")
    
  }
  
  return(estimate_list)
  
}

#' Calculate Bivariate Ripley's K function for IF data
#'
#' @description This function calculates Ripley's K function of IF data for  
#'   two markers.
#' @param mif An MIF object
#' @param id Character string of variable name for subject ID in TMA data.
#' @param mnames A list of character strings containing two marker names
#' @param wshape Character string of window shape. Potenital values are
#'  'rectangle" for rectangular window or "circle" for
#'  circular window. Default is circle.
#' @param r_range Numeric vector of potential r values to estimate K at. 
#' @param csr_calculation Character value indicating the method of calculating Ripley's K. 
#'  Options can be either "theoretical" or "permutation". 
#' @param edge_correction Character value indicating the type of edge correction 
#'  to use. Options include "theoretical", "translation", "isotropic" or "border". 
#'  Various edges corrections are most appropriate in different settings. Default
#'  is "none". 
#' @param kestimation Logical value determining the type estimation performed.
#'  TRUE estimates Ripley's reduced second moment function while FALSE 
#'  estimates Besags's transformation of Ripley's K.
#' @param num_permutations Numeric value indicating the number of permutations used. 
#'  Default is 50.   
#' @param keep_perm_dis Logical value determining whether or not to keep the full 
#'  distribution of permuted K values
#' 
#' @return Returns a list of data frames 
#'    \item{sample}{Subject ID in TMA data}
#'    \item{marker}{Ripley's K estimate using translation edge correction}
#'    \item{theoretical_estimate}{theoretical value of k}
#'    \item{observed_estimate}{observed estimate of k}
#' @export
#'
bi_ripleys_k <- function(mif,
                         id,
                         mnames, 
                         wshape = "circle",
                         r_range = seq(0, 100, 50),
                         # pick permutation vs theoretical 
                         csr_calculation = "permutation",
                         # permutation number 
                         num_permutations = 50,
                         # edge correction 
                         edge_correction = "translation",
                         # k or l 
                         kestimation = TRUE,
                         keep_perm_dis = FALSE) {
  
  data <- mif[["spatial"]]
  
  all_mnames <- unlist(mnames)
  
  mnames <- full_list_combinations(mnames)
  
  # check if any/all provided marker names are not present in the data
  if (all(!all_mnames %in% colnames(data[[1]]))) {
    stop("No marker names are in the data")
  } else if (any(!all_mnames %in% colnames(data[[1]]))) {
    stop("Marker names: `", 
         paste(all_mnames[!all_mnames %in% colnames(data[[1]])], collapse = ", "),
         "` are not in the data")
  }
  
  # check if provided window shape is valid 
  if (!wshape %in% c("circle", "rectangle"))
    stop("invalid window shape name")
  
  # determine calc type
  if (!csr_calculation %in% c("permutation", "theoretical"))
    stop("invalid calculation type")
  
  # determine edge correction
  if (!edge_correction %in% c("none", "translation", "isotropic", "border"))
    stop("invalid edge correction")
  
  # check if provided window shape is valid 
  if  (csr_calculation == "theoretical" & keep_perm_dis == TRUE)
    stop("Permutation distributions not available for theoretical K/L calculations")
  
  # progress bar for k estimation
  pb <- dplyr::progress_estimated(length(data))
  
  if (csr_calculation == "theoretical") {
    
    #   # update progress bar
    #   pb$tick()$print()
    
    estimate_list <- purrr::map(data, bivariate_ripleys_k, id, mnames, 
                                wshape, r_range, edge_correction, kestimation) 
    
    estimate_list <- data.frame(matrix(unlist(estimate_list), ncol=4, byrow=T))
    colnames(estimate_list) <- c("sample", "marker", "theoretical_estimate",
                                 "observed_estimate")
    
  } else {
    
    estimate_list <- lapply(data, function(data){
      
      # update progress bar
      pb$tick()$print()
      
      perms <- modelr::permute(data, n = num_permutations, unlist(mnames))
      
      perms_df <- lapply(perms$perm, as.data.frame)
      
      ripleys_estimates <- lapply(perms_df, function(perm_data){
        
        perm_k <- bivariate_ripleys_k(perm_data, id, mnames, wshape, r_range,
                                       edge_correction, kestimation)
        
        return(perm_k)
        
      })
      
      # ripleys_estimates <- purrr::map(perms_df, univariate_ripleys_k, id, mnames,
      #                                        wshape, r_range, edge_correction, kestimation)
      
      
      results_list <- ripleys_estimates %>% 
        unlist() %>% 
        matrix(ncol = 5, byrow = TRUE) %>% 
        tibble::as_tibble() 
      
      if (keep_perm_dis == TRUE){
        results_list
      } else {
        results_list <- results_list %>% 
          dplyr::group_by(.data$V1, .data$V2) %>%
          dplyr::rename(sample = .data$V1) %>% 
          dplyr::rename(markers = .data$V2) %>%
          dplyr::summarise(avg_theoretical = mean(as.numeric(.data$V3), na.rm = TRUE),
                           avg_observed = mean(as.numeric(.data$V4), na.rm = TRUE))
      }
      
      # results_list <- plyr::ldply(results_list, data.frame)
      
      # results_list <- dplyr::bind_rows(results_list, .id = "V1")
      
      return(results_list)
      
    })
    
    estimate_list <- dplyr::bind_rows(estimate_list, .id = "sample")
    
  }
  
  return(estimate_list)
  
}


