#' Calculate Univariate Ripley's K function for IF data
#'
#' @description This function calculates Ripley's K function of IF data to 
#'   characterize correlation of spatial point process using tranlation and
#'   isotropic edge correction method.
#' @param mif An MIF object
#' @param id Character string of variable name for subject ID in TMA data.
#' @param mnames Character vector of marker names to calculate Ripley's K on.
#' @param r_range Numeric vector of potential r values to estimate K at. 
#' @param edge_correction Character value indicating the type of edge correction 
#'  to use. Options include "translation", "isotropic" or "border". 
#'  Various edges corrections are most appropriate in different settings. Default
#'  is "translation". 
#' @param kestimation Logical value determining the type estimation performed.
#'  TRUE estimates Ripley's reduced second moment function while FALSE 
#'  estimates Besags's transformation of Ripley's K.
#' @param num_permutations Numeric value indicating the number of permutations used. 
#'  Default is 50.   
#' @param keep_perm_dis Logical value determining whether or not to keep the full 
#'  distribution of permuted K values
#' @param mlabels Character vector of label for marker names to display in the marker column.
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
                      r_range = seq(0, 100, 50),
                      # permutation number 
                      num_permutations = 50,
                      # edge correction 
                      edge_correction = "translation",
                      # k or l 
                      kestimation = TRUE,
                      keep_perm_dis = FALSE,
                      mlabels = NULL) {
  
  data <- mif[["spatial"]]
  
  # check if any/all provided marker names are not present in the data
  if (all(!mnames %in% colnames(data[[1]]))) {
    stop("No marker names are in the data")
  } else if (any(!mnames %in% colnames(data[[1]]))) {
    stop("Marker names: `", 
         paste(mnames[!mnames %in% colnames(data)], collapse = ", "),
         "` are not in the data")
  }
  
  # determine calc type
  # if (!csr_calculation %in% c("permutation", "observed"))
  #   stop("invalid calculation type")
  
  # determine edge correction
  if (!edge_correction %in% c("none", "translation", "isotropic", "border"))
    stop("invalid edge correction")
  
  # # check if provided window shape is valid 
  # if  (csr_calculation == "observed" & keep_perm_dis == TRUE)
  #   stop("Permutation distributions not available for observed K/L calculations")
  # 
  # progress bar for k estimation
  pb <- dplyr::progress_estimated(length(data))
  
  # if (csr_calculation == "observed") {
  #   
  #   #   # update progress bar
  #   #   pb$tick()$print()
  #     
  #   estimate_list <- purrr::map(data, univariate_ripleys_k, id, mnames, 
  #                               r_range, edge_correction, kestimation, mlabels) 
  #   
  #   estimate_list <- dplyr::bind_rows(estimate_list)  %>% 
  #     dplyr::rename('Theoretical CSR' = .data$csr_theoretical,
  #                    'Observed K' = .data$observed_estimate) %>%
  #     mutate('Degree of Clustering' = `Observed K` - `Theoretical CSR`)
  #   
  # } else {
    
    estimate_list <- lapply(data, function(data){

      # update progress bar
      pb$tick()$print()

      perms <- modelr::permute(data, n = num_permutations, mnames)

      perms_df <- lapply(perms$perm, as.data.frame)

      ripleys_estimates <- lapply(perms_df, function(perm_data){
        
        perm_k <- univariate_ripleys_k(perm_data, id, mnames, r_range,
                                       edge_correction, kestimation, mlabels)
        
        perm_k <- dplyr::bind_rows(perm_k)

        return(perm_k)

      })
      
      results_list <- dplyr::bind_rows(ripleys_estimates)
      
      if (keep_perm_dis == TRUE){
        results_list<- results_list %>% 
          dplyr::rename('Permuted CSR' = .data$observed_estimate,
                        'Theoretical CSR' = .data$csr_theoretical) 
        
      } else {
        results_list <- results_list %>% 
          dplyr::rename(csr_permuted = .data$observed_estimate) %>%
          dplyr::group_by(.data[[id]], .data$marker, .data$r_value) %>%
          dplyr::summarise("Permuted CSR" = mean(as.numeric(.data$csr_permuted),
                                                   na.rm = TRUE),
                           "Theoretical CSR" = mean(as.numeric(.data$csr_theoretical),
                                                      na.rm = TRUE)) 
      }
      
      return(results_list)
      
    })
    
    #The commented version replaced the sample column with 1,2,3,....
    #estimate_list <- dplyr::bind_rows(estimate_list, .id = "sample")
    estimate_list <- dplyr::bind_rows(estimate_list)
    observed_list <- purrr::map(data, univariate_ripleys_k, id, mnames,
                                r_range, edge_correction, kestimation, mlabels)
    observed_list <- dplyr::bind_rows(observed_list)
    
    estimate_list <- estimate_list %>%
      dplyr::left_join(observed_list %>%
                         dplyr::select(.data[[id]], .data$marker,
                                       .data$r_value, .data$observed_estimate) %>%
                         dplyr::rename(`Observed K` = .data$observed_estimate),
                       by = c(id, "marker", "r_value")) %>%
      dplyr::mutate(`Degree of Clustering` = .data$`Observed K` - .data$`Permuted CSR`)
    
  #}
  
  return(estimate_list)
  
}

#' Calculate Bivariate Ripley's K function for IF data
#'
#' @description This function calculates Ripley's K function of IF data for  
#'   two markers.
#' @param mif An MIF object
#' @param id Character string of variable name for subject ID in TMA data.
#' @param mnames A list of character strings containing two marker names
#' @param r_range Numeric vector of potential r values to estimate K at. 
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
#' @param mlabels A list of character strings containing two marker labels
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
                         r_range = seq(0, 100, 50),
                         # permutation number 
                         num_permutations = 50,
                         # edge correction 
                         edge_correction = "translation",
                         # k or l 
                         kestimation = TRUE,
                         keep_perm_dis = FALSE,
                         mlabels = NULL) {
  
  data <- mif[["spatial"]]
  
  all_mnames <- unlist(mnames)
  
  mnames <- full_list_combinations(mnames)
  
  if(is.null(mlabels)){
    mlabels = mnames
  }
  
  all_mlabels <- unlist(mlabels)
  
  mlabels <- full_list_combinations(mlabels)
  
  # check if any/all provided marker names are not present in the data
  if (all(!all_mnames %in% colnames(data[[1]]))) {
    stop("No marker names are in the data")
  } else if (any(!all_mnames %in% colnames(data[[1]]))) {
    stop("Marker names: `", 
         paste(all_mnames[!all_mnames %in% colnames(data[[1]])], collapse = ", "),
         "` are not in the data")
  }
  
  # # determine calc type
  # if (!csr_calculation %in% c("permutation", "observed"))
  #   stop("invalid calculation type")
  # 
  # determine edge correction
  if (!edge_correction %in% c("none", "translation", "isotropic", "border"))
    stop("invalid edge correction")
  
  # # check if provided window shape is valid 
  # if  (csr_calculation == "observed" & keep_perm_dis == TRUE)
  #   stop("Permutation distributions not available for observed K/L calculations")
  # 
  # progress bar for k estimation
  pb <- dplyr::progress_estimated(length(data))
  
  # if (csr_calculation == "observed") {
  #   
  #   estimate_list <- purrr::map(data, bivariate_ripleys_k, id, mnames, 
  #                               r_range, edge_correction, kestimation) 
  #   
  #   estimate_list <- dplyr::bind_rows(estimate_list)
  #   
  #   
  # } else {
  #   
    estimate_list <- lapply(data, function(data){
      
      # update progress bar
      pb$tick()$print()
      
      perms <- modelr::permute(data, n = num_permutations, unique(unlist(mnames)))
      
      perms_df <- lapply(perms$perm, as.data.frame)
      
      ripleys_estimates <- lapply(perms_df, function(perm_data){
        
        perm_k <- bivariate_ripleys_k(perm_data, id, mnames, r_range,
                                       edge_correction, kestimation, mlabels)
        
        perm_k <- dplyr::bind_rows(perm_k)
        
        return(perm_k)
        
      })
      
      results_list <- dplyr::bind_rows(ripleys_estimates)
      
      if (keep_perm_dis == TRUE){
        results_list <- results_list %>% 
          dplyr::rename(`Permuted CSR` = .data$observed_estimate,
                        `Theoretical CSR` = .data$csr_theoretical) 
      } else {
        results_list <- results_list %>% 
          dplyr::rename(csr_permuted = .data$observed_estimate) %>%
          dplyr::group_by(.data[[id]], .data$anchor_marker, .data$comparison_marker,
                          .data$r_value) %>%
          dplyr::summarise(`Permuted CSR` = mean(as.numeric(.data$csr_permuted),
                                                 na.rm = TRUE),
                           `Theoretical CSR` = mean(as.numeric(.data$csr_theoretical),
                                                    na.rm = TRUE))
      }
      # results_list <- plyr::ldply(results_list, data.frame)
      
      # results_list <- dplyr::bind_rows(results_list, .id = "V1")
      
      return(results_list)
      
    })
    
    estimate_list <- dplyr::bind_rows(estimate_list)
    
    observed_list <- purrr::map(data, bivariate_ripleys_k, id, mnames, 
                                r_range, edge_correction, kestimation,
                                mlabels) 
    observed_list <- dplyr::bind_rows(observed_list)
    
    estimate_list <- estimate_list %>%
      dplyr::left_join(observed_list %>% 
                         dplyr::select(.data[[id]], .data$anchor_marker,
                                       .data$comparison_marker, .data$r_value,
                                       .data$observed_estimate) %>% 
                         dplyr::rename(`Observed K` = .data$observed_estimate),
                       by = c(id, "anchor_marker", "comparison_marker", "r_value") )%>%
      dplyr::mutate(`Degree of Clustering` = .data$`Observed K` - .data$`Permuted CSR`)
    
  #}
  
  return(estimate_list)
  
}


