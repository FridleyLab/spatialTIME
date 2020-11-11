#' Calculate Ripley's K function for IF data
#'
#' @description This function calculates Ripley's K function of IF data to 
#'   characterize correlation of spatial point process using tranlation and
#'   isotropic edge correction method.
#' @param mif An MIF object
#' @param id Character string of variable name for subject ID in TMA data.
#' @param mnames Character vector of marker names to calculate Ripley's K on.
#' @param wshape Character string of window shape. Potenital values are
#'  'rectangle" for rectanglular window or "circle" for
#'  circular window. Default is circle.
#' @param r_range Numeric vector of potential r values to estimate K at. 
#' @param edge_correction Character value indicating the type of edge correction 
#'  to use. Options include "theoretical", "translation", "isotropic" or "border". 
#'  Various edges corrections are most appropriate in different settings. Default
#'  is "theroretical". 
#' @param  kestimation Logical value determining the type estimation performed.
#'  TRUE estimates Ripley's reduced second moment function while FALSE 
#'  estimates Besags's transformation of Ripley's K.
#' 
#' @return Returns a data frame
#'    \item{r}{Subject ID in TMA data}
#'    \item{theo}{Ripley's K estimate using translation edge correction}
#'    \item{border}{Ripley's K estimate using isotropic edge correction}
#'    \item{trans}{Theoritical value of Ripley's K}
#'    \item{iso}{Intensity of TMA data}
#'    
#' @export
#'
ripleys_k <- function(mif,
                      id,
                      mnames, 
                      wshape = c("circle", "rectangle"),
                      r_range = seq(0, 100, 50),
                      # pick permutation vs theoretical 
                      # permutation number 
                      # edge correction 
                      edge_correction = c("theoretical", "translation", "isotropic", "border"),
                      # k or l 
                      kestimation = TRUE) {
  
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
  
  # progress bar for k estimation
  pb <- dplyr::progress_estimated(length(data))
  
  estimate_list <- lapply(data, function(data){
  # x and y coordinates for cells
  X <- data %>% 
    janitor::clean_names() %>% 
    dplyr::mutate(xloc = (x_min + x_max) / 2) %>%
    dplyr::mutate(yloc = (y_min + y_max) / 2) %>%
    dplyr::mutate(positive_cell = rowSums(dplyr::select(., !!mnames)) > 0) 
  
  w <- spatstat::convexhull.xy(x = X$xloc, y = X$yloc)
  if (wshape == "circle") {
    w <- spatstat::boundingcircle(w)
  } 
  if(wshape == "rectangle") {
    w <- spatstat::boundingbox(w)
  }
  
  X <- X %>%
    # data with positive marker cell only
    dplyr::filter(positive_cell == TRUE) 
  
  # point pattern object
  p <- spatstat::ppp(x = X$xloc, y = X$yloc, window = w)
  
  # estimate K for variety of distances (r)
  # k_est <- spatstat::Kest(p, r = r_range)
  # we need the function to eventually return K and L estimates 
  if (kestimation == TRUE) {
    est <- spatstat::Kest(p, r = r_range)
  } else {
    est <- spatstat::Lest(p, r = r_range)
  }
  
  # need to figure out dist information from chris 
  if(edge_correction == "theoretical") {
    # possion process - therotical
    k_value <- mean(est$theo) 
  } else if (edge_correction == "isotropic") {
    # isotropic edge correction, good for small number of points
    # if stationary process, trans and iso should be similar
    k_value <- mean(est$iso)  
  } else if (edge_correction == "translation") {
    # translation edge correction, good for small number of points
    k_value <- mean(est$trans)  
  } else {
    k_value <- mean(est$border)  
  }
  
  # intensity 
  int_est <- spatstat::intensity(p)

  results_list <- list(
    `sample_id` = unique(X[[id]]), 
    `estimate` = k_value,
    `intensity` = int_est, 
    `num_points` = sum(X %>% dplyr::pull(positive_cell))
    )
  
  return(results_list)
  
  })
  
}

