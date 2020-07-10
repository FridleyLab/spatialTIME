#' Calculate Ripley's K function for IF data
#'
#' @description This function calculates Ripley's K function of IF data to 
#'   characterize correlation of spatial point process using tranlation and
#'   isotropic edge correction method.
#' @param data A TMA data frame.
#' @param id Character string of variable name for subject ID in TMA data.
#' @param mnames Character vector of marker names to calculate Ripley's K on.
#' @param wshape Character string of window shape. Potenital values are
#'  "r" or "rectangle" for rectanglular window, "c" or "circle" for
#'  circular window and "i" or "irregular" for a custom polygon.
#' @param r_range Numeric vector of potential r values to estimate K at. 
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
#' @examples
#' ripleys_k(if_data[[1]], .id = "subid.x", mnames = marker_names)
#'
ripleys_k <- function(data, id, mnames, 
                      wshape = "irregular", r_range = seq(0, 100, 50)) {
  # check if any/all provided marker names are not present in the data
  if (all(!mnames %in% colnames(data))) {
    stop("No marker names are in the data")
  } else if (any(!mnames %in% colnames(data))) {
    stop("Marker names: `", 
         paste(mnames[!mnames %in% colnames(data)], collapse = ", "),
         "` are not in the data")
  }
  
  # check if provided window shape is valid 
  if (!wshape %in% c("circle", "c", "rectangle", "r", "i", "irregular"))
    stop("invalid window shape name")
  
  # x and y coordinates for cells
  X <- data %>% 
    dplyr::mutate(xloc = (XMin + XMax) / 2) %>%
    dplyr::mutate(yloc = (YMin + YMax) / 2) %>%
    dplyr::mutate(positive_cell = rowSums(dplyr::select(., !!mnames)) > 0) 
  
  w <- spatstat::convexhull.xy(x = X$xloc, y = X$yloc)
  if (wshape %in% c("circle", "c")) {
    w <- boundingcircle(w)
  } 
  if(wshape %in% c("rectangle", "r")) {
    w <- boundingbox(w)
  }
  
  X <- X %>%
    # data with positive marker cell only
    dplyr::filter(positive_cell == TRUE) 
  
  # point pattern object
  p <- spatstat::ppp(x = X$xloc, y = X$yloc, window = w)
  
  # estimate K for variety of distances (r)
  # k_est <- spatstat::Kest(p, r = r_range)
  # we need the function to eventually return K and L estimates 
  l_est <- spatstat::Lest(p, r = r_range)
  return(l_est)
  
  # # border edge correction, not good for small number of points
  # border <- mean(k_est$border[round(k_est$r) == dist]) 
  # 
  # # translation edge correction, good for small number of points
  # translate <- mean(k_est$trans[round(k_est$r) == dist])  
  # 
  # # isotropic edge correction, good for small number of points
  # # if stationary process, trans and iso should be similar
  # isotropic <- mean(k_est$iso[round(k_est$r) == dist])  
  # 
  # # possion process - therotical
  # theo <- mean(k_est$theo[round(k_est$r) == dist])  
  # 
  # # intensity 
  # int_est <- spatstat::intensity(p)
  # 
  # list(
  #   `Subject ID` = unique(X[[.id]]), `Translation K` = translate,
  #   `Isotropic K` = isotropic, `Theoritical K` = theo,
  #   `Intensity estimate` = int_est, `N points` = sum(X[["positive_cell"]]))
}

