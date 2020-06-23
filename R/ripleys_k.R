#' Calculate Ripley's K function for IF data
#' 
#' @description This function calculates Ripley's K function of IF data to 
#'   characterize correlation of spatial point process using tranlation and
#'   isotropic edge correction method.
#' @param data A TMA data frame.
#' @param id Character string of variable name for subject ID in TMA data.
#' @param mnames Character vector of marker names to calculate Ripley's K on.
#' @param wshape Character string of window shape. Potenital values are
#'  "r" or "rectangle" for rectanglular window and "c" or "circle" for
#'  circular window.
#' @param dist Distance at which to calculate Ripley's K
#' @param inside Buffer around circular window 
#' 
#' @return Returns a list of components
#'    \item{Subject ID}{Subject ID in TMA data}
#'    \item{Translation K}{Ripley's K estimate using translation edge correction}
#'    \item{Isotropic K}{Ripley's K estimate using isotropic edge correction}
#'    \item{Theoritical K}{Theoritical value of Ripley's K}
#'    \item{Intensity estimate}{Intensity of TMA data}
#'    \item{N points}{Number of cells where the markers are positive}
#'    
#' @export
#' 
#' @examples
#' ripleys_k(if_data[[1]], .id = "subid.x", mnames = marker_names, wshape = "r")
#'
ripleys_k <- function(data, id, mnames, wshape, dist = 200, inside = 75) {
  # check if any/all provided marker names are not present in the data
  if (all(!mnames %in% colnames(data))) {
    stop("No marker names are in the data")
  } else if (any(!mnames %in% colnames(data))) {
    stop("Marker names: `", 
         paste(mnames[!mnames %in% colnames(data)], collapse = ", "),
         "` are not in the data")
  }
  
  # check if provided window shape is valid 
  if (!wshape %in% c("circle", "c", "rectangle", "r"))
    stop("invalid window shape name")
  
  ripleys_data <- data %>%
    dplyr::mutate(xloc = (XMin + XMax) / 2) %>% 
    dplyr::mutate(yloc = (YMin + YMax) / 2) %>%
    # one positive_cell per marker or overall - currently overall
    # dplyr::mutate(dplyr::across(dplyr::contains(mnames),
    #                             .fns = list(positive_cell = 
    #                                           ~rowSums(.))) %>% 
    dplyr::mutate(positive_cell = rowSums(select(., !!mnames)) > 0) %>% 
    dplyr::mutate(w = )
  
  # x and y coordinates for cells
  X <- data %>% 
    dplyr::mutate(
      xloc = (XMin + XMax) / 2,
      yloc = (YMin + YMax) / 2,
      positive_cell = rowSums(select(., !!mnames)) > 0) 
  
  # window
  if (wshape %in% c("circle", "c")) {
    # radius for x and y direction
    xrad <- (max(X[["xloc"]]) - min(X[["xloc"]])) / 2 
    yrad <- (max(X[["yloc"]]) - min(X[["yloc"]])) / 2
    
    # taking radius to make circle to be max radius
    rad <- max(xrad, yrad) + inside
    
    # center coordinate of image
    xcen <- min(X[["xloc"]]) + xrad
    ycen <- min(X[["yloc"]]) + yrad
    
    # window
    w <- spatstat::disc(radius = rad, centre = c(xcen, ycen))
  } else {
    # range for rectangular plot
    xrange <- range(X[["xloc"]])
    yrange <- range(X[["yloc"]])
    
    # data with positive marker cell only
    Xpos <- X[X[["positive_cell"]], ]
    
    # window
    w <- spatstat::owin(xrange, yrange)
  }
  
  # data with positive marker cell only
  Xpos <- X[X[["positive_cell"]], ]
  
  # point pattern object
  p <- spatstat::ppp(x = Xpos$xloc, y = Xpos$yloc, window = w)
  
  # estimate K for variety of distances (r)
  k_est <- spatstat::Kest(p)
  
  # border edge correction, not good for small number of points
  border <- mean(k_est$border[round(k_est$r) == dist]) 
  
  # translation edge correction, good for small number of points
  translate <- mean(k_est$trans[round(k_est$r) == dist])  
  
  # isotropic edge correction, good for small number of points
  # if stationary process, trans and iso should be similar
  isotropic <- mean(k_est$iso[round(k_est$r) == dist])  
  
  # possion process - therotical
  theo <- mean(k_est$theo[round(k_est$r) == dist])  
  
  # intensity 
  int_est <- spatstat::intensity(p)
  
  list(
    `Subject ID` = unique(X[[.id]]), `Translation K` = translate,
    `Isotropic K` = isotropic, `Theoritical K` = theo,
    `Intensity estimate` = int_est, `N points` = sum(X[["positive_cell"]]))
}

