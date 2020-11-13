# Dark theme for plotting
theme_dark_mode <- function () { 
  ggdark::dark_mode() %+replace% 
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
}

# positive cells/stroma for simulating IF data 
modes_stroma <- function(k, xmin = 0, xmax = 10, ymin = 0, ymax = 10, sdmin = 1/2, sdmax = 2){ 
  center_peak <- data.frame(x = stats::runif(k,xmin,xmax), 
                            y = stats::runif(k,ymin,ymax),
                            sd_x = stats::runif(k,sdmin,sdmax), 
                            sd_y = stats::runif(k, sdmin, sdmax))
  
  center_peak <- center_peak %>% 
    dplyr::rowwise() %>%
    dplyr::mutate(rho = stats::runif(1, -sd_x * sd_y, sd_x * sd_y))
  
  return(center_peak)
  
}

# Generate random holes for simulating IF data  
hole <- function(xmin = xmin, xmax = xmax, ymin = ymin, 
                ymax = ymax, sdmin = sdmin, sdmax = sdmax){
  #Random select the percent to total area missing
  percent_area_missing <- stats::runif(1, 0.2, 0.35)
  #Random select the number of holes, can range from 1 to floor(percent_area_missing*10)
  num_holes <- sample(2:floor(percent_area_missing*10), 1)
  if(percent_area_missing < 0.2 | num_holes == 1){
    num_holes <- 1
    area_allocation <- percent_area_missing
  }else{
    area_allocation <- Surrogate::RandVec(a = 0.1, 
                                          b = percent_area_missing, 
                                          s = percent_area_missing,
                                          n = num_holes)$RandVecOutput
  }
  
  center_hole <- data.frame(x = stats::runif(num_holes, xmin, xmax),
                            y = stats::runif(num_holes, ymin, ymax),
                            # rho = 0,
                            sd_x = stats::runif(num_holes, sdmin, sdmax), 
                            sd_y = stats::runif(num_holes, sdmin, sdmax))
  center_hole <- center_hole %>%
    dplyr::rowwise() %>%
    dplyr::mutate(rho = stats::runif(1, -sd_x * sd_y, sd_x * sd_y)) %>% 
    dplyr::bind_cols(area = area_allocation[,1]) %>% 
    dplyr::mutate(max_dist = sqrt(area * (xmax - xmin) * (ymax - ymin) / pi))
  
  return(center_hole)
  
}

dist_function <- function(a, other_data){
  diff <- c((pp$x[i] - other_data[a,1]), (pp$y[i] - other_data[a,2]))
  diff <- unlist(diff)
  sigma <- matrix(c(other_data$sd_x[a]^2, rep(other_data$rho[a], 2),
                    other_data$sd_y[a]^2), nrow = 2, ncol = 2)
  z <- t(diff) %*% solve(sigma) %*% diff
  closest <- exp(-z)
}

univariate_ripleys_k <- function(data,
                                 id,
                                 mnames, 
                                 wshape = c("circle", "rectangle"),
                                 r_range = seq(0, 100, 50),
                                 edge_correction = c("theoretical", "translation", "isotropic", "border"),
                                 kestimation = TRUE) {
  
  # estimate_list <- lapply(data, function(data){
  
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
    
  # })
    
}

