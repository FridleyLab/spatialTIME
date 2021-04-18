#' @importFrom ggplot2 %+replace%
#' @importFrom plyr .

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
    dplyr::mutate(rho = stats::runif(1, -.data$sd_x * .data$sd_y, .data$sd_x * .data$sd_y))
  
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
    dplyr::mutate(rho = stats::runif(1, -.data$sd_x * .data$sd_y, .data$sd_x * .data$sd_y)) %>% 
    dplyr::bind_cols(area = area_allocation[,1]) %>% 
    dplyr::mutate(max_dist = sqrt(.data$area * (.data$xmax - .data$xmin) * (.data$ymax - .data$ymin) / pi))
  
  return(center_hole)
  
}

# dist_function <- function(a, other_data){
#   diff <- c((pp$x[i] - other_data[a,1]), (pp$y[i] - other_data[a,2]))
#   diff <- unlist(diff)
#   sigma <- matrix(c(other_data$sd_x[a]^2, rep(other_data$rho[a], 2),
#                     other_data$sd_y[a]^2), nrow = 2, ncol = 2)
#   z <- t(diff) %*% solve(sigma) %*% diff
#   closest <- exp(-z)
# }

univariate_ripleys_k <- function(data,
                                 id,
                                 mnames, 
                                 r_range = seq(0, 100, 50),
                                 edge_correction = c("translation", "isotropic", "border"),
                                 kestimation = TRUE,  mlabels=NULL) {
  if(is.null(mlabels)){
    mlabels = mnames
    }
  mlabels = stats::setNames(object = mlabels, nm = mnames)
    estimate_list <- lapply(mnames, function(mnames){
      mnames_clean = janitor::make_clean_names(mnames)
      mlabels = mlabels[mnames]
    # x and y coordinates for cells
    ## Jordan why do the names need to be cleaned? I see
    ## Is there away to only change a couple of columns (i.e. leave out mnames)
    X <- data %>% 
      janitor::clean_names() %>% 
      dplyr::mutate(xloc = (.data$x_min + .data$x_max) / 2) %>%
      dplyr::mutate(yloc = (.data$y_min + .data$y_max) / 2) %>%
      dplyr::mutate(positive_cell = rowSums(dplyr::select(., !!mnames_clean)) > 0) 
    
    sample_name <- data %>% 
      dplyr::slice(1) %>% 
      dplyr::pull(!!id)
    
    #Something changed this function is no longer the spatstat namespace
    w <- spatstat.geom::convexhull.xy(x = X$xloc, y = X$yloc)
    
    X <- X %>%
      # data with positive marker cell only
      dplyr::filter(.data$positive_cell == TRUE) 
    
    # double check with chris that the function should return NA 
    ## How do you find the column that the sample and spatial columns will be 
    ## merged by? sample should have the same column name as that merge variable
    ## To do any downstream analysis this obj will be merged with clinical and summary data
    if (nrow(X) <= 1) {
      results_list <- data.frame(
        id = sample_name,
        marker = unname(mlabels),
        r_value = r_range,
        observed_estimate = NA,
        csr_theoretical = NA,
        degree_of_spatial_diff = NA 
      )
    } else {
      
      # point pattern object
      p <- suppressWarnings(spatstat.geom::ppp(x = X$xloc, y = X$yloc, window = w))

      # we need the function to eventually return K and L estimates 
      if (kestimation == TRUE) {
        est <- spatstat.core::Kest(p, r = r_range, correction="all")
      } else {
        est <- spatstat.core::Lest(p, r = r_range, correction="all")
      }
      
      if (edge_correction == "isotropic") {
        # isotropic edge correction, good for small number of points
        # if stationary process, trans and iso should be similar
        # k_value <- mean(est$iso)  
        k_value <- est$iso  
      } else if (edge_correction == "translation") {
        # translation edge correction, good for small number of points
        # k_value <- mean(est$trans)  
        k_value <- est$trans
      } else {
        # k_value <- mean(est$border) 
        k_value <- est$border
      }
      
      results_list <- data.frame(
        id = sample_name,
        marker = unname(mlabels),
        r_value = r_range,
        observed_estimate = k_value,
        csr_theoretical = est$theo
      )
      
    }
    colnames(results_list)[1] = id 
    return(results_list)
    
  })
}

bivariate_ripleys_k <- function(data,
                                id,
                                mnames, 
                                r_range = seq(0, 100, 50),
                                edge_correction = c("translation", "isotropic", "border"),
                                kestimation = TRUE, mlabels = NULL) {
  if(is.null(mlabels)){
    mlabels = mnames
  }
  counter = 0
  estimate_list <- lapply(mnames, function(mnames){
    counter <<- counter + 1
    mnames_clean = lapply(mnames,janitor::make_clean_names)
    # x and y coordinates for cells
    X <- data %>% 
      janitor::clean_names() %>% 
      dplyr::mutate(xloc = (.data$x_min + .data$x_max) / 2) %>%
      dplyr::mutate(yloc = (.data$y_min + .data$y_max) / 2) %>%
      dplyr::mutate(type1_cell = rowSums(dplyr::select(., .data[[mnames_clean[[1]]]]) > 0)) %>% 
      dplyr::mutate(type2_cell = rowSums(dplyr::select(., .data[[mnames_clean[[2]]]]) > 0)) %>% 
      dplyr::mutate(overall_type = dplyr::case_when(
        .data$type1_cell == 1 ~ "type_one",
        .data$type2_cell == 1 ~ "type_two",
        TRUE ~ "neither"
      ))
    
    if(nrow(X[X$type1_cell == 1,]) == 0){
      warning("No cells positive for ", mlabels[[counter]][[1]],
              " were found - returning NA")
      
      results_list <- data.frame(
        sample = X[[id]][1],
        anchor_marker = mlabels[[counter]][[1]],
        comparison_marker = mlabels[[counter]][[2]],
        r_value = r_range,
        csr_theoretical = NA,
        observed_estimate = NA 
      )
      colnames(results_list)[1] = id
      return(results_list)
    }
    
    if(nrow(X[X$type2_cell == 1,]) == 0){
      warning("No cells positive for ", mlabels[[counter]][[2]],
              " were found - returning NA")
      
      results_list <- data.frame(
        sample = X[[id]][1],
        anchor_marker = mlabels[[counter]][[1]],
        comparison_marker = mlabels[[counter]][[2]],
        r_value = r_range,
        csr_theoretical = NA,
        observed_estimate = NA 
      )
      colnames(results_list)[1] = id
      return(results_list)
    }
    
    if(nrow(X[X$type1_cell == 1 & X$type2_cell ==1,]) >= 1){
      
      warning(paste0(nrow(X[X$type1_cell == 1 & X$type2_cell ==1,]), 
                     " cells removed due to being positive for both ",
                     mlabels[[counter]][[1]], " and ", mlabels[[counter]][[2]]))
      
      X <- X %>% 
        dplyr::filter(!(.data$type1_cell == 1 & .data$type2_cell ==1))
      
    }
    
    sample_name <- X %>% 
      dplyr::slice(1) %>% 
      dplyr::pull(!!id)
    
    w <- spatstat.geom::convexhull.xy(x = X$xloc, y = X$yloc)
    if (nrow(X) <= 1 | nrow(X[X$overall_type=="type_one",]) <=1 |
                            nrow(X[X$overall_type=="type_two",]) <=1 ) {
      
      results_list <- data.frame(
        sample = sample_name,
        anchor_marker = mlabels[[counter]][[1]],
        comparison_marker = mlabels[[counter]][[2]],
        r_value = r_range,
        csr_theoretical = NA,
        observed_estimate = NA 
      )
      
    } else {
      
      X <- X %>% 
        dplyr::filter(.data$overall_type != "neither") %>% 
        dplyr::select(.data$xloc, .data$yloc, .data$overall_type)
      
      #Make a marked point process with the window from above,
      #cell locations, and a marks as defined by cell type
      pp_cross = spatstat.geom::ppp(x = X$xloc, y = X$yloc, 
                               window = w, marks = factor(X$overall_type)) 
      
      
      if (kestimation == TRUE) {
        suppressWarnings({est <- spatstat.core::Kcross(X = pp_cross, i = "type_one", j = "type_two", r = r_range)})
      } else {
        suppressWarnings({est <- spatstat.core::Lcross(pp_cross, i = "type_one", j = "type_two", r = r_range)})
      }
      
      if (edge_correction == "isotropic") {
        # k_value <- mean(est$iso)  
        k_value <- est$iso
      } else if (edge_correction == "translation") {
        # k_value <- mean(est$trans)  
        k_value <- est$trans
      } else {
        # k_value <- mean(est$border)  
        k_value <- est$border 
      }
      
      results_list <- data.frame(
        sample = sample_name,
        anchor_marker = mlabels[[counter]][[1]],
        comparison_marker = mlabels[[counter]][[2]],
        r_value = r_range,
        csr_theoretical = est$theo,
        observed_estimate = k_value 
      )
    }
    
    colnames(results_list)[1] = id
    return(results_list)
    
  })
  
}

full_list_combinations <- function(list_pairs){
  
  expanded_list <- list_pairs
  
  for (i in 1:length(list_pairs)){
    
    expanded_list[[length(list_pairs) + i]] <- list(unlist(list_pairs[[i]][2]),
                                                    unlist(list_pairs[[i]][1]))
    # expanded_list[length(list_pairs) + i] <- list(list_pairs[[i]][2], list_pairs[[i]][1])
    
  }
  
  return(expanded_list)
}
