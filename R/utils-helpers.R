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
      
    if(nrow(data) > 0){s
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
    
  }else{
    results_list <- data.frame(
      id = NA,
      marker = unname(mlabels),
      r_value = r_range,
      observed_estimate = NA,
      csr_theoretical = NA,
      degree_of_spatial_diff = NA 
    )
    return(results_list)
  }
  })
    return(estimate_list)
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
        sample = data[[id]][1],
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
        sample = data[[id]][1],
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
    
       
    w <- spatstat.geom::convexhull.xy(x = X$xloc, y = X$yloc)
    if (nrow(X) <= 1 | nrow(X[X$overall_type=="type_one",]) <=1 |
                            nrow(X[X$overall_type=="type_two",]) <=1 ) {
      
      results_list <- data.frame(
        sample = data[[id]][1],
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
        sample = data[[id]][1],
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

#################################################################################
#Updates Chris 7/8/2021

K_out = function(data, marker, id, iter, correction,r_value, win){
  #Does the actual computation of Ripley' K
  K_obs = spatstat.geom::ppp(x = data$xloc, y = data$yloc, window = win) %>%
    spatstat.core::Kest(r = r_value, correction = correction) %>%
    data.frame() %>%
    dplyr::filter(r != 0) %>%
    dplyr::mutate(Marker = marker,
           label = data[[id]][1],
           iter = iter) %>%
    dplyr::select(iter, label, Marker, r, theo, correction) %>%
    dplyr::mutate(label = as.character(label))
  return(K_obs)
}

perm_data = function(data, markers){
  #Generate data uses permutation for all markers
  data_pos = data %>%
    dplyr::group_by(Marker) %>%
    dplyr::mutate(Positive = sample(Positive, nrow(data)/length(markers), replace = FALSE))
  return(data_pos)
}

## Univariate K Helper functions

uni_K = function(data = data, iter, marker, id, correction, r_value, win){
  #Filters down to the positive cells and then computes Ripley's K
  data_pos = data %>% 
    dplyr::filter(Marker == marker, Positive == 1)
  
  if(nrow(data_pos)==0){
    K = data.frame(iter = iter,label = data[[id]][1],
                   Marker = marker, r = r_value, 
                   theo = NA, correction = NA) %>%
      dplyr::filter(r > 0) %>%
      dplyr::mutate(label = as.character(label))
    colnames(K)[6] = correction
  }else{
    K = K_out(data = data_pos, marker = marker, 
              id = id, iter = iter, correction = correction,
              r_value = r_value, win = win)
  }
  return(K)
}


uni_Rip_K = function(data, markers, id, num_iters, correction = 'trans', method = 'K', 
                     perm_dist, r){
  #Main function
  #Check if the method is selected from K, L, M
  if(!(method %in% c('K',"L","M"))){
    stop("Did not provide a valid method.")
  }
  
  #Notice that this follows spatstat's notation and argument name
  if(!(correction %in% c('trans', "iso", "border", 'translation', 'isotropic'))){
    stop("Did not provide a valid edge correcion method.")
  }
  
  #These next two if statements deal with the fact that spatstat::Kest accepts,
  #multiple spelling of the edge correction, but they not necessarily match the
  #column name in the Kest output.
  if(correction == 'translation'){
    correction = 'trans'
  }
  
  if(correction == 'isotropic'){
    correction = 'iso'
  }
  
  #Use set the cell location as the center of the cell
  data = data %>% 
    dplyr::mutate(xloc = (XMin + XMax)/2,
           yloc = (YMin + YMax)/2
    )
  
  #Create the region that the point process exists. This only needs to be done
  #once per image
  win = spatstat.geom::convexhull.xy(x = data$xloc, y = data$yloc)
  
  #Make the data a long format
  data = data %>%
    tidyr::pivot_longer(cols = all_of(markers), names_to = 'Marker', values_to = 'Positive')
  
  grid = expand.grid(markers, 1:num_iters) %>%
    dplyr::mutate(Var1 = as.character(Var1))
  
  perms = purrr::map_df(.x = 1:nrow(grid),~{
    data_new = perm_data(data, markers)
    return(uni_K(data = data_new, iter = grid[.x,2], 
                 marker = grid[.x,1], id = id, r_value = r,
                 correction = correction, win = win))})
  colnames(perms)[c(2,5,6)] = c(id, 'Theoretical CSR','Permuted K')
  
  obs = purrr::map_df(.x = markers, ~uni_K(data = data, iter = 'Observed', 
                                    marker = .x, correction = correction,
                                    id = id, r_value = r,
                                    win = win)) %>%
    dplyr::select(-iter) 
  colnames(obs)[c(1,4,5)] = c(id, 'Theoretical CSR', 'Observed K')
  
  final = suppressMessages(dplyr::left_join(perms, obs))
  
  if(method == 'M'){
    final = final %>% dplyr::mutate_at(c("Theoretical CSR","Permuted K","Observed K"),
                                ~./(pi*r^2)) %>%
      dplyr::mutate('Degree of Clustering Permutation' = `Permuted K`,
             'Degree of Clustering Theoretical' = `Theoretical CSR`) %>%
      dplyr::rename('Permuted M' = `Permuted K`,
             'Theoretical CSR' = `Theoretical CSR`,
             'Observed M' = `Observed K`)
  }
  
  #Conduct the necessary transformation
  if(method == 'L'){
    final = final %>% dplyr::mutate_at(c("Theoretical CSR","Permuted K","Observed K"),
                                ~sqrt(./pi)) %>%
      dplyr::mutate('Degree of Clustering Permutation' = `Observed K` - `Permuted K`,
             'Degree of Clustering Theoretical' = `Observed K` - `Theoretical CSR`)  %>%
      dplyr::rename('Permuted L' = `Permuted K`,
             'Theoretical CSR' = `Theoretical CSR`,
             'Observed L' = `Observed K`)
  }
  
  if(method == 'K'){
    final = final %>% 
      dplyr::mutate('Degree of Clustering Permutation' = `Observed K` - `Permuted K`,
             'Degree of Clustering Theoretical' = `Observed K` - `Theoretical CSR`)
  }
  
  if(!perm_dist){
    final = final %>% 
      dplyr::mutate(id = get(id)) %>%
      dplyr::select(-(1:2)) %>%
      dplyr::group_by(id,Marker, r) %>%
      #dplyr::summarize(`Theoretical CSR` = mean(`Theoretical CSR`, na.rm = TRUE),
      #          `Permuted CSR` = mean(.[[grep('Permuted', colnames(.), value = TRUE)]], na.rm = TRUE),
      #          `Observed` = mean(.[[grep('Observed', colnames(.), value = TRUE)]], na.rm = TRUE),
      #          `Degree of Clustering Theoretical` = mean(`Degree of Clustering Theoretical`, na.rm = TRUE),
      #          `Degree of Clustering Permutation` =  mean(`Degree of Clustering Permutation`, na.rm = TRUE)
      #          )
      dplyr::summarize_all(~mean(., na.rm = TRUE))
  }
  colnames(final)[1] = id
  return(final)
}

## Bivariate Helper Functions
bi_K = function(data, mark_pair, r, correction, id, iter, win){
  mark_pair = unname(mark_pair)
  #Keeps all positive cell for either marker, and then discards the repeated locations
  data_new = data %>% 
    dplyr::filter(Marker %in% c(mark_pair[1], mark_pair[2]),
           Positive == 1) %>%
    dplyr::group_by(xloc, yloc) %>%
    dplyr::filter(dplyr::n() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Marker = factor(Marker, levels = c(mark_pair[1], mark_pair[2])))
  
  
  if(any(table(data_new$Marker) == 0)){
    #weird error here :row names were found from a short variable and have been discarded
    K = data.frame(r = r, theo = NA, trans = NA, 
                   anchor = levels(data_new$Marker)[1], 
                   counted = levels(data_new$Marker)[2],
                   label = data[[id]][1], 
                   iter = iter)  %>% 
      dplyr::select(iter, label, anchor, counted,  r, theo, correction) %>% #hardcoded
      dplyr::filter(r>0)
  }else{
    X = spatstat.geom::ppp(x = data_new$xloc, y = data_new$yloc, window = win, marks = data_new$Marker)
    K = spatstat.core::Kcross(r =r, X = X, i = levels(data_new$Marker)[1], 
               j = levels(data_new$Marker)[2], correction = 'translation') %>% #Hard coded
      data.frame() %>%
      dplyr::filter(r!=0) %>%
      dplyr::mutate(anchor = levels(data_new$Marker)[1],
             counted = levels(data_new$Marker)[2],
             iter = iter,
             label = data[[id]][1]) %>% 
      dplyr::select(iter, label, anchor, counted, r, theo, correction) #Hardcoded
  }
  return(K)
}

bi_Rip_K = function(data, markers, id, num_iters, correction = 'trans', 
                    method, perm_dist, r, exhaustive){
  #Main function
  #Check if the method is selected from K, L, M
  if(!(method %in% c('K',"L","M"))){
    stop("Did not provide a valid method.")
  }
  
  #Check if exhaustive is FALSE, then markers must be a data.frame
  if(exhaustive == FALSE & class(markers) != 'data.frame'){
    stop("If exhaustive == FALSE, then markers must be a data.frame.")
  }
  if(exhaustive == TRUE & class(markers) != 'character'){
    stop("If exhaustive == TRUE, then markers must be a character vector")
  }
  
  #Notice that this follows spatstat's notation and argument name
  if(!(correction %in% c('trans', "iso", "border", 'translation', 'isotropic'))){
    stop("Did not provide a valid edge correcion method.")
  }
  
  #These next two if statements deal with the fact that spatstat::Kest accepts,
  #multiple spelling of the edge correction, but they not necessarily match the
  #column name in the Kest output.
  if(correction == 'translation'){
    correction = 'trans'
  }
  
  if(correction == 'isotropic'){
    correction = 'iso'
  }
  
  #Use set the cell location as the center of the cell
  data = data %>% 
    dplyr::mutate(xloc = (XMin + XMax)/2,
           yloc = (YMin + YMax)/2
    )
  
  #Create the region that the point process exists. This only needs to be done
  #once per image
  win = spatstat.geom::convexhull.xy(x = data$xloc, y = data$yloc)
  
  #Make the data a long format
  
  #Enumerates all possible combination of markers, and removes the ones where
  #marker 1 and marker 2 are the same
  if(exhaustive){
    data = data %>%
      tidyr::pivot_longer(cols = all_of(markers), 
                          names_to = 'Marker', values_to = 'Positive')
    grid = expand.grid(markers, markers, 1:num_iters) %>%
      dplyr::mutate(Var1 = as.character(Var1),
           Var2 = as.character(Var2)) %>%
      dplyr::filter(Var1 != Var2)
  }else{
    grid = expand_grid(markers, 1:num_iters) %>%
      data.frame() %>%
      dplyr::mutate(Var1 = as.character(Var1),
             Var2 = as.character(Var2)) %>%
      dplyr::filter(Var1 != Var2)
    markers = grid %>% dplyr::select(Var1, Var2) %>%
      unlist() %>% unique()
    data = data %>%
      tidyr::pivot_longer(cols = all_of(markers), 
                          names_to = 'Marker', values_to = 'Positive')
  }
  
  perm = purrr::map_df(.x = 1:nrow(grid), 
                ~{
                  data_new = perm_data(data, markers)
                  return(bi_K(data = data_new, mark_pair = grid[.x,1:2], r = r, 
                              correction = 'trans', id = id, iter = grid[.x,3], win = win))
                }
  )
  colnames(perm)[c(2,6,7)] = c(id, 'Theoretical CSR', 'Permuted K')
  
  obs = purrr::map_df(.x = 1:nrow(grid[!duplicated(grid[,1:2]),]), 
               ~bi_K(data = data, mark_pair = grid[.x,1:2], r = r, 
                     correction = 'trans', id = id, iter = 1, win = win)) %>%
    dplyr::select(-iter)
  colnames(obs)[c(1,5,6)] = c(id, 'Theoretical CSR', 'Observed K')
  
  final = suppressMessages(dplyr::left_join(perm, obs)) 
  
  if(method == 'M'){
    final = final %>% dplyr::mutate_at(c("Theoretical CSR","Permuted K","Observed K"),
                                ~./(pi*r^2)) %>%
      dplyr::mutate('Degree of Clustering Permutation' = `Permuted K`,
             'Degree of Clustering Theoretical' = `Theoretical CSR`) %>%
      dplyr::rename('Permuted M' = `Permuted K`,
             'Theoretical CSR' = `Theoretical CSR`,
             'Observed M' = `Observed K`)
  }
  
  #Conduct the necessary transformation
  if(method == 'L'){
    final = final %>% dplyr::mutate_at(c("Theoretical CSR","Permuted K","Observed K"),
                                ~sqrt(./pi)) %>%
      dplyr::mutate('Degree of Clustering Permutation' = `Observed K` - `Permuted K`,
             'Degree of Clustering Theoretical' = `Observed K` - `Theoretical CSR`)  %>%
      dplyr::rename('Permuted L' = `Permuted K`,
             'Theoretical CSR' = `Theoretical CSR`,
             'Observed L' = `Observed K`)
  }
  
  if(method == 'K'){
    final = final %>%
      dplyr::mutate('Degree of Clustering Permutation' = `Observed K` - `Permuted K`,
             'Degree of Clustering Theoretical' = `Observed K` - `Theoretical CSR`)
  }
  
  if(!perm_dist){
    final = final %>% 
      dplyr::mutate(id = get(id)) %>%
      dplyr::select(-c(1,2)) %>%
      dplyr::group_by(id, anchor, counted, r) %>%
      #dplyr::summarize(`Theoretical CSR` = mean(`Theoretical CSR`),
      #          `Permuted CSR` = mean(.[[grep('Permuted', colnames(.), value = TRUE)]],
      #                                na.rm = TRUE),
      #          `Observed` = mean(.[[grep('Observed', colnames(.), value = TRUE)]],
      #                            na.rm = TRUE),
      #          `Degree of Clustering Theoretical` = mean(`Degree of Clustering Theoretical`,
      #                                                    na.rm = TRUE),
      #          `Degree of Clustering Permutation` =  mean(`Degree of Clustering Permutation`,
      #          na.rm = TRUE))
      dplyr::summarize_all(~mean(.,na.rm = TRUE))
  }
  colnames(final)[1] = id
  return(final)
}

## Nearest Neighbor

G_out = function(data, marker, id, iter, correction,r_value, win){
  #Does the actual computation of Ripley' K
  G_obs = spatstat.geom::ppp(x = data$xloc, y = data$yloc, window = win) %>%
    spatstat.core::Gest(r = r_value, correction = correction) %>%
    data.frame() %>%
    dplyr::filter(r != 0) %>%
    dplyr::mutate(Marker = marker,
           label = data[[id]][1],
           iter = iter) %>%
    dplyr::select(iter, label, Marker, r, theo, correction) %>%
    dplyr::mutate(label = as.character(label))
  return(G_obs)
}


uni_G = function(data = data, iter, marker, id, correction, r_value, win){
  #Filters down to the positive cells and then computes  G
  data_pos = data %>% dplyr::filter(Marker == marker,  Positive == 1)
  
  if(nrow(data_pos)==0){
    G = data.frame(iter = iter,label = data[[id]][1],
                   Marker = marker, r = r_value, 
                   theo = NA, correction = NA) %>%
      dplyr::filter(r > 0) %>%
      dplyr::mutate(label = as.character(label))
    colnames(G)[6] = correction
  }else{
    G = G_out(data = data_pos, marker = marker, 
              id = id, iter = iter, correction = correction,
              r_value = r_value, win = win)
  }
  return(G)
}


uni_NN_G = function(data, markers, id, num_iters, correction,
                    perm_dist, r){
  #Main function
  
  #Notice that this follows spatstat's notation and argument name
  if(!(correction %in% c("rs", "km", "han"))){
    stop("Did not provide a valid edge correcion method.")
  }
  
  #Use set the cell location as the center of the cell
  data = data %>% 
    dplyr::mutate(xloc = (XMin + XMax)/2,
           yloc = (YMin + YMax)/2
    )
  
  #Create the region that the point process exists. This only needs to be done
  #once per image
  win = spatstat.geom::convexhull.xy(x = data$xloc, y = data$yloc)
  
  #Make the data a long format
  data = data %>%
    tidyr::pivot_longer(cols = all_of(markers), names_to = 'Marker', values_to = 'Positive')
  
  grid = expand.grid(markers, 1:num_iters) %>%
    dplyr::mutate(Var1 = as.character(Var1))
  
  perms = purrr::map_df(.x = 1:nrow(grid),~{
    data_new = perm_data(data, markers)
    return(uni_G(data = data_new, iter = grid[.x,2], 
                 marker = grid[.x,1], id = id, r_value = r,
                 correction = correction, win = win))})
  colnames(perms)[c(2,5,6)] = c(id, 'Theoretical CSR','Permuted CSR')
  
  obs = purrr::map_df(.x = markers, ~uni_G(data = data, iter = 'Observed', 
                                    marker = .x, correction = correction,
                                    id = id, r_value = r,
                                    win = win)) %>%
   dplyr::select(-iter) 
  
  colnames(obs)[c(1,4,5)] = c(id, 'Theoretical CSR', 'Observed')
  
  final = suppressMessages(dplyr::left_join(perms, obs)) %>%
    dplyr::mutate(
      `Degree of Clustering Theoretical` = (`Observed`) - (`Theoretical CSR`),
      `Degree of Clustering Permutation` = (`Observed`) - (`Permuted CSR`)) %>%
    dplyr::select(-iter)
  
  if(!perm_dist){
    final = final %>% 
      dplyr::mutate(id = .data[[id]]) %>%
      dplyr::select(-1) %>%
      dplyr::group_by(id, Marker, r) %>%
      #dplyr::summarize(`Theoretical CSR` = mean(`Theoretical CSR`,na.rm = TRUE),
      #                 `Permuted CSR` = mean(.[[grep('Permuted', colnames(.), value = TRUE)]],
      #                                       na.rm = TRUE),
      #                 `Observed` = mean(.[[grep('Observed', colnames(.), value = TRUE)]],
      #                                   na.rm = TRUE),
      #                 `Degree of Clustering Theoretical` = mean(`Degree of Clustering Theoretical`,
      #                                                           na.rm = TRUE),
      #                 `Degree of Clustering Permutation` =  mean(`Degree of Clustering Permutation`,
      #                                                            na.rm = TRUE))
      dplyr::summarize_all(~mean(.,na.rm = TRUE))
    colnames(final)[1] = id
  }
  
  
  return(final) 
}

bi_G = function(data, mark_pair, r, correction, id, iter, win){
  mark_pair = unname(mark_pair)
  #Keeps all positive cell for either marker, and then discards the repeated locations
  data_new = data %>% 
    dplyr::filter(Marker %in% c(mark_pair[1], mark_pair[2]),
           Positive == 1) %>%
    dplyr::group_by(xloc, yloc) %>%
    dplyr::filter(dplyr::n() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Marker = factor(Marker, levels = c(mark_pair[1], mark_pair[2])))
  
  
  if(any(table(data_new$Marker) == 0)){
    #weird error here :row names were found from a short variable and have been discarded
    G = data.frame(r = r, theo = NA, trans = NA, 
                   anchor = levels(data_new$Marker)[1], 
                   counted = levels(data_new$Marker)[2],
                   label = data[[id]][1], 
                   iter = iter,
                   theo = NA,
                   correction = NA)  %>% 
      dplyr::select(iter, label, anchor, counted,  r, theo, correction) %>% 
      dplyr::filter(r>0)
    colnames(G)[7] = correction
  }else{
    X = spatstat.geom::ppp(x = data_new$xloc, y = data_new$yloc, window = win,
                           marks = data_new$Marker)
    G = spatstat.core::Gcross(r = r, X = X, i = levels(data_new$Marker)[1], 
               j = levels(data_new$Marker)[2], correction = correction) %>% 
      data.frame() %>%
      dplyr::filter(r!=0) %>%
      dplyr::mutate(anchor = levels(data_new$Marker)[1],
             counted = levels(data_new$Marker)[2],
             iter = iter,
             label = data[[id]][1]) %>% 
      dplyr::select(iter, label, anchor, counted, r, theo, correction) 
  }
  return(G)
}

bi_NN_G_sample = function(data, markers, id, num_iters, correction, 
                   perm_dist, r, exhaustive){
  #Main function
  #Notice that this follows spatstat's notation and argument name
  if(!(correction %in% c('rs', 'hans'))){
    stop("Did not provide a valid edge correcion method.")
  }
  
  #Check if exhaustive is FALSE, then markers must be a data.frame
  if(exhaustive == FALSE & class(markers) != 'data.frame'){
    stop("If exhaustive == FALSE, then markers must be a data.frame.")
  }
  if(exhaustive == TRUE & class(markers) != 'character'){
    stop("If exhaustive == TRUE, then markers must be a character vector")
  }
  
  
  #Use set the cell location as the center of the cell
  data = data %>% 
    dplyr::mutate(xloc = (XMin + XMax)/2,
           yloc = (YMin + YMax)/2
    )
  
  #Create the region that the point process exists. This only needs to be done
  #once per image
  win = spatstat.geom::convexhull.xy(x = data$xloc, y = data$yloc)

  #Enumerates all possible combination of markers, and removes the ones where
  #marker 1 and marker 2 are the same
  if(exhaustive){
    data = data %>%
      tidyr::pivot_longer(cols = all_of(markers), 
                          names_to = 'Marker', values_to = 'Positive')
    grid = expand.grid(markers, markers, 1:num_iters) %>%
      dplyr::mutate(Var1 = as.character(Var1),
                    Var2 = as.character(Var2)) %>%
      dplyr::filter(Var1 != Var2)
  }else{
    grid = expand_grid(markers, 1:num_iters) %>%
      data.frame() %>%
      dplyr::mutate(Var1 = as.character(Var1),
                    Var2 = as.character(Var2)) %>%
      dplyr::filter(Var1 != Var2)
    markers = grid %>% dplyr::select(Var1, Var2) %>%
      unlist() %>% unique()
    data = data %>%
      tidyr::pivot_longer(cols = all_of(markers), 
                          names_to = 'Marker', values_to = 'Positive')
  }
  
  
  perm = purrr::map_df(.x = 1:nrow(grid), 
                ~{
                  data_new = perm_data(data, markers)
                  return(bi_G(data = data_new, mark_pair = grid[.x,1:2], r = r, 
                              correction = correction, id = id, iter = grid[.x,3], win = win))
                }
  )
  colnames(perm)[c(2,6,7)] = c(id, 'Theoretical CSR', 'Permuted G')
  
  obs = purrr::map_df(.x = 1:nrow(grid[!duplicated(grid[,1:2]),]), 
               ~bi_G(data = data, mark_pair = grid[.x,1:2], r = r, 
                     correction = correction, id = id, iter = 1, win = win)) %>%
    dplyr::select(-iter)
  colnames(obs)[c(1,5,6)] = c(id, 'Theoretical CSR', 'Observed G')
  
  final = suppressMessages(dplyr::left_join(perm, obs)) %>%
    dplyr::mutate(`Degree of Clustering Permutation` = ifelse(`Permuted G` == 0, NA, (`Observed G`)-(`Permuted G`)),
           `Degree of Clustering Theoretical` = ifelse(`Theoretical CSR` == 0, NA, (`Observed G`)-(`Theoretical CSR`)))
  
  if(!perm_dist){
    final = final %>% 
      dplyr::mutate(id = .data[[id]]) %>%
      dplyr::select(-c(1,2)) %>%
      dplyr::group_by(id, anchor, counted, r) %>%
      #dplyr::summarize(`Theoretical CSR` = mean(`Theoretical CSR`,na.rm = TRUE),
      #          `Permuted CSR` = mean(.[[grep('Permuted', colnames(.), value = TRUE)]],
      #                                na.rm = TRUE),
      #          `Observed` = mean(.[[grep('Observed', colnames(.), value = TRUE)]],
      #                            na.rm = TRUE),
      #          `Degree of Clustering Theoretical` = mean(`Degree of Clustering Theoretical`,
      #                                                    na.rm = TRUE),
      #          `Degree of Clustering Permutation` =  mean(`Degree of Clustering Permutation`,
      #                                                     na.rm = TRUE))
      dplyr::summarize_all(~mean(.,na.rm = TRUE))

    colnames(final)[1] = id
  }
  
  
  return(final)
}

