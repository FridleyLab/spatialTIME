K_out = function(data, marker, id, iter, correction,r_value, win){
  un = NULL
  #Does the actual computation of Ripley' K
  K_obs = spatstat.geom::ppp(x = data$xloc, y = data$yloc, window = win) %>%
    spatstat.explore::Kest(r = r_value, correction = correction) %>%
    data.frame() %>%
    dplyr::filter(r != 0) %>%
    dplyr::mutate(Marker = marker,
           label = data[[id]][1],
           iter = iter)
  if(correction=="none"){
    K_obs = K_obs %>%
      dplyr::select(iter, label, Marker, r, theo, un) %>%
      dplyr::mutate(label = as.character(label))
    return(K_obs)
  } else {
    K_obs = K_obs %>%
      dplyr::select(iter, label, Marker, r, theo, correction) %>%
      dplyr::mutate(label = as.character(label))
    return(K_obs)
  }
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
                     perm_dist, r,xloc, yloc){
  #Main function
  #Check if the method is selected from K, L, M
  if(!(method %in% c('K',"L","M"))){
    stop("Did not provide a valid method.")
  }
  
  #Notice that this follows spatstat's notation and argument name
  if(!(correction %in% c('trans', "iso", "border", 'translation', 'isotropic', 'none'))){
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
  
  
  if(is.null(xloc) + is.null(yloc) == 1){
    stop("Both xloc and yloc must be either NULL or a defined column")
  }
  
  if(!is.null(xloc) && !is.null(yloc)){
    data = data %>%
      dplyr::mutate(XMin = get(xloc),
             XMax = get(xloc),
             YMin = get(yloc),
             YMax = get(yloc))
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
    tidyr::pivot_longer(cols = tidyselect::all_of(markers), names_to = 'Marker', values_to = 'Positive')
  
  grid = expand.grid(markers, 1:num_iters) %>%
    dplyr::mutate(Var1 = as.character(Var1))
  
  perms = purrr::map_df(.x = 1:nrow(grid),~{
    data_new = perm_data(data, markers)
    print(which(data_new$CD3..CD8. == 1))
    return(uni_K(data = data_new, iter = grid[.x,2], 
                 marker = grid[.x,1], id = id, r_value = r,
                 correction = correction, win = win))})
  if(correction == "none"){
    perms = perms %>%
      select(-!!correction) %>%
      rename(!!correction := un)
  }
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
      dplyr::group_by(id, Marker, r) %>%
      #dplyr::summarize(`Theoretical CSR` = mean(`Theoretical CSR`, na.rm = TRUE),
      #          `Permuted CSR` = mean(.[[grep('Permuted', colnames(.), value = TRUE)]], na.rm = TRUE),
      #          `Observed` = mean(.[[grep('Observed', colnames(.), value = TRUE)]], na.rm = TRUE),
      #          `Degree of Clustering Theoretical` = mean(`Degree of Clustering Theoretical`, na.rm = TRUE),
      #          `Degree of Clustering Permutation` =  mean(`Degree of Clustering Permutation`, na.rm = TRUE)
      #          )
      dplyr::summarize_all(~mean(., na.rm = TRUE))
      colnames(final)[1] = id
  }else{
    colnames(final)[2] = id
  }
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
    K = spatstat.explore::Kcross(r =r, X = X, i = levels(data_new$Marker)[1], 
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
                    method, perm_dist, r, exhaustive, xloc, yloc){
  #Main function
  #Check if the method is selected from K, L, M
  if(!(method %in% c('K',"L","M"))){
    stop("Did not provide a valid method.")
  }
  
  #Check if exhaustive is FALSE, then markers must be a data.frame
  if(exhaustive == FALSE & !is(markers, 'data.frame')){
    stop("If exhaustive == FALSE, then markers must be a data.frame.")
  }
  if(exhaustive == TRUE & !is(markers, 'character')){
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
  
  
  if(is.null(xloc) + is.null(yloc) == 1){
    stop("Both xloc and yloc must be either NULL or a defined column")
  }
  
  if(!is.null(xloc) && !is.null(yloc)){
    data = data %>%
      dplyr::mutate(XMin = get(xloc),
                    XMax = get(xloc),
                    YMin = get(yloc),
                    YMax = get(yloc))
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
      tidyr::pivot_longer(cols = tidyselect::all_of(markers), 
                          names_to = 'Marker', values_to = 'Positive')
    grid = expand.grid(markers, markers, 1:num_iters) %>%
      dplyr::mutate(Var1 = as.character(Var1),
           Var2 = as.character(Var2)) %>%
      dplyr::filter(Var1 != Var2)
  }else{
    grid = tidyr::expand_grid(markers, 1:num_iters) %>%
      data.frame() %>%
      dplyr::mutate(Var1 = as.character(Var1),
             Var2 = as.character(Var2)) %>%
      dplyr::filter(Var1 != Var2)
    markers = grid %>% dplyr::select(Var1, Var2) %>%
      unlist() %>% unique()
    data = data %>%
      tidyr::pivot_longer(cols = tidyselect::all_of(markers), 
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
  colnames(final)[2] = id
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
  }else{
    colnames(final)[2] = id
  }
  return(final)
}

## Nearest Neighbor

G_out = function(data, marker, id, iter, correction,r_value, win){
  #Does the actual computation of Ripley' K
  G_obs = spatstat.geom::ppp(x = data$xloc, y = data$yloc, window = win) %>%
    spatstat.explore::Gest(r = r_value, correction = correction) %>%
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
                    perm_dist, r, xloc, yloc){
  #Main function
  
  #Notice that this follows spatstat's notation and argument name
  if(!(correction %in% c("rs", "km", "han"))){
    stop("Did not provide a valid edge correcion method.")
  }
  
  if(is.null(xloc) + is.null(yloc) == 1){
    stop("Both xloc and yloc must be either NULL or a defined column")
  }
  
  if(!is.null(xloc) && !is.null(yloc)){
    data = data %>%
      dplyr::mutate(XMin = get(xloc),
                    XMax = get(xloc),
                    YMin = get(yloc),
                    YMax = get(yloc))
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
    tidyr::pivot_longer(cols = tidyselect::all_of(markers), names_to = 'Marker', values_to = 'Positive')
  
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
      `Degree of Clustering Permutation` = (`Observed`) - (`Permuted CSR`))
  
  if(!perm_dist){
    final = final %>% 
      dplyr::mutate(id = .data[[id]]) %>%
      dplyr::select(-(1:2)) %>%
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
  }else{
    colnames(final)[2] = id
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
    G = spatstat.explore::Gcross(r = r, X = X, i = levels(data_new$Marker)[1], 
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
                   perm_dist, r, exhaustive, xloc, yloc){
  #Main function
  #Notice that this follows spatstat's notation and argument name
  if(!(correction %in% c('rs', 'hans'))){
    stop("Did not provide a valid edge correcion method.")
  }
  
  #Check if exhaustive is FALSE, then markers must be a data.frame
  if(exhaustive == FALSE & !is(markers, 'data.frame')){
    stop("If exhaustive == FALSE, then markers must be a data.frame.")
  }
  if(exhaustive == TRUE & !is(markers, 'character')){
    stop("If exhaustive == TRUE, then markers must be a character vector")
  }

  
  if(is.null(xloc) + is.null(yloc) == 1){
    stop("Both xloc and yloc must be either NULL or a defined column")
  }
  
  if(!is.null(xloc) && !is.null(yloc)){
    data = data %>%
      dplyr::mutate(XMin = get(xloc),
                    XMax = get(xloc),
                    YMin = get(yloc),
                    YMax = get(yloc))
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
      tidyr::pivot_longer(cols = tidyselect::all_of(markers), 
                          names_to = 'Marker', values_to = 'Positive')
    grid = expand.grid(markers, markers, 1:num_iters) %>%
      dplyr::mutate(Var1 = as.character(Var1),
                    Var2 = as.character(Var2)) %>%
      dplyr::filter(Var1 != Var2)
  }else{
    grid = tidyr::expand_grid(markers, 1:num_iters) %>%
      data.frame() %>%
      dplyr::mutate(Var1 = as.character(Var1),
                    Var2 = as.character(Var2)) %>%
      dplyr::filter(Var1 != Var2)
    markers = grid %>% dplyr::select(Var1, Var2) %>%
      unlist() %>% unique()
    data = data %>%
      tidyr::pivot_longer(cols = tidyselect::all_of(markers), 
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
  colnames(final)[2] = id
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
  }else{
    colnames(final)[2] = id
  }
  
  
  return(final)
}

list.append = function(list, new){
  new_list = list
  new_list[[length(new_list) + 1]] = new
  return(new_list)
}


#special thanks to Simon Vandekar and Julia Wrobel
get_exactK = function(pp_obj,
                     mark1,
                     mark2 = NULL,
                     r_vec = NULL,
                     ...){
  
  n = pp_obj$n
  K = spatstat.explore::Kest(pp_obj, r = r_vec, ...)
  sumW = K * n * (n-1)
  
  if(is.null(mark2)){
    X1 = ifelse(pp_obj$marks == mark1, 1, 0)
    sX1 = sum(X1)
    #v1 = sX1/n
    v2 = (sX1 * (sX1 - 1))/(n * (n-1))
    lambda2 = (sX1 *(sX1-1))
  }else{
    X1 = ifelse(pp_obj$marks == mark1, 1, 0)
    X2 = ifelse(pp_obj$marks == mark2, 1, 0)
    sX1 = sum(X1)
    sX2 = sum(X2)
    sX1X2 = sum(X1*X2)
    
    v1 = sX1X2/n
    v2 = (sX1*sX2-sX1X2)/(n*(n-1))
    lambda2 = sX1*sX2/(areapp^2)
  }
  
  Kperm = v2 * sumW/lambda2
  tidyr::as_tibble(Kperm[,-2])
}

#get tile counts
getTile = function(slide, l, size){
  if(slide == 1){
    v = rep(TRUE, l)
    return(list(v))
  }
  lapply(1:slide, function(s){
    if(s == 1){
      w = 1:size
    }
    if(s != 1 & s != slide){
      w = (((s-1)*size)+1):((s)*size)
    }
    if(s == slide){
      w = (((s-1)*size)+1):(l)
    }
    v = rep(FALSE, l)
    v[w] = TRUE
    v
  })
}

calculateK = function(i_dat, j_dat, anchor, counted, area, win, big, r_range, edge_correction, cores){
  li = nrow(i_dat)
  lj = nrow(j_dat)
  #find intensity of each marker
  lambdai = li/area
  lambdaj = lj/area
  #point pattern for anchor for distance matrix
  ppi = spatstat.geom::ppp(x = i_dat[,1],
                           y = i_dat[,2],
                           window = win)
  #point pattern for counted for distance matrix
  ppj = spatstat.geom::ppp(x = j_dat[,1],
                           y = j_dat[,2],
                           window = win)
  if(li > big | lj > big){
    #find the number of tiles of distance matrix in columns and rows
    i_slides = ceiling(li / big)
    j_slides = ceiling(lj / big)
    #produce vectors of length ns, create T/F vector of which to subset
    i_ranges = getTile(slide = i_slides, l = li, size = big)
    j_ranges = getTile(slide = j_slides, l = lj, size = big)
    #count up those within the specified distance
    counts = parallel::mclapply(i_ranges, function(i_section){
      j_out = parallel::mclapply(j_ranges, function(j_section){
        #subset to tile
        i_tmp = ppi[i_section]
        j_tmp = ppj[j_section]
        
        #if the correction method is set to border, then run the num and den from spatstat Kmulti
        if(edge_correction %in% c("border")){
          bI = spatstat.geom::bdist.points(i_tmp)
          bcloseI = bI[spatstat.geom::crosspairs(i_tmp, j_tmp, max(r_range), what = "ijd")$i]
          RS = spatstat.explore::Kount(spatstat.geom::crosspairs(i_tmp, j_tmp, max(r_range), what = "ijd")$d,
                                    spatstat.geom::crosspairs(i_tmp, j_tmp, max(r_range), what = "ijd")$i,
                                    bI, spatstat.geom::handle.r.b.args(r_range, breaks=NULL, win, rmaxdefault = max(r_range)))
          return(RS)
        }
        #calculate distances
        dists = spatstat.geom::crossdist(i_tmp, j_tmp)
        rmv_i = rowSums(dists < max(r_range)) != 0
        rmv_j = colSums(dists < max(r_range)) != 0
        i_tmp = i_tmp[rmv_i] #I
        j_tmp = j_tmp[rmv_j] #J
        dists = spatstat.geom::crossdist(i_tmp, j_tmp)
        if(0 %in% dim(dists)){
          return(rep(0, length(r_range)))
        }
        #calculate edge correcion
        if(edge_correction %in% c("trans", "translation")){
          edge = spatstat.explore::edge.Trans(i_tmp, j_tmp)
          #count edge correction matrix for cells within range r in distance matrix
          counts = sapply(r_range, function(r){sum(edge[which(dists < r)])})
          #remove large distance and edge correction matrix to keep ram usage down
          rm(dists, edge)
          #return counts for tile
          return(counts)
        }
        if(edge_correction %in% c("none")){
          counts = cumsum(spatstat.geom::whist(spatstat.geom::crosspairs(i_tmp, j_tmp, max(r_range), what = "ijd")$d,
                                               spatstat.geom::handle.r.b.args(r_range, breaks=NULL, win, rmaxdefault = max(r_range))$val))
          return(counts)
        }
      }) 
      
      if(edge_correction == "border"){
        num = lapply(j_out, function(j_big){
          j_big[[1]]
        }) %>%
          do.call(rbind.data.frame, .) %>%
          colSums() %>%
          unname()
        den = j_out[[1]][[2]]
        return(list(num = num, den = den))
      }
      j_out  %>% #use 1 core per tile
        #bind all j tiles to data frame
        do.call(rbind.data.frame, .) %>%
        #take column sums and return
        colSums()
    },mc.cores = cores, mc.preschedule = FALSE, mc.allow.recursive = TRUE)
    
    if(edge_correction == "border"){
      num = lapply(counts, function(j_big){
        j_big[[1]]
      }) %>%
        do.call(rbind.data.frame, .) %>%
        colSums() %>%
        unname()
      den = ppi$n
      estimatedK = num / (lambdaj * den)
    } else {
      counts = counts %>%
        #bind i tile counts and sum
        do.call(rbind.data.frame, .) %>%
        colSums() %>%
        unname()
      #calulate clustering from counted edge corrections with intensities of i and j
      estimatedK = (1/(lambdai * lambdaj * area)) * counts
    }
  }
  #if there are less than 10,000 j or i just compute matrix
  if(!(li > big | lj > big)){
    #calculate distance matrix
    sp_tmp2 = as.data.frame(rbind(i_dat, j_dat)) %>% 
      dplyr::mutate(Marker = c(rep(anchor, nrow(i_dat)), rep(counted, nrow(j_dat))))
    pp_obj = spatstat.geom::ppp(x = sp_tmp2$xloc, y = sp_tmp2$yloc, window = win, marks = factor(sp_tmp2$Marker))
    estimatedK = spatstat.explore::Kcross(pp_obj,
                                          i = unique(sp_tmp2$Marker)[1],
                                          j = unique(sp_tmp2$Marker)[2],
                                          r = r_range,
                                          correction = edge_correction) %>%
      data.frame() %>%
      dplyr::pull(3)
  }
  
  return(estimatedK)
}
