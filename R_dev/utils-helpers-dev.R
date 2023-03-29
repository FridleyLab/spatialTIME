#utils-helpers dev

#calculate dixons s Z table
dix_s_z = function(data, markers, num_permutations, xloc, yloc){
  #in function that calls dix_s_classifier needs to assign or find x and y locations
  #see uni_rip_k function 
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
  #identify the different classifier levels
  #loop through for markers
  dixon_s = purrr::map(.x = 1:nrow(markers), 
                       ~{
                         marker = markers[.x,] %>% unlist() %>% as.character()
                         #print(marker)
                         tmp = data %>% 
                           dplyr::filter(get(!!marker[1]) == 1 | get(!!marker[2]) == 1, get(!!marker[1]) != get(!!marker[2])) %>% 
                           dplyr::select(xloc, yloc, !!marker) %>% 
                           tidyr::gather("Marker", "Positive", -xloc, -yloc) %>% 
                           dplyr::filter(Positive == 1)
                         #check if marker is only in one tissue type or there are zero rows
                         if(length(unique(tmp[["Marker"]])) == 1){
                           dat = expand.grid("From" = marker, "To" = marker) %>%
                             dplyr::bind_cols(data.frame("    Obs.Count" = rep(NA, nrow(.)),
                                                         "    Exp. Count" = rep(NA, nrow(.)),
                                                         "S " = rep(NA, nrow(.)),
                                                         "Z " = rep(NA, nrow(.)),
                                                         "  p-val.Z" = rep(NA, nrow(.)),
                                                         "  p-val.Nobs" = rep(NA, nrow(.)), check.names = F))
                           ns = tmp %>% dplyr::group_by(Marker) %>% dplyr::summarise(dplyr::n()) %>% dplyr::pull(2, 1)
                           ns_mat = matrix(data = rep(ns, 4), nrow = 4, byrow = T) %>% data.frame(check.names=F)
                           colnames(ns_mat) = names(ns)
                           dat = dplyr::bind_cols(dat, ns_mat)
                           #if there are zero can skip finding how many there are
                         } else if(nrow(tmp) == 0){
                           dat = expand.grid("From" = marker, "To" = marker) %>%
                             dplyr::bind_cols(data.frame("    Obs.Count" = rep(NA, nrow(.)),
                                                         "    Exp. Count" = rep(NA, nrow(.)),
                                                         "S " = rep(NA, nrow(.)),
                                                         "Z " = rep(NA, nrow(.)),
                                                         "  p-val.Z" = rep(NA, nrow(.)),
                                                         "  p-val.Nobs" = rep(NA, nrow(.)), check.names = F))
                           ns = c(0, 0)
                           ns_mat = matrix(data = rep(ns, 4), nrow = 4, byrow = T) %>% data.frame(check.names=F)
                           names(ns_mat) = marker
                           dat = dplyr::bind_cols(dat, ns_mat)
                           #if there less than 2 in all classifiers
                         } else if(TRUE %in% ((tmp %>% dplyr::group_by(Marker) %>% dplyr::summarise(dplyr::n()) %>%
                                               dplyr::pull(2, 1)) < 3)){
                           dat = expand.grid("From" = marker, "To" = marker) %>%
                             dplyr::bind_cols(data.frame("    Obs.Count" = rep(NA, nrow(.)),
                                                         "    Exp. Count" = rep(NA, nrow(.)),
                                                         "S " = rep(NA, nrow(.)),
                                                         "Z " = rep(NA, nrow(.)),
                                                         "  p-val.Z" = rep(NA, nrow(.)),
                                                         "  p-val.Nobs" = rep(NA, nrow(.)), check.names = F))
                           ns = tmp %>% dplyr::group_by(Marker) %>% dplyr::summarise(dplyr::n()) %>% dplyr::pull(2, 1)
                           ns_mat = matrix(data = rep(ns, 4), nrow = 4, byrow = T) %>% data.frame(check.names=F)
                           colnames(ns_mat) = names(ns)
                           dat = dplyr::bind_cols(dat, ns_mat)
                           #if dixons s can be computed do so
                         } else {
                           invisible(capture.output(dat <- dixon::dixon(tmp, nsim = 1)$tablaZ))
                           ns = tmp %>% dplyr::group_by(Marker) %>% dplyr::summarise(dplyr::n()) %>% dplyr::pull(2, 1)
                           ns_mat = matrix(data = rep(ns, 4), nrow = 4, byrow = T) %>% data.frame(check.names=F)
                           colnames(ns_mat) = names(ns)
                           dat = dplyr::bind_cols(dat, ns_mat)
                         }
                         return(dat)
                       }
  ) %>%
    #bind the resulting dfs together
    do.call(dplyr::bind_rows, .)
  
  marker_combination = dixon_s %>% dplyr::mutate(grouping = rep(1:(nrow(.)/4),  each = 4)) %>%
    dplyr::group_by(grouping) %>% 
    dplyr::arrange(From, To, .by_group = T) %>%
    dplyr::select(1:2, grouping) %>% 
    dplyr::slice(2)
  Obs_table = dixon_s %>% 
    dplyr::mutate(grouping = rep(1:(nrow(.)/4), each = 4)) %>%
    dplyr::group_by(grouping) %>% 
    dplyr::arrange(From, To, .by_group = T) %>% 
    dplyr::select(grouping, `    Obs.Count`) %>%
    dplyr::mutate(Key = c("Naa", "Nab", "Nba", "Nbb")) %>%
    dplyr::mutate(`    Obs.Count` = ifelse(is.na(`    Obs.Count`), 0, `    Obs.Count`)) %>%
    tidyr::spread("Key", "    Obs.Count") %>%
    dplyr::mutate(Na = Naa + Nab, Nb = Nba + Nbb,
                  S_A_star = log(((Naa+1)/(Nab+1))/((Na+1)/(Nb+2))),
                  S_B_star = log(((Nbb+1)/(Nba+1))/((Nb+1)/(Na+2))))
  final_table = dplyr::full_join(marker_combination, Obs_table, by = "grouping")
  return("Adjusted Dixon" = dixon_s)
}

#calculate dixons s C table
dix_s_c = function(data, markers, classifier_label, num_permutations, xloc, yloc){
  #in function that calls dix_s_classifier needs to assign or find x and y locations
  #see uni_rip_k function 
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
  #identify the different classifier levels
  
  #loop through for markers
  dixon_s = purrr::map(.x = 1:nrow(markers), 
                       ~{
                         marker = markers[.x,] %>% unlist() %>% as.character()
                         #print(marker)
                         tmp = data %>% 
                           dplyr::filter(get(!!marker[1]) == 1 | get(!!marker[2]) == 1, get(!!marker[1]) != get(!!marker[2])) %>% 
                           dplyr::select(xloc, yloc, !!marker) %>% 
                           tidyr::gather("Marker", "Positive", -xloc, -yloc) %>% 
                           dplyr::filter(Positive == 1)
                         #check if marker is only in one tissue type or there are zero rows
                         if(length(unique(tmp[["Marker"]])) == 1){
                           dat = data.frame("  df " = rep(NA, length(marker) + 1),
                                            "Chi-sq" = rep(NA, length(marker) + 1),
                                            "P.asymp" = rep(NA, length(marker) + 1),
                                            "  P.rand" = rep(NA, length(marker) + 1), check.names = F) %>%
                             mutate(Segregation = c("Overall segregation", paste("From", marker)), .before = 1)
                           ns = tmp %>% dplyr::group_by(Marker) %>% dplyr::summarise(dplyr::n()) %>% pull(2, 1)
                           ns = c(ns, 0)
                           names(ns)[2] = marker[!(marker %in% names(ns))]
                           ns_mat = matrix(data = rep(ns, length(marker) + 1), nrow = length(marker) + 1, byrow = T) %>% data.frame(check.names=F)
                           colnames(ns_mat) = names(ns)
                           dat = dplyr::bind_cols(dat, ns_mat)
                           #if there are zero can skip finding how many there are
                         } else if(nrow(tmp) == 0){
                           dat = data.frame("  df " = rep(NA, length(marker) + 1),
                                            "Chi-sq" = rep(NA, length(marker) + 1),
                                            "P.asymp" = rep(NA, length(marker) + 1),
                                            "  P.rand" = rep(NA, length(marker) + 1), check.names = F) %>%
                             mutate(Segregation = c("Overall segregation", paste("From", marker)), .before = 1)
                           ns = c(0, 0)
                           ns_mat = matrix(rep(ns, length(marker) + 1), nrow = length(marker) + 1, byrow = T) %>% data.frame(check.names = F)
                           colnames(ns_mat) = marker
                           dat = dplyr::bind_cols(dat, ns_mat)
                           #if there less than 2 in all classifiers
                         } else if(TRUE %in% ((tmp %>% dplyr::group_by(Marker) %>% dplyr::summarise(dplyr::n()) %>%
                                               dplyr::pull(2, 1)) < 2)){
                           dat = data.frame("  df " = rep(NA, length(marker) + 1),
                                            "Chi-sq" = rep(NA, length(marker) + 1),
                                            "P.asymp" = rep(NA, length(marker) + 1),
                                            "  P.rand" = rep(NA, length(marker) + 1), check.names = F) %>%
                             mutate(Segregation = c("Overall segregation", paste("From", marker)), .before = 1)
                           ns = tmp %>% dplyr::group_by(Marker) %>% dplyr::summarise(dplyr::n()) %>% dplyr::pull(2, 1)
                           ns_mat = matrix(data = rep(ns, length(marker) + 1), nrow = length(marker) + 1, byrow = T) %>% data.frame(check.names=F)
                           colnames(ns_mat) = names(ns)
                           dat = dplyr::bind_cols(dat, ns_mat)
                           #if dixons s can be computed do so
                         } else {
                           invisible(capture.output(dat <- dixon::dixon(tmp, nsim = num_permutations)$tablaC))
                           dat = dat %>%
                             tibble::rownames_to_column("Segregation") %>%
                             dplyr::mutate(Segregation = Segregation %>% gsub("\\  ", " ", .) %>% gsub("\\  *$", "", .))
                           ns = tmp %>% dplyr::group_by(Marker) %>% dplyr::summarise(dplyr::n()) %>% dplyr::pull(2, 1)
                           ns_mat = matrix(data = rep(ns, length(marker) + 1), nrow = length(marker) + 1, byrow = T) %>% data.frame(check.names=F)
                           colnames(ns_mat) = names(ns)
                           dat = dplyr::bind_cols(dat, ns_mat)
                         }
                         return(dat)
                       }
  ) %>%
    #bind the resulting dfs together
    do.call(dplyr::bind_rows, .)
  return(dixon_s)
}
