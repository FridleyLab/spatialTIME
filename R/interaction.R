#' Bivariate Interaction Variable
#' 
#' Single-cell spatial-protein metric introduce by Steinhart et al in https://doi.org/10.1158/1541-7786.mcr-21-0411
#'
#' @param mif object of class `mif`
#' @param mnames a character vector or table with 2 columns indicating the from-to markers to assess
#' @param r_range numeric vector of radii for which to calculate the interaction variable at
#' @param num_permutations integer for how many permutations to use to derive the interaction estimate under CSR
#' @param keep_permutation_distribution boolean for whether or not to keep all permutation results or average them
#' @param workers integer for the number of CPU cores to use for permutations, markers, and spatial samples
#' @param overwrite boolean for whether to overwrite existing interaction variable results
#' @param xloc column name in spatial files containing the x location - if left NULL will average columns XMin and XMax
#' @param yloc column name in spatial files containing the y location - if left NULL will average columns YMin and YMax
#'
#' @return object of class mif with the interaction variable derive slot filled
#' @export
#'
interaction_variable = function(mif, 
                                mnames,
                                r_range = NULL,
                                num_permutations = 100,
                                keep_permutation_distribution = FALSE,
                                workers = 1,
                                overwrite = FALSE,
                                xloc = NULL,
                                yloc = NULL){
  #make sure that the range satisfies requirements for spatstat
  if(!(0 %in% r_range) & !is.null(r_range))
    r_range = c(0, r_range)
  #make sure that the mif object is, a mif object
  if(!inherits(mif, "mif"))
    stop("mIF should be of class `mif` created with function `createMIF()`\n\tTo check use `inherits(mif, 'mif')`")
  if(!inherits(mnames, 'character') & !inherits(mnames, 'data.frame'))
    stop("mnames must either be character vector or data.frame")
  #make marker table
  if(inherits(mnames, 'character')){
    mnames =  expand.grid(Anchor = mnames,
                          Counted = mnames) %>%
      dplyr::filter(Anchor != Counted)
  }
  #spatial file loop
  out = parallel::mclapply(mif$spatial, function(spat){
    if(FALSE %in% unique(unlist(mnames)) %in% colnames(spat))
      stop("Marker name not in spatial data")
    #get center of the cells
    if(is.null(xloc) | is.null(yloc)){
      spat = spat %>%
        dplyr::mutate(xloc = (XMax + XMin)/2,
                      yloc = (YMax + YMin)/2)
    } else {
      #rename columns to follow xloc and yloc names
      spat = spat %>%
        dplyr::rename("xloc" = !!xloc, 
                      "yloc" = !!yloc)
    }
    #select only needed columns
    spat = spat %>%
      dplyr::select(!!mif$sample_id, xloc, yloc, dplyr::any_of(!!unique(unlist(mnames))))
    #ppp window
    win = spatstat.geom::convexhull.xy(spat$xloc, spat$yloc)
    
    res = parallel::mclapply(1:nrow(mnames), function(marker){
      markers = unname(unlist(mnames[marker,])) %>% as.character()
      sample_ppp = spatstat.geom::ppp(spat$xloc, spat$yloc,
                                      window = win)
      cells = get_bi_rows(spat, markers)
      #if not enough cells
      if(length(unique(cells$Marker)) < 2){
        if(!keep_permutation_distribution){
          df = data.frame(From = markers[1], To = markers[2],
                          r = r_range)
        } else {
          df = data.frame(From = markers[1], To = markers[2],
                          r = r_range) %>%
            dplyr::full_join(expand.grid(r = r_range,
                                  iter = seq(num_permutations)))
        }
        return(df)
      }
      
      #subset ppp to contain right cells
      ps = subset(sample_ppp, cells$cell)
      spatstat.geom::marks(ps) = cells$Marker
      dists = spatstat.geom::crossdist(subset(ps, marks == markers[1]),
                                       subset(ps, marks == markers[2]))
      #so we first identify the nearest neighbor for each of the anchor/"from" cells
      #then check which range or bin the nearest neighbor distance is in
      #finally do cumulative sum and normalize to total number of from/to cells
      ITvals = data.frame(r = r_range,
                          `Observed Interaction` = c(0, cumsum(table(cut(apply(dists, 1, min), breaks=r_range, include.lowest = T)))/nrow(cells)*100) %>% unname(),
                          check.names = FALSE)
      
      obs = data.frame(r = r_range,
                            From = markers[1],
                            To = markers[2]) %>%
        #mutate(!!mif$sample_id := unique(spat[[mif$sample_id]])) %>%
        dplyr::full_join(ITvals, by = dplyr::join_by(r))
      
      perms = parallel::mclapply(seq(num_permutations), function(p){
        tmp_ps = subset(sample_ppp, sample(1:nrow(spat), ps$n, replace = FALSE))
        spatstat.geom::marks(tmp_ps) = cells$Marker
        tmp_dists = spatstat.geom::crossdist(subset(tmp_ps, marks == markers[1]),
                                         subset(tmp_ps, marks == markers[2]))
        data.frame(r = r_range,
                   From = markers[1],
                   To = markers[2],
                   `Permuted Interaction` = c(0, cumsum(table(cut(apply(tmp_dists, 1, min), breaks=r_range, include.lowest = T)))
                                              /nrow(cells)*100) %>% unname(),
                   iter = p,
                   check.names = FALSE)
      }, mc.preschedule = FALSE,
      mc.allow.recursive = TRUE) %>%
        do.call(dplyr::bind_rows, .)
      
      #calculate number of permutations with greater interaction than observed
      dat = dplyr::full_join(obs,
                             perms, by = dplyr::join_by(r, From, To)) %>%
        dplyr::group_by(From, To, r) %>%
        dplyr::mutate(`Permuted_larger_than_Observed` = sum(`Permuted Interaction` > unique(`Observed Interaction`), na.rm = TRUE))
      
      #collapse if not needing permutations
      if(!keep_permutation_distribution){
        dat %>%
          dplyr::select(-iter) %>%
          dplyr::group_by(Permuted_larger_than_Observed, .add = TRUE) %>%
          dplyr::summarise_all(mean, na.rm = TRUE)
      } else {
        dat
      }
      
    }, mc.preschedule = FALSE,
    mc.allow.recursive = TRUE) %>%
      do.call(dplyr::bind_rows, .)
    
    #add in sample ID
    res %>%
      dplyr::mutate(!!mif$sample_id := unique(spat[[mif$sample_id]]),
                    .before = 1)
  }, mc.preschedule = FALSE,
  mc.cores = workers,
  mc.allow.recursive = TRUE)
  
  #collapse and claculate degree of interaction from permutations
  out = out %>%
    do.call(dplyr::bind_rows, .)
  out$`Degree of Interaction Permuted` = out$`Observed Interaction` - out$`Permuted Interaction`
  
  #add back to mIF object
  if(overwrite | !exists("interaction_variable", where = mif$derived)){
    mif$derived$interaction_variable = out %>%
      dplyr::mutate(Run = 1)
  } else {
    n_run = max(mif$derived$interaction_variable$Run) + 1
    mif$derived$interaction_variable = dplyr::bind_rows(mif$derived_intraction_variable,
                                                        out %>% dplyr::mutate(Run = n_run))
  }
  return(mif)
}
