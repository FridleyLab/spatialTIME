#' Bivariate Pair Correlation Function
#'
#' @param mif object of class `mif`
#' @param mnames character vector or dataframe with 2 columns containing markers/marker combinations to run
#' @param r_range numeric vector radii to measure
#' @param num_permutations integer for the number of permutations to run
#' @param edge_correction character string for which edge correction to implement for Ripley's K
#' @param keep_permutation_distribution boolean whether to summarise the permutations or keep all
#' @param workers integer for number of cores to use when calculating
#' @param overwrite boolean for whether to overwrite existing bivariate pair correlation results
#' @param xloc x location column in spatial files
#' @param yloc y location column in spatial files
#' @param ... other variables to pass to `[spatstat.explore::pcfcross]`
#' 
#' @return `mif` object with the bivariate_pair_correlation slot filled
#' @export
#'
bi_pair_correlation = function(mif,
                               mnames,
                               r_range = NULL,
                               num_permutations= 100,
                               edge_correction = "translation",
                               keep_permutation_distribution = FALSE,
                               workers = 1,
                               overwrite = FALSE,
                               xloc = NULL,
                               yloc = NULL,
                               ...){
  #make sure that the range satisfies requirements for spatstat
  if(!(0 %in% r_range) & !is.null(r_range))
    r_range = c(0, r_range)
  #check edge correction
  if(length(edge_correction)!=1)
    stop("edge_correction must be of length 1")
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
  
  #over spatial failes
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
    #window
    win = spatstat.geom::convexhull.xy(spat$xloc, spat$yloc)
    
    #over markers
    res = parallel::mclapply(1:nrow(mnames), function(marker){
      #bivariate is pcfcross
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
      
      #ps = subset(sample_ppp, cells$cell)
      #spatstat.geom::marks(ps) = cells %>% arrange(cellid) %>% pull(Marker)
      ps = spatstat.geom::ppp(cells$xloc, cells$yloc, window = win, marks = cells$Marker)
      
      obs = spatstat.explore::pcfcross(ps, markers[1], markers[2],
                                  r = r_range, correction = edge_correction,
                                  ...) %>%
        data.frame() %>%
        rename("Observed g" = 3)
      
      #run permutations
      perms = parallel::mclapply(seq(num_permutations), function(p){
        tmp_ps = subset(sample_ppp, sample(1:nrow(spat), ps$n, replace = FALSE))
        spatstat.geom::marks(tmp_ps) = cells$Marker
        spatstat.explore::pcfcross(tmp_ps, markers[1], markers[2],
                              r = r_range, correction = edge_correction,
                              ...) %>%
          data.frame() %>%
          dplyr::rename("Permuted g" = 3) %>%
          mutate(iter = p)
      }, mc.preschedule = FALSE,
      mc.allow.recursive = TRUE) %>%
        do.call(dplyr::bind_rows, .)
      
      dat = dplyr::full_join(obs,
                             perms) %>%
        dplyr::rename("Theoretical g" = 2) %>%
        dplyr::mutate(From = markers[1], To = markers[2], .before = 1) %>%
        dplyr::group_by(From, To, r) %>%
        dplyr::mutate(`Permuted_larger_than_Observed` = sum(`Permuted g` > unique(`Observed g`), na.rm = TRUE))
      
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
    
    res %>%
      dplyr::mutate(!!mif$sample_id := unique(spat[[mif$sample_id]]),
                    .before = 1)
  }, mc.preschedule = FALSE,
  mc.cores = workers,
  mc.allow.recursive = TRUE)
  
  out = out %>%
    do.call(dplyr::bind_rows, .)
  out$`Degree of Correlation Theoretical` = out$`Observed g` - out$`Theoretical g`
  out$`Degree of Correlation Permuted` = out$`Observed g` - out$`Permuted g`
  
  if(overwrite | !exists("bivariate_pair_correlation", where = mif$derived)){
    mif$derived$bivariate_pair_correlation = out %>%
      dplyr::mutate(Run = 1)
  } else {
    n_run = max(mif$derived$bivaraite_pair_correlation$Run)+1
    mif$derived$bivaraite_pair_correlation = dplyr::bind_rows(mif$derived$bivaraite_pair_correlation,
                                                               out %>% dplyr::mutate(Run = n_run))
  }
  
  return(mif)
}
