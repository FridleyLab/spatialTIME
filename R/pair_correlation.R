#' Univariate Pair Correlation Function
#' 
#' Implementation of the univariate pair correlation function from spatstat
#'
#' @param mif object of class `mif`
#' @param mnames character vector of marker names
#' @param r_range numeric vector including 0. If ignored, `spatstat` will decide range
#' @param num_permutations integer indicating how many permutations to run to determine CSR estimate
#' @param edge_correction character string of edge correction to apply to Ripley's K estimation
#' @param keep_permutation_distribution boolean for whether to keep the permutations or not
#' @param workers integer for number of threads to use when calculating metrics
#' @param overwrite boolean whether to overwrite existing results in the univariate_pair_correlation slot
#' @param xloc column name of single x value
#' @param yloc column name of single y value
#' @param ... other parameters to provide `spatstat::pcf` 
#' 
#' The Pair Correlation Function uses the derivative of Ripley's K so it does take slightly longer to calculate
#' 
#' `xloc` and `yloc`, if NULL, will be calculated from columns `XMax`, `XMin`, `YMax`, and `YMin`.
#'
#' @return mif object with with the univariate_pair_correlation derived slot filled or appended to
#' @export
#'
pair_correlation = function(mif,
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
  if(!(0 %in% r_range) & !is.null(r_range)){
    r_range = c(0, r_range)
  }
  
  #check edge correction
  if(length(edge_correction)!=1){
    stop("edge_correction must be of length 1")
  }
  #make sure that the mif object is, a mif object
  if(!inherits(mif, "mif")){
    stop("mIF should be of class `mif` created with function `createMIF()`\n\tTo check use `inherits(mif, 'mif')`")
  }
  #over spatial failes
  out = parallel::mclapply(mif$spatial, function(spat){
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
      dplyr::select(!!mif$sample_id, xloc, yloc, !!mnames)
    #window
    win = spatstat.geom::convexhull.xy(spat$xloc, spat$yloc)
    #over markers
    res = parallel::mclapply(mnames, function(marker){
      #bivariate is pcfcross
      sample_ppp = spatstat.geom::ppp(spat$xloc, spat$yloc,
                                      window = win)
      ps = subset(sample_ppp, spat[[marker]] == 1)
      obs = spatstat.explore::pcf(subset(sample_ppp, spat[[marker]] == 1),
                r = r_range, correction = edge_correction, ...) %>%
        data.frame() %>%
        rename("Observed g" = 3)
      
      #run permutations
      perms = parallel::mclapply(seq(num_permutations), function(p){
        spatstat.explore::pcf(subset(sample_ppp, as.logical(sample(spat[[marker]]))),
                              r = r_range, correction = edge_correction, ...) %>%
          data.frame() %>%
          dplyr::rename("Permuted g" = 3) %>%
          dplyr::mutate(iter = p)
      }, mc.preschedule = FALSE,
      mc.allow.recursive = TRUE) %>%
        do.call(dplyr::bind_rows, .)
      
      dat = dplyr::full_join(obs,
                       perms) %>%
        dplyr::rename("Theoretical g" = 2) %>%
        dplyr::mutate(Marker = marker, .before = 1) %>%
        dplyr::group_by(Marker, r) %>%
        dplyr::mutate(`Permuted_larger_than_Observed` = sum(`Permuted g` > unique(`Observed g`), na.rm = TRUE))
      
      #collapse if not needing permutations
      if(!keep_permutation_distribution){
        dat %>%
          select(-iter) %>%
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
  
  if(overwrite | !exists("univariate_pair_correlation", where = mif$derived)){
    mif$derived$univariate_pair_correlation = out %>%
      dplyr::mutate(Run = 1)
  } else {
    n_run = max(mif$derived$univaraite_pair_correlation$Run)+1
    mif$derived$univaraite_pair_correlation = dplyr::bind_rows(mif$derived$univaraite_pair_correlation,
                                                               out %>% mutate(Run = n_run))
  }
  
  return(mif)
}
