#' Univariate Nearest Neighbor G(r)
#'
#' @param mif object of class `mif` created by function `create_mif()`
#' @param mnames character vector of column names within the spatial files, indicating whether a cell row is positive for a phenotype
#' @param r_range numeric vector of radii around marker positive cells which to use for G(r)
#' @param num_permutations integer number of permutations to use for estimating core specific complete spatial randomness (CSR)
#' @param edge_correction character vector of edge correction methods to use: "rs", "km" or "han"
#' @param keep_perm_dis boolean for whether to summarise permutations to a single value or maintain each permutations result
#' @param workers integer number for the number of CPU cores to use in parallel to calculate all samples/markers
#' @param overwrite boolean whether to overwrite previous run of NN G(r) or increment "RUN" and maintain  previous measurements
#' @param xloc,yloc the x and y location columns in the spatial files that indicate the center of the respective cells
#' 
#' @return object of class `mif` containing a new slot under `derived` got nearest neighbor distances
#' @export
#'
#' @examples
#' library(dplyr)
#' x <- spatialTIME::create_mif(clinical_data = spatialTIME::example_clinical %>% 
#'   dplyr::mutate(deidentified_id = as.character(deidentified_id)),
#'   sample_data = spatialTIME::example_summary %>% 
#'   dplyr::mutate(deidentified_id = as.character(deidentified_id)),
#'   spatial_list = spatialTIME::example_spatial,
#'   patient_id = "deidentified_id", 
#'   sample_id = "deidentified_sample")
#'     
#' mnames_good <- c("CD3..Opal.570..Positive","CD8..Opal.520..Positive",
#'   "FOXP3..Opal.620..Positive","PDL1..Opal.540..Positive",
#'   "PD1..Opal.650..Positive","CD3..CD8.","CD3..FOXP3.")
#'   
#' x2 = NN_G(mif = x, mnames = mnames_good[1:2], 
#' r_range = 0:100, num_permutations = 10, 
#' edge_correction = "rs", keep_perm_dis = FALSE, 
#' workers = 1, overwrite = TRUE)
NN_G = function(mif,
                 mnames,
                 r_range = 0:100,
                 num_permutations = 50,
                 edge_correction = "rs",
                 keep_perm_dis = FALSE,
                 workers = 1,
                 overwrite = FALSE,
                 xloc = NULL,
                 yloc = NULL){
  if(!inherits(mif, "mif")){
    stop("Please submit a mif created with `create_mif`")
  }
  if(!(0 %in% r_range)){
    r_range = c(0, r_range)
  }
  out = parallel::mclapply(mif$spatial, function(spat){
    if(is.null(xloc)){
      spat$xloc = (spat$XMax + spat$XMin)/2
    } else {
      spat$xloc = spat[[xloc]]
    }
    if(is.null(yloc)){
      spat$yloc = (spat$YMax + spat$YMin)/2
    } else {
      spat$yloc = spat[[yloc]]
    }
    core = spat[1, mif$sample_id]
    spat = spat %>%
      dplyr::select(xloc, yloc, dplyr::any_of(mnames)) %>% 
      as.matrix()
    win = spatstat.geom::convexhull.xy(spat[,"xloc"], spat[,"yloc"])
    #using spatstat
    exact_G = spatstat.explore::nearest.neighbour(spatstat.geom::ppp(x = spat[,"xloc"],
                                                                     y = spat[,"yloc"],
                                                                     window = win),
                                                  r = r_range, correction = edge_correction) %>%
      data.frame()
    
    res = parallel::mclapply(mnames, function(marker){
      df = spat[spat[,marker] == 1,]
      if(sum(spat[,marker] == 1) < 3){
        perms = data.frame(r = r_range,
                           `Theoretical G` = NA,
                           `Permuted G` = NA, check.names = F) %>%
          dplyr::full_join(expand.grid(r = r_range, iter = seq(num_permutations)), by = "r") %>%
          mutate(`Observed G` = NA,
                 Marker = marker)
        return(perms)
      }
      pp_obj = spatstat.geom::ppp(x = df[,"xloc"],
                                  y = df[,"yloc"],
                                  window= win)
      
      G = spatstat.explore::nearest.neighbour(pp_obj, r = r_range, correction = edge_correction) %>%
        data.frame() %>%
        dplyr::rename("Theoretical G" = 2, "Observed G" = 3)
      
      perms = parallel::mclapply(seq(num_permutations), function(n){
        df = spat[sample(1:nrow(spat), nrow(df), replace =F),]
        pp_obj = spatstat.geom::ppp(x = df[,"xloc"],
                                    y = df[,"yloc"],
                                    window= win)
        spatstat.explore::nearest.neighbour(pp_obj, r = r_range, correction = edge_correction) %>%
          data.frame() %>%
          dplyr::rename("Theoretical G" = 2, "Permuted G" = 3) %>%
          dplyr::mutate(iter = n)
      }) %>%
        do.call(dplyr::bind_rows, .) %>%
        dplyr::full_join(G, by = c("r", "Theoretical G")) %>%
        dplyr::mutate(Marker = marker)
      
      return(perms)
    }) %>%
      do.call(dplyr::bind_rows, .) %>%
      dplyr::mutate(!!mif$sample_id := core)
    res = res[,c(7, 6, 4, 1, 2, 5, 3)]
    if(keep_perm_dis){
      return(res)
    }
    res = res[,-3]
    res %>% dplyr::group_by(across(1:3)) %>% dplyr::summarise_all(~mean(., na.rm = T))
  }, mc.cores = workers, mc.preschedule = FALSE, mc.allow.recursive = T) %>%
    do.call(dplyr::bind_rows, .) %>%
    #calculate the degree of clustering from both the theoretical and permuted
    dplyr::mutate(`Degree of Clustering Permutation` = `Observed G` - `Permuted G`,
                  `Degree of Clustering Theoretical` = `Observed G` - `Theoretical G`)
  
  if(overwrite){
    mif$derived$univariate_NN = out %>%
      dplyr::mutate(Run = 1)
  }
  if(!overwrite){
    mif$derived$univariate_NN = mif$derived$univariate_NN %>%
      dplyr::bind_rows(out %>%
                         dplyr::mutate(Run = ifelse(!exists("univariate_NN", mif$derived),
                                                    1,
                                                    max(mif$derived$univariate_NN$Run) + 1)))
  }
  return(mif)
}
