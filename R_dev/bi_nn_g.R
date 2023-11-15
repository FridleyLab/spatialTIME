#' Bivariate Nearest Neighbor G(r)
#'
#' @param mif object of class `mif` created by function `create_mif()`
#' @param mnames character vector of column names within the spatial files, indicating whether a cell row is positive for a phenotype
#' @param r_range numeric vector of radii around marker positive cells which to use for G(r)
#' @param num_permutations integer number of permutations to use for estimating core specific complete spatial randomness (CSR)
#' @param edge_correction character vector of edge correction methods to use: "rs", "km" or "han"
#' @param keep_permutation_distribution boolean for whether to summarise permutations to a single value or maintain each permutations result
#' @param workers integer number for the number of CPU cores to use in parallel to calculate all samples/markers
#' @param overwrite boolean whether to overwrite previous run of NN G(r) or increment "RUN" and maintain  previous measurements
#' @param xloc,yloc the x and y location columns in the spatial files that indicate the center of the respective cells
#' 
#' @return object of class `mif` containing a new slot under `derived` got nearest neighbor distances
#' @export
#'
#' @examples
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
#' x2 = bi_NN_G2(mif = x, mnames = mnames_good, r_range = 0:100, num_permutations = 25, edge_correction = "rs", keep_permutation_distribution = FALSE, workers = 1, overwrite = TRUE)
bi_NN_G2 = function(mif,
                 mnames,
                 r_range = 0:100,
                 num_permutations = 50,
                 edge_correction = "rs",
                 keep_permutation_distribution = FALSE,
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
  if(length(mnames) == 1){
    stop("Please use univariate NN G(r) function for only a single marker")
  }
  if(length(mnames) == 0){
    stop("Need at least 2 markers provided to `mnames` in order to do bivariate NN G(r)")
  }
  #run bivar nn g
  out = parallel::mclapply(mif$spatial, function(spat){
    if(is.null(xloc)){
      spat$xloc = (spat$XMax + spat$XMin)/2
    }
    if(is.null(yloc)){
      spat$yloc = (spat$YMax + spat$YMin)/2
    }
    #get name of sample, make spatial a matrix and build sample window
    core = spat[1, mif$sample_id]
    spat = spat %>%
      dplyr::select(xloc, yloc, any_of(mnames)) %>% 
      as.matrix()
    win = spatstat.geom::convexhull.xy(spat[,"xloc"], spat[,"yloc"])
    areaW = spatstat.geom::area(win)
    #get the different marker combinations
    marker_combos = expand.grid(anchor = markers,
                                counted = markers) %>%
      dplyr::filter(anchor != counted) %>%
      dplyr::mutate_all(as.character)
    marker_list = split(unlist(marker_combos), rep(seq(nrow(marker_combos)), 2))
    #compute bivariate G(r)
    res = parallel::mclapply(marker_list, function(marks){ #marks is a charter vector of 2L
      df = spat[,c("xloc", "yloc", marks)]
      df = df[df[[1]] != 1 & df[[2]] != 1,]
      
      if(TRUE %in% (colSums(df[,marks]) < 3)){
        perms = data.frame(samp = core,
                           Anchor = unname(marks[1]),
                           Counted = unname(marks[2]),
                           r = r_range,
                           `Theoretical CSR` = NA,
                           `Permuted G` = NA,
                           `Observed G` = NA,
                           check.names = F) %>%
          dplyr::full_join(expand.grid(r = r_range, iter = seq(num_permutations)), by = "r")
        colnames(perms)[1] = mif$sample_id
        return(perms)
      }
      
      #total number of points
      npts = nrow(df)
      lamJ = sum(df[,marks[2]])/areaW
      rmax = max(r_range)
      
      pp_obj = spatstat.geom::ppp(x = tmp[,"xloc"],
                                  y = tmp[,"yloc"],
                                  window= win,
                                  marks = tmp[,"Marks"])
      #area = spatstat.geom::area(win)
      #n = spatstat.geom::npoints(pp_obj)
      #theo = 1 - exp(-(n/area) * pi * r^2)
      #dists = as.matrix(dist(df[,1:2]))
      #dists[dists == 0] = NA
      #nnd = apply(dists, 1, min, na.rm=T)
      #d = (nnd <= spatstat.geom::bdist.points(pp_obj))
      #x = nnd[d]
      #a=spatstat.geom::eroded.areas(win, r_range)
      #this is for hanish
      #G = cumsum(hist(x[x<=max(r_range)],plot =F, breaks = r_range)$counts/100)
      #G/max(G[is.finite(G)])
      G = spatstat.explore::nearest.neighbour(pp_obj, r = r_range, correction = edge_correction) %>%
        data.frame() %>%
        dplyr::rename("Theoretical G" = 2, "Observed G" = 3)
      
      if(permute){
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
      } else {
        perms = data.frame(r = r_range,
                           `Theoretical G` = 1 - exp(-(nrow(spat)/spatstat.geom::area(win)) * pi * r_range^2))
      }
      
      return(perms)
    }) %>%
      do.call(dplyr::bind_rows, .) %>%
      dplyr::mutate(!!mif$sample_id := core)
    res = res[,c(7, 6, 4, 1, 2, 5, 3)]
    if(keep_permutation_distribution){
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
