#' Bivariate Nearest Neighbor G(r)
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
#'
#' @examples
#' x <- spatialTIME::create_mif(clinical_data = spatialTIME::example_clinical %>% 
#'   dplyr::mutate(deidentified_id = as.character(deidentified_id)),
#'   sample_data = spatialTIME::example_summary %>% 
#'   dplyr::mutate(deidentified_id = as.character(deidentified_id)),
#'   spatial_list = spatialTIME::example_spatial[1:2],
#'   patient_id = "deidentified_id", 
#'   sample_id = "deidentified_sample")
#'     
#' mnames_good <- c("CD3..Opal.570..Positive","CD8..Opal.520..Positive",
#'   "FOXP3..Opal.620..Positive","PDL1..Opal.540..Positive",
#'   "PD1..Opal.650..Positive","CD3..CD8.","CD3..FOXP3.")
#' \dontrun{
#' x2 = bi_NN_G(mif = x, mnames = mnames_good[1:2], 
#'       r_range = 0:100, num_permutations = 10, 
#'       edge_correction = "rs", keep_perm_dis = FALSE, 
#'       workers = 1, overwrite = TRUE)
#' }
#' 
#' @export
bi_NN_G = function(mif,
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
  if(length(mnames) == 1){
    stop("Please use univariate NN G(r) function for only a single marker")
  }
  if(length(mnames) == 0){
    stop("Need at least 2 markers provided to `mnames` in order to do bivariate NN G(r)")
  }
  #run bivar nn g
  out = pbmcapply::pbmclapply(mif$spatial, function(spat){
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
    #get name of sample, make spatial a matrix and build sample window
    core = unlist(spat[1, mif$sample_id])
    spat = spat %>%
      dplyr::select(xloc, yloc, dplyr::any_of(mnames)) %>% 
      as.matrix()
    win = spatstat.geom::convexhull.xy(spat[,"xloc"], spat[,"yloc"])
    areaW = spatstat.geom::area(win)
    #get the different marker combinations
    marker_combos = expand.grid(anchor = mnames,
                                counted = mnames) %>%
      dplyr::filter(anchor != counted) %>%
      dplyr::mutate_all(as.character)
    marker_list = split(unlist(marker_combos), rep(seq(nrow(marker_combos)), 2))
    #compute bivariate G(r)
    res = parallel::mclapply(marker_list, function(marks){ #marks is a charter vector of 2L
      df = spat[,c("xloc", "yloc", marks)]
      df = df[!(df[,marks[1]] == 1 & df[,marks[2]] == 1),]
      
      mark_tabs = colSums(df[,marks])
      if(TRUE %in% (mark_tabs < 3)){
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
      zeroes <- numeric(length(r_range))
      G_cross_df <- data.frame(r = r_range, 
                               `Theoretical CSR` = 1 - exp(-lamJ * pi * 
                                                           r_range^2),
                               check.names = F)
      dists = as.matrix(dist(df[,1:2]))
      diag(dists) = NA
      #for observed
      obs_nnd = apply(dists[df[,marks[1]] == 1,df[,marks[2]] == 1], 1, min, na.rm = T)
      obs_bdry = spatstat.geom::bdist.points(spatstat.geom::ppp(x = df[df[,marks[1]] == 1,"xloc"],
                                                            y = df[df[,marks[1]] == 1,"yloc"],
                                                            window = win))
      obs_d = (obs_nnd <= obs_bdry)
      
      if(edge_correction == "none"){
        G_cross_df = cbind(G_cross_df, `Observed G` = c(0,unname(cumsum(table(cut(obs_nnd, r_range)))/length(obs_nnd))))
        
        G_cross_df2 = lapply(seq(1:num_permutations), function(perm_num){
          df_rows = sample(1:nrow(df), sum(colSums(df[,marks])), replace = F)
          names(df_rows) = rep(names(mark_tabs), mark_tabs)
          #for permutation
          obs_nnd = apply(dists[df_rows[names(df_rows) == marks[1]],df_rows[names(df_rows) == marks[2]]], 1, min, na.rm = T)
          obs_bdry = spatstat.geom::bdist.points(spatstat.geom::ppp(x = df[df_rows[names(df_rows) == marks[1]],"xloc"],
                                                                    y = df[df_rows[names(df_rows) == marks[1]],"yloc"],
                                                                    window = win))
          obs_d = (obs_nnd <= obs_bdry)
          
          data.frame(r = r_range,
                     `Permuted G` = c(0,unname(cumsum(table(cut(obs_nnd, r_range)))/length(obs_nnd))),
                     iter = perm_num, check.names = F)
        }) %>%
          do.call(dplyr::bind_rows, .) %>%
          dplyr::full_join(G_cross_df, ., by = "r")
        
      } else if(edge_correction == "han"){
        x = obs_nnd[obs_d]
        a = spatstat.geom::eroded.areas(win, r_range)
        G = unname(cumsum(c(0, table(cut(x, r_range)))/a))
        G_cross_df = cbind(G_cross_df, `Observed G` = G/max(G[is.finite(G)]))
        #permutations
        G_cross_df2 = lapply(seq(1:num_permutations), function(perm_num){
          df_rows = sample(1:nrow(df), sum(colSums(df[,marks])), replace = F)
          names(df_rows) = rep(names(mark_tabs), mark_tabs)
          #for permutation
          obs_nnd = apply(dists[df_rows[names(df_rows) == marks[1]],df_rows[names(df_rows) == marks[2]]], 1, min, na.rm = T)
          obs_bdry = spatstat.geom::bdist.points(spatstat.geom::ppp(x = df[df_rows[names(df_rows) == marks[1]],"xloc"],
                                                                    y = df[df_rows[names(df_rows) == marks[1]],"yloc"],
                                                                    window = win))
          obs_d = (obs_nnd <= obs_bdry)
          x = obs_nnd[obs_d]
          a = spatstat.geom::eroded.areas(win, r_range)
          G = unname(cumsum(c(0, table(cut(x, r_range)))/a))
          
          data.frame(r = r_range,
                     `Permuted G` = G/max(G[is.finite(G)]),
                     iter = perm_num, check.names = F)
        }) %>%
          do.call(dplyr::bind_rows, .) %>%
          dplyr::full_join(G_cross_df, ., by = "r")
        
      } else if(edge_correction == "rs"){
        #for edge correction rs
        o <- pmin.int(obs_nnd, obs_bdry)
        result = spatstat.univar::km.rs(o, obs_bdry, obs_d,
                                        spatstat.geom::handle.r.b.args(r_range, NULL, W, 
                                                                       rmaxdefault = spatstat.explore::rmax.rule("G", win, lamJ)))
        G_cross_df = cbind(G_cross_df, `Observed G` = result$rs)
        #permutations
        G_cross_df2 = lapply(seq(1:num_permutations), function(perm_num){
          df_rows = sample(1:nrow(df), sum(colSums(df[,marks])), replace = F)
          names(df_rows) = rep(names(mark_tabs), mark_tabs)
          #for permutation
          obs_nnd = apply(dists[df_rows[names(df_rows) == marks[1]],df_rows[names(df_rows) == marks[2]]], 1, min, na.rm = T)
          obs_bdry = spatstat.geom::bdist.points(spatstat.geom::ppp(x = df[df_rows[names(df_rows) == marks[1]],"xloc"],
                                                                    y = df[df_rows[names(df_rows) == marks[1]],"yloc"],
                                                                    window = win))
          obs_d = (obs_nnd <= obs_bdry)
          o <- pmin.int(obs_nnd, obs_bdry)
          result = spatstat.univar::km.rs(o, obs_bdry, obs_d,
                                          spatstat.geom::handle.r.b.args(r_range, NULL, W, 
                                                                         rmaxdefault = spatstat.explore::rmax.rule("G", win, lamJ)))
          
          data.frame(r = r_range,
                     `Permuted G` = result$rs,
                     iter = perm_num, check.names = F)
        }) %>%
          do.call(dplyr::bind_rows, .) %>%
          dplyr::full_join(G_cross_df, ., by = "r")
      }
      
      return(G_cross_df2 %>%
               dplyr::mutate(!!mif$sample_id := core,
                             Anchor = marks[1],
                             Counted = marks[2],
                             .before = 1))
    }, mc.allow.recursive = T) %>%
      do.call(dplyr::bind_rows, .)
      
    if(keep_perm_dis){
      return(res)
    }
    
    res %>% 
      dplyr::select(-iter) %>%
      dplyr::group_by(across(1:4)) %>% 
      dplyr::summarise_all(~mean(., na.rm = T))
  }, mc.cores = workers, mc.preschedule = FALSE, mc.allow.recursive = T) %>%
    do.call(dplyr::bind_rows, .) %>%
    #calculate the degree of clustering from both the theoretical and permuted
    dplyr::mutate(`Degree of Clustering Permutation` = `Observed G` - `Permuted G`,
                  `Degree of Clustering Theoretical` = `Observed G` - `Theoretical CSR`)
  
  if(overwrite){
    mif$derived$bivariate_NN = out %>%
      dplyr::mutate(Run = 1)
  }
  if(!overwrite){
    mif$derived$bivariate_NN = mif$derived$bivariate_NN %>%
      dplyr::bind_rows(out %>%
                         dplyr::mutate(Run = ifelse(!exists("bivariate_NN", mif$derived),
                                                    1,
                                                    max(mif$derived$bivariate_NN$Run) + 1)))
  }
  return(mif)
}
