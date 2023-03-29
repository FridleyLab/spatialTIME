#' Calculate Ripley's K 
#'
#' @param mif object of class `mif` created with `create_mif`
#' @param mnames cell phenotype markers to calculate Ripley's K for
#' @param r_range radius range (including 0)
#' @param num_permutations number of permutations to use to estimate CSR. If `keep_perm_dis` is set to FALSE, this will be ignored
#' @param edge_correction edge correction method to pass to `Kest`. can take one of "translation", "isotropic", "none"
#' @param method not used currently
#' @param permute whether to use CSR estimate or use permutations to determine CSR
#' @param keep_permutation_distribution whether to find mean of permutation distribution or each
#' permutation calculation
#' @param workers number of cores to use for calculations
#' @param overwrite whether to overwrite the `univariate_Count` slot within `mif$derived`
#' @param xloc the location of the center of cells. If left `NULL`, `XMin`, `XMax`, `YMin`, and `YMax` must be present.
#' @param yloc the location of the center of cells. If left `NULL`, `XMin`, `XMax`, `YMin`, and `YMax` must be present.
#' @param big the number of cells at which to flip from an edge correction method other than 'none' to 'none' due to size
#' 
#' @description 
#' ripleys_k() calculates the emperical Ripley's K measurement for the cell types specified by mnames in the mIF object. This
#' is very useful when exploring the spatial clustering of single cell types on TMA cores or ROI spots following proccessing
#' with a program such as HALO for cell phenotyping.
#' 
#' In the `ripleys_k` function, there is the ability to perform permutations in order to assess whether the clustering
#' of a cell type is significant, or the ability to derive the exact CSR and forgo permutations for much faster sample
#' processing. Permutations can be helpful if the significance of clustering wasnts to be identified - run 1000 permutations 
#' and if observed is outside 95-percentile then significant clustering. We, however, recommend using the exact CSR estimate
#' due to speed.
#' 
#' Some things to be aware of when computing the exact Ripley's K estimate, if your spatial file is greater than 
#' the `big` size, the edge correction will be converted to 'none' in order to save on resources and compute time. 
#' Due to the introduction of Whole Slide Imaging (WSI), this can easily be well over 1,000,000 cells, and calculating 
#' edge correction for these spatial files will not succeed when attempting to force an edge correction on it.
#' 
#' @return object of class `mif`
#' @export
#'
#' @examples
#' x <- spatialTIME::create_mif(clinical_data =spatialTIME::example_clinical %>% 
#'   dplyr::mutate(deidentified_id = as.character(deidentified_id)),
#'   sample_data = spatialTIME::example_summary %>% 
#'   dplyr::mutate(deidentified_id = as.character(deidentified_id)),
#'   spatial_list = spatialTIME::example_spatial,
#'   patient_id = "deidentified_id", 
#'   sample_id = "deidentified_sample")
#' mnames = x$spatial[[1]] %>%
#'   colnames() %>%
#'   grep("Pos|CD", ., value =TRUE) %>%
#'   grep("Cyto|Nucle", ., value =TRUE, invert =TRUE)
#' x2 = ripleys_k(mif = x, 
#'   mnames = mnames[1], 
#'   r_range = seq(0, 100, 1), 
#'   num_permutations = 100,
#'   edge_correction = "translation", 
#'   method = "K", 
#'   permute = FALSE,
#'   keep_permutation_distribution =FALSE, 
#'   workers = 1, 
#'   overwrite =TRUE)
ripleys_k = function(mif, 
                     mnames,
                     r_range = seq(0, 100, 1), 
                     num_permutations = 50, 
                     edge_correction = "translation",
                     method = "K", 
                     permute = FALSE,
                     keep_permutation_distribution = FALSE,
                     workers = 1,
                     overwrite = FALSE,
                     xloc = NULL,
                     yloc = NULL,
                     big = 10000){
  
  if(keep_permutation_distribution == TRUE & permute == FALSE){
    stop("Conflicting `perm` and `keep_permutation_distribution` parameters. Using estimate\n
         \tIf wanting to use permutatations, set `permute = TRUE`")
  }
  #check if 0 is contained within the range of r
  if(!(0 %in% r_range)){
    r_range = c(0, r_range)
  }
  if(!inherits(mif, "mif")){
    stop("mIF should be of class `mif` created with function `createMIF()`\n\tTo check use `inherits(mif, 'mif')`")
  }
  
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
    if(nrow(spat) > big){
      edge_correction = 'none'
    }
    #select only needed columns
    spat = spat %>%
        dplyr::select(!!mif$sample_id, xloc, yloc, !!mnames)
    #window
    win = spatstat.geom::convexhull.xy(spat$xloc, spat$yloc)
    #begin calculating the ripley's k and the permutation/estimate
    if(nrow(spat)<10000 & permute == TRUE){ #have to calculate the permutation distribution
      dists = as.matrix(dist(spat[,c("xloc", "yloc")]))
      area = spatstat.geom::area(win)
      
      #calculate the edge corrections
      if(edge_correction %in% c("translation", "trans")){
        edge = spatstat.explore::edge.Trans(spatstat.geom::ppp(x = spat$xloc, y = spat$yloc, window = win), W = win)
      } else if(edge_correction %in% c("isotropic", "iso")){
        edge = spatstat.explore::edge.Ripley(spatstat.geom::ppp(x = spat$xloc, y = spat$yloc, window = win))
      } else if(edge_correction == "none"){
        edge = matrix(nrow = nrow(spat), ncol = nrow(spat), data = 1)
      }
      
      #dists = fast_mm(dists, edge)
      res = parallel::mclapply(mnames, function(marker){
        #find the rows that are positive for marker
        pos = which(spat[marker] == 1)
        #if there are less than 3 cells then just return NA table because cannot calculate
        if(length(pos) < 3){
          d = data.frame(iter = as.character(seq(num_permutations)),
                         Label = spat[1,mif$sample_id],
                         Marker = marker,
                         `Observed K` = NA,
                         `Permuted CSR` = NA,
                         `Exact CSR` = NA,
                         check.names =FALSE)
          d = dplyr::full_join(d, expand.grid(iter = as.character(seq(num_permutations)),
                                              r = r_range), by = "iter") %>%
            dplyr::mutate(`Theoretical CSR` = pi * r^2)
          return(d)
        }
        
        edge_pos = edge[pos, pos]
        counts = sapply(r_range, function(r) sum(edge_pos[which(dists[pos, pos] > 0 & dists[pos, pos] < r)]))
        obs = (counts * area)/(length(pos)*(length(pos)-1))
        theo = pi * r_range^2
        
        perms = parallel::mclapply(seq(num_permutations), function(x){
          n_pos = sample(1:nrow(spat), length(pos), replace =FALSE)
          edge_pos = edge[n_pos, n_pos]
          counts = sapply(r_range, function(r) sum(edge_pos[which(dists[n_pos, n_pos] > 0 & dists[n_pos, n_pos] < r)]))
          permed = (counts * area)/(length(pos)*(length(pos)-1))
          theo = pi * r_range^2
          return(data.frame(iter = as.character(x),
                            r = r_range,
                            `Theoretical CSR` = theo,
                            `Permuted CSR` = permed,
                            `Exact CSR` = NA,
                            check.names = FALSE))
        }, mc.allow.recursive = TRUE, mc.preschedule = FALSE) %>%
          do.call(dplyr::bind_rows, .)
        final = dplyr::full_join(data.frame(r = r_range,
                                            `Theoretical CSR` = theo,
                                            `Observed K` = obs,
                                            check.names = FALSE),
                                 perms, by = c("r", "Theoretical CSR"))
        final$Marker = marker
        final$Label = spat[1,1] #hard coded for example
        final[,c(4,8, 7, 1,2,3,5,6)]
      }, mc.allow.recursive = TRUE, mc.preschedule = FALSE) %>% #collapse all markers for spat
        do.call(dplyr::bind_rows, .) %>%
        dplyr::mutate(`Degree of Clustering Permutation` = `Observed K` - `Permuted CSR`,
                      `Degree of Clustering Theoretical` = `Observed K` - `Theoretical CSR`,
                      `Degree of Clustering Exact` = `Observed K` - `Exact CSR`)
    } else if(permute == TRUE & nrow(spat)>10000){ #if we need perm distribution but too many cells
      res = parallel::mclapply(mnames, function(marker){
        #select the center of cells and marker column
        dat = spat %>%
          dplyr::select(xloc, yloc, !!marker)
        #select columns that are postive for marker
        dat2 = dat %>% dplyr::filter(get(marker) != 0)
        #calculate the observed K
        kobs = spatstat.geom::ppp(dat2$xloc, dat2$yloc, window = win) %>%
          spatstat.explore::Kest(r = r_range, correction = edge_correction) %>%
          data.frame() %>%
          dplyr::rename("Theoretical CSR" = 2, "Observed K" = 3) %>%
          dplyr::mutate(Label = unique(spat[[mif$sample_id]]),
                        Marker = marker, .before=1)
        #calculate the permuted CSR
        kperms = parallel::mclapply(seq(num_permutations), function(perm){
          dat2 = dat[sample(seq(nrow(dat)), sum(dat[[marker]]), replace=F),]
          spatstat.geom::ppp(dat2$xloc, dat2$yloc, window = win) %>%
            spatstat.explore::Kest(r = r_range, correction = edge_correction) %>%
            data.frame() %>%
            dplyr::rename("Theoretical CSR" = 2, "Permuted CSR" = 3) %>%
            dplyr::mutate(Label = unique(spat[[mif$sample_id]]),
                          Marker = marker,.before=1) %>%
            dplyr::mutate(iter = as.character(perm), .before = 1)
        }, mc.allow.recursive = TRUE, mc.preschedule = FALSE) %>%
          do.call(dplyr::bind_rows, .) %>% 
          mutate(`Exact CSR` = NA)
        K = dplyr::full_join(kobs, kperms,
                             by = c("Label", "Marker", "r", "Theoretical CSR")) %>%
          relocate(iter, .before = 1)
      }, mc.allow.recursive = TRUE, mc.preschedule = FALSE) %>%
        do.call(dplyr::bind_rows, .) %>%
        dplyr::mutate(`Degree of Clustering Permutation` = `Observed K` - `Permuted CSR`,
                      `Degree of Clustering Theoretical` = `Observed K` - `Theoretical CSR`,
                      `Degree of Clustering Exact` = `Observed K` - `Exact CSR`)
    } else if(permute == FALSE){
      marker_res = parallel::mclapply(mnames, function(marker){
        #select the center of cells and marker column
        dat = spat %>%
          dplyr::select(xloc, yloc, !!marker)
        #select columns that are postive for marker
        dat2 = dat %>% dplyr::filter(get(marker) != 0)
        
        if(nrow(dat2) < 3){
          return(data.frame(iter = "Estimator", 
                            Label = unique(spat[[mif$sample_id]]),
                            Marker = marker,
                            r = r_range,
                            `Theoretical CSR` = pi * r_range^2,
                            `Observed K` = NA,
                            check.names =FALSE))
        }
        #calculate the observed K
        kobs = spatstat.geom::ppp(dat2$xloc, dat2$yloc, window = win) %>%
          spatstat.explore::Kest(r = r_range, correction = edge_correction) %>%
          data.frame() %>%
          dplyr::rename("Theoretical CSR" = 2, "Observed K" = 3) %>%
          dplyr::mutate(Label = unique(spat[[mif$sample_id]]),
                        Marker = marker, .before=1,
                        r = round(r))
        return(kobs)
      }, mc.allow.recursive = TRUE, mc.preschedule = FALSE) %>% #collapse all markers for spat
        do.call(dplyr::bind_rows, .)
      
      if(nrow(spat) < big){
        pp_obj = spatstat.geom::ppp(x = spat$xloc, y = spat$yloc, window = win)
        k_est = spatstat.explore::Kest(pp_obj, r = r_range, correction = edge_correction)
        gc(full=T)
        k_est2 = k_est %>%
          data.frame(check.names = FALSE) %>%
          dplyr::select(1,3) %>%
          dplyr::rename("Exact CSR" = 2) %>%
          dplyr::mutate(r = round(r),
                        `Permuted CSR` = NA, .before = `Exact CSR`)
      } else {
        ns = nrow(spat)
        slide = ceiling(ns / big)
        ranges = getTile(slide = slide, l = ns, size = big)
        counts = parallel::mclapply(ranges, function(i_range){
          parallel::mclapply(ranges, function(j_range){
              i_tmp = spatstat.geom::ppp(x = spat$xloc[i_range], y = spat$yloc[i_range], window = win)
              j_tmp = spatstat.geom::ppp(x = spat$xloc[j_range], y = spat$yloc[j_range], window = win)
              dists = spatstat.geom::crossdist(i_tmp, j_tmp)
              dists[dists == 0 | dists > max(r_range)] = NA
            
              if(edge_correction == "none"){
                counts = cumsum(spatstat.geom::whist(dists, 
                                                     spatstat.geom::handle.r.b.args(r_range, breaks=NULL, win, rmaxdefault = max(r_range))$val))
              } else {            
                
                edge = spatstat.explore::edge.Trans(i_tmp, j_tmp)
                counts = sapply(r_range, function(r){sum(edge[which(dists < r)])})
              }
              rm(dists)
            return(counts)
          }, mc.allow.recursive = TRUE, mc.preschedule = FALSE)%>%
            do.call(cbind, .) %>%
            rowSums()
        }, mc.allow.recursive = TRUE, mc.preschedule = FALSE)
        counts = counts %>%
          do.call(cbind, .) %>%
          rowSums()
        k = (counts * spatstat.geom::area(win))/(ns * (ns-1))
        k_est2 = data.frame(r = r_range,
                            `Permuted CSR` = NA,
                            `Exact CSR` = k, check.names = FALSE)
      }
      
      res = dplyr::full_join(marker_res, k_est2) %>%
        dplyr::mutate(iter = "Estimater", .before = 1) %>%
        dplyr::mutate(`Degree of Clustering Permutation` = NA,
                      `Degree of Clustering Theoretical` = `Observed K` - `Theoretical CSR`,
                      `Degree of Clustering Exact` = `Observed K` - `Exact CSR`)
      
    }
    #return final results
    return(res)
  }, mc.cores = workers, mc.allow.recursive =TRUE, mc.preschedule =FALSE) %>% #set mclapply params
    do.call(dplyr::bind_rows, .)
  out = out%>% #collapse all samples to single data frame
    dplyr::rename(!!mif$sample_id := Label)
  
  if(permute == TRUE & keep_permutation_distribution == FALSE){
    out = out %>%
      dplyr::mutate(iter = "Permuted") %>%
      dplyr::group_by(iter, across(mif$sample_id), Marker, r) %>%
      dplyr::summarise_all(~mean(., na.rm=TRUE))
  }
  
  if(overwrite){ #overwrite existing data in univariate_Count slot
    mif$derived$univariate_Count = out %>%
      dplyr::mutate(Run = 1)
  } else { #dont overwrite
    mif$derived$univariate_Count = mif$derived$univariate_Count %>%
      dplyr::bind_rows(out %>%
                         dplyr::mutate(Run = ifelse(exists("univariate_Count", mif$derived),
                                                    max(mif$derived$univariate_Count$Run) + 1,
                                                    1)))
  }
  
  return(mif)
}
