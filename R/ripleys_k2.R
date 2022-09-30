#' Calculate Ripley's K 
#'
#' @param mif object of class `mif` created with `create_mif`
#' @param mnames cell phenotype markers to calculate Ripley's K for
#' @param r_range radius range (including 0)
#' @param num_permutations number of permutations to use to estimate CSR. If `keep_perm_dis` is set to FALSE, this will be ignored
#' @param edge_correction edge correction method to pass to `Kest`. can take one of "translation", "isotropic", "none"
#' @param method not used currently
#' @param keep_perm_dis whether to keep the permutation results. Setting to FALSE will cause skipping of permutations and use CSR estimate
#' @param workers number of cores to use for calculations
#' @param overwrite whether to overwrite the `univariate_Count` slot within `mif$derived`
#' @param xloc, yloc the location of the center of cells. If left `NULL`, `XMin`, `XMax`, `YMin`, and `YMax` must be present.
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
#'   grep("Pos|CD", ., value = T) %>%
#'   grep("Cyto|Nucle", ., value = T, invert = T)
#' x2 = ripleys_k2(mif = x, 
#'   mnames = mnames, 
#'    r_range = seq(0, 100, 1), 
#'     num_permutations = 100,
#'      edge_correction = "translation", 
#'       method = "K", 
#'       keep_perm_dis = F, 
#'       workers = 1, 
#'        overwrite = T)
ripleys_k2 = function(mif, 
                     mnames,
                     r_range = seq(0, 100, 1), 
                     num_permutations = 50, 
                     edge_correction = "translation",
                     method = "K", 
                     keep_perm_dis = FALSE,
                     workers = 1,
                     overwrite = F,
                     xloc = NULL,
                     yloc = NULL){
  
  if(keep_perm_dis == FALSE){
    message("With 'keep_perm_dis' equal to FALSE, CSR Estimator will be used.")
    message("If wanting to specify permutations, set 'keep_perm_dis' to TRUE.")
  }
  #check if 0 is contained within the range of r
  if(!(0 %in% r_range)){
    r_range = c(0, r_range)
  }
  
  out = parallel::mclapply(mif$spatial, function(spat){
    #get center of the cells
    if(is.null(xloc) | is.null(yloc)){
      spat = spat %>%
        dplyr::mutate(xloc = (XMax + XMin)/2,
                      yloc = (YMax + YMin)/2)
    }
    #select only needed columns
    spat = spat %>%
        dplyr::select(!!mif$sample_id, xloc, yloc, !!mnames)
    #window
    win = spatstat.geom::convexhull.xy(spat$xloc, spat$yloc)
    #begin calculating the ripley's k and the permutation/estimate
    if(nrow(spat)<10000 & keep_perm_dis == T){ #have to calculate the permutation distribution
      #get distances between cells and area
      dists = as.matrix(dist(spat[,c("xloc", "yloc")]))
      area = spatstat.geom::area(win)
      
      #calculate the edge corrections
      if(edge_correction %in% c("translation", "trans")){
        edge = spatstat.core::edge.Trans(spatstat.geom::ppp(x = spat$xloc, y = spat$yloc, window = win), W = win)
      } else if(edge_correction %in% c("isotropic", "iso")){
        edge = spatstat.core::edge.Ripley(spatstat.geom::ppp(x = spat$xloc, y = spat$yloc, window = win))
      } else if(edge_correction == "none"){
        edge = matrix(nrow = nrow(spat), ncol = nrow(spat), data = 1)
      }
      
      #dists = fast_mm(dists, edge)
      res = parallel::mclapply(mnames, function(marker){
        print(marker)
        #find the rows that are positive for marker
        pos = which(spat[marker] == 1)
        #if there are less than 3 cells then just return NA table because cannot calculate
        if(length(pos) < 3){
          d = data.frame(iter = as.character(seq(num_permutations)), 
                     Label = spat[1,"deidentified_sample"],
                     Marker = marker,
                     `Observed K` = NA,
                     `Permuted K` = NA,
                     check.names = F)
          d = dplyr::full_join(d, expand.grid(iter = as.character(seq(num_permutations)),
                                              r = r_range)) %>%
            dplyr::mutate(`Theoretical K` = pi * r^2)
          return(d)
        }
        #for each r
        parallel::mclapply(r_range, function(r){
          #if r = 0 return 0s for cells - cannot be positive for itself
          if(r == 0){
            return(data.frame(iter = as.character(seq(num_permutations)), 
                              label = spat[1,"deidentified_sample"],
                              Marker = marker,
                              r = 0,
                              theo = 0,
                              obs = 0,
                              permed = 0))
          }
          #look to see if other true positives are within r and not itself
          in_range = dists[pos, pos] <= r & dists[pos,pos]>0
          #adjust positives with the corresponding weights
          in_range_adj = (in_range * 1) * edge[pos, pos]
          #count all adjusted positives
          counts = sum(in_range_adj)
          #calculate K
          obs = (counts * area)/(length(pos)*(length(pos)-1))
          theo = pi * r^2
          #prep permutations by selecting rnadom rows of positive length
          perms = sapply(seq(num_permutations), function(x){
            sample(1:nrow(spat), length(pos), replace = F)
          }) %>% t() %>% data.frame() %>% dplyr::distinct() %>% t()
          #perform above calculations on permuted positives
          permed = apply(perms, 2, function(x){
            in_range = dists[x, x] <= r & dists[x, x]>0
            in_range_adj = (in_range * 1) * edge[x, x]
            counts = sum(in_range_adj)
            (counts * area)/(length(x)*(length(x)-1))
          })
          #return permuted ripk for marker at r
          return(data.frame(iter = as.character(seq(num_permutations)),
                            label = spat[1, "deidentified_sample"],
                            Marker = marker,
                            r = r,
                            theo = theo,
                            obs = obs,
                            permed = permed))
        }) %>% #collapse all rs for marker
          do.call(dplyr::bind_rows, .) %>%
          dplyr::rename(Label = 2, `Theoretical K` = 5, `Observed K` = 6, `Permuted K` = 7)
      }) %>% #collapse all markers for spat
        do.call(dplyr::bind_rows, .) %>%
        dplyr::mutate(`Degree of Clustering Permutation` = `Observed K` - `Permuted K`,
                      `Degree of Clustering Theoretical` = `Observed K` - `Theoretical K`)
    } else if(nrow(spat) >= 10000 & keep_perm_dis == T){ #if we need perm distribution but too many cells
      res = parallel::mclapply(mnames, function(marker){
        #select the center of cells and marker column
        dat = spat %>%
          dplyr::select(xloc, yloc, !!marker)
        #select columns that are postive for marker
        dat2 = dat %>% dplyr::filter(get(marker) != 0)
        #calculate the observed K
        kobs = spatstat.geom::ppp(dat2$xloc, dat2$yloc, window = win) %>%
          spatstat.core::Kest(r = r_range, correction = edge_correction) %>%
          data.frame() %>%
          dplyr::rename("Theoretical K" = 2, "Observed K" = 3) %>%
          dplyr::mutate(Label = unique(spat[[mif$sample_id]]),
                        Marker = marker, .before=1)
        #calculate the permuted K
        kperms = parallel::mclapply(seq(num_permutations), function(perm){
          dat2 = dat[sample(seq(nrow(dat)), sum(dat[[marker]]), replace=F),]
          spatstat.geom::ppp(dat2$xloc, dat2$yloc, window = win) %>%
            spatstat.core::Kest(r = r_range, correction = edge_correction) %>%
            data.frame() %>%
            dplyr::rename("Theoretical K" = 2, "Permuted K" = 3) %>%
            dplyr::mutate(Label = unique(spat[[mif$sample_id]]),
                          Marker = marker,.before=1) %>%
            dplyr::mutate(iter = as.character(perm), .before = 1)
        }) %>%
          do.call(dplyr::bind_rows, .)
        K = dplyr::full_join(kobs, kperms,
                             by = c("Label", "Marker", "r", "Theoretical K"))
      }) %>%
        do.call(dplyr::bind_rows, .) %>%
        dplyr::mutate(`Degree of Clustering Permutation` = `Observed K` - `Permuted K`,
                      `Degree of Clustering Theoretical` = `Observed K` - `Theoretical K`)
    } else if(keep_perm_dis == F){
      res = parallel::mclapply(mnames, function(marker){
        #select the center of cells and marker column
        dat = spat %>%
          dplyr::select(xloc, yloc, !!marker)
        #select columns that are postive for marker
        dat2 = dat %>% dplyr::filter(get(marker) != 0)
        
        if(nrow(dat2) < 3){
          return(data.frame(iter = "Estimator", 
                            Label = spat[1,"deidentified_sample"],
                            Marker = marker,
                            r = r_range,
                            `Theoretical K` = pi * r_range^2,
                            `Observed K` = NA,
                            `Permuted K` = NA,
                            check.names = F))
        }
        #calculate the observed K
        kobs = spatstat.geom::ppp(dat2$xloc, dat2$yloc, window = win) %>%
          spatstat.core::Kest(r = r_range, correction = edge_correction) %>%
          data.frame() %>%
          dplyr::rename("Theoretical K" = 2, "Observed K" = 3) %>%
          dplyr::mutate(Label = unique(spat[[mif$sample_id]]),
                        Marker = marker, .before=1,
                        r = round(r))
        spat2 = spat %>%
          dplyr::select(1:3, !!marker) %>%
          dplyr::mutate(background = ifelse(get(marker) == 1, 0, 1)) %>%
          tidyr::gather("marks", "positive", -c(1:3)) %>%
          dplyr::filter(positive == 1 )
        pp_obj = spatstat.geom::ppp(x = spat2$xloc, y = spat2$yloc, window = win, marks = spat2$marks)
        k_est = get_kperm(pp_obj, mark1 = marker, r_vec = r_range, correction = edge_correction) %>%
          dplyr::rename("Permuted K" = 2) %>%
          dplyr::mutate(r = round(r))
        final = dplyr::full_join(kobs, k_est) %>%
          dplyr::mutate(iter = "Estimater", .before = 1)
        return(final)
      }) %>% #collapse all markers for spat
        do.call(dplyr::bind_rows, .) %>%
        dplyr::mutate(`Degree of Clustering Permutation` = `Observed K` - `Permuted K`,
                      `Degree of Clustering Theoretical` = `Observed K` - `Theoretical K`)
    }
    #return final results
    return(res)
  }, mc.cores = workers, mc.allow.recursive = T, mc.preschedule = F) %>% #set mclapply params
    do.call(dplyr::bind_rows, .) #collapse all samples to single data frame
  
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
