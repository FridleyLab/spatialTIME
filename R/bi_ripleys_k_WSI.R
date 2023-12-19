#' Bivariate Ripley's K for Whole Slide Images
#'
#' @param mif mIF object with spatial data frames, clinical, and per-sample summary information
#' @param mnames vector of column names for phenotypes or data frame of marker combinations
#' @param r_range vector range of radii to calculate co-localization *K*
#' @param edge_correction character edge_correction method, one of "translation", or none" 
#' @param num_permutations integer number of permutations to estimate CSR
#' @param permute whether or not to use permutations to estimate CSR (TRUE) or to calculate exact CSR (FALSE)
#' @param keep_permutation_distribution boolean as to whether to summarise permutations to mean
#' @param overwrite boolean as to whether to replace existing bivariate_Count if exists
#' @param workers integer number of CPU workers to use
#' @param big integer used as the threshold for subsetting large samples, default is 1000 either *i* or *j*
#' @param nlarge number of cells in either *i* or *j* to flip to no edge correction - at small (relative to whole spatial region) *r* values differences in results between correction methods is negligible so running a few samples is recommended. Perhaps compute outweighs small differences in correction methods.
#' @param xloc the x and y positions that correspond to cells. If left as NULL, XMin, XMax, YMin, and YMax must be present in the spatial files
#' @param yloc the x and y positions that correspond to cells. If left as NULL, XMin, XMax, YMin, and YMax must be present in the spatial files
#'
#' @return mif object with bivariate Ripley's K calculated
#' 
#' @description
#' Bivariate Ripley's K function within spatialTIME, `bi_ripleys_k` is a function that takes in a `mIF` object, along with 
#' some parameters like marker names of interest and range of radii in which to assess bivariate clustering or colocalization.
#' In 1.3.3.3 we have introduced the ability to forsgo the need for permutations with the implementation of the exact CSR estimate.
#' This is both faster and being the exact CSR, produces an exact degree of clustering in the spatial files.
#' 
#' Due to the availability of whole slide images (WSI), there's a possibility users will be running bivariate Ripley's K on samples
#' that have millions of cells. When doing this, keep in mind that a nearest neighbor matrix with *n* cell is *n* by *n* in size and 
#' therefore easily consumers high performance compute levels of RAM. To combat this, we have implemented a tiling method that performs
#' counts for small chunks of the distance matrix at a time before finally calculating the bivariate Ripley's K value on the total counts.
#' When doing this there are now 2 import parameters to keep in mind. The `big` parameter is the size of the tile to use. We have found
#' 1000 to be a good number that allows for high number of cores while maintaining low RAM usage. The other important parameter when
#' working with WSI is nlarge which is the fall over for switching to no edge correction. The spatstat.explore::Kest univariate 
#' Ripley's K uses a default of 3000 but we have defaulted to 1000 to keep compute minimized as edge correction uses large amounts
#' of RAM over 'none'.
#' 
#' @export
#'
#' @examples
#' x <- spatialTIME::create_mif(clinical_data = spatialTIME::example_clinical %>% 
#'                                dplyr::mutate(deidentified_id = as.character(deidentified_id)),
#'                              sample_data = spatialTIME::example_summary %>% 
#'                                dplyr::mutate(deidentified_id = as.character(deidentified_id)),
#'                              spatial_list = spatialTIME::example_spatial,
#'                              patient_id = "deidentified_id", 
#'                              sample_id = "deidentified_sample")
#' mnames_good <- c("CD3..Opal.570..Positive","CD8..Opal.520..Positive",
#'                  "FOXP3..Opal.620..Positive","PDL1..Opal.540..Positive",
#'                  "PD1..Opal.650..Positive","CD3..CD8.","CD3..FOXP3.")
#' x2 = bi_ripleys_k_WSI(mif = x, mnames = mnames_good[1:2], 
#'                    r_range = 0:100, edge_correction = "none", permute = FALSE,
#'                    num_permutations = 50, keep_permutation_distribution = FALSE, 
#'                    workers = 1, big = 1000)
bi_ripleys_k_WSI = function(mif,
                         mnames,
                         r_range = 0:100,
                         edge_correction = "translation",
                         num_permutations = 50,
                         permute = FALSE, #redo for permutation or estimate
                         keep_permutation_distribution = FALSE,
                         overwrite = TRUE,
                         workers = 6,
                         big = 1000,
                         nlarge = 1000,
                         xloc = NULL,
                         yloc = NULL){
  Label = Anchor = Counted = `Exact CSR` = NULL
  #check whether the object assigned to mif is of class mif
  if(!inherits(mif, "mif")){
    stop("Please use a mIF object for mif")
  }
  #check whether mnames is either a character vector a data frame
  if(!inherits(mnames, "character") & !inherits(mnames, "data.frame")){
    stop("Please use either a character vector or data frame of marker combinations for mnames")
  }
  #r_range has to have 0 for use with AUC (0,0)
  if(!(0 %in% r_range)){
    r_range = c(0, r_range)
  }
  if(!(edge_correction %in% c("trans", "none", "translation")))
    stop("provide either translation or none for border correction")
  #split mif into jobs for spatial files
  #split mif into jobs for spatial files
  out = parallel::mclapply(names(mif$spatial), function(spatial_name){
    #prepare spatial data with x and y location (cell centers)
    spat = mif$spatial[[spatial_name]]
    
    if(is.null(xloc) & is.null(yloc)){
      spat = spat %>%
        dplyr::mutate(xloc = (XMin + XMax)/2,
                      yloc = (YMin + YMax)/2)
    } else {
      spat = spat %>%
        dplyr::rename('xloc' := xloc,
                      'yloc' := yloc)
    }
    #find the window of the point process
    win = spatstat.geom::convexhull.xy(spat$xloc, spat$yloc)
    #calculate area of the window
    area = spatstat.geom::area(win)
    #matrix operations are WAY faster than data frame
    #since now all numeric, easy enough to use matrix
    spat = as.matrix(spat[,c("xloc", "yloc", as.character(unique(unlist(mnames))))])
    #get the combinations data frame
    if(inherits(mnames, "data.frame")){
      m_combos = mnames
    }
    if(inherits(mnames, "character")){
      m_combos = expand.grid(anchor = mnames,
                             counted = mnames) %>%
        dplyr::filter(anchor != counted)
    }
    #calculating exact K works now!
    if(!permute){
      #calculate exact K
      exact_K = spatstat.explore::Kest(spatstat.geom::ppp(x = spat[,'xloc'],
                                                          y = spat[,'yloc'],
                                                          window = win),
                                       r = r_range,
                                       correction = edge_correction) %>%
        data.frame() %>%
        dplyr::pull(3)
    }
    
    #for the combinations of markers, do bivark and permutations
    res = parallel::mclapply(1:nrow(m_combos), function(combo){
      #pull anchor and counted marker from combos data frame
      anchor = m_combos[combo, ] %>% dplyr::pull(anchor) %>% as.character()
      counted = m_combos[combo, ] %>% dplyr::pull(counted) %>% as.character()
      cat(spatial_name, "\t", combo, "\t", anchor, "\t", counted, "\n")
      #remove rows that are positive for both counted and anchor
      spat_tmp = spat[!(spat[,anchor] == 1 & spat[,counted] == 1),c("xloc", "yloc", anchor, counted)]
      
      #get data for both anchor and counted
      i_dat = spat_tmp[spat_tmp[,3] == 1,]
      j_dat = spat_tmp[spat_tmp[,4] == 1,]
      
      #if the number of positive cells for either counted or anchor is less than 2, return empty K
      if(sum(spat_tmp[,3]) < 2 | sum(spat_tmp[,4]) < 2){
        final = data.frame(Label = spatial_name,
                           r = r_range,
                           Anchor = anchor,
                           Counted = counted,
                           `Theoretical CSR` = pi*r_range^2,
                           `Observed K` = NA,
                           `Permuted CSR` = NA,
                           `Exact CSR` = NA,
                           check.names=FALSE) 
        if(permute){
          final = final %>%
          dplyr::full_join(expand.grid(r = r_range,
                                       iter = seq(num_permutations)),
                           by = "r")
        } else {
          final$iter = 1
        }
        return(final)
      }
      #make empty data frame to begin the final K table
      K_obs = data.frame(r = r_range,
                         `Theoretical CSR` = pi * r_range^2,
                         check.names = FALSE)
      #get observed K
      K_obs$`Observed K` = calculateK(i_dat = i_dat, 
                                      j_dat = j_dat,
                                      anchor = anchor,
                                      counted = counted, 
                                      area = area, win = win, big = big, 
                                      r_range = r_range,
                                      edge_correction = edge_correction,
                                      cores = 1)
      #set the anchor and counted in final table
      K_obs$Anchor = anchor
      K_obs$Counted = counted
      
      if(permute){
        #randomly sample the rows of possible cell locations for permuting
        perm_rows = lapply(seq(num_permutations), function(x){
          sample(1:nrow(spat), sum(nrow(i_dat), nrow(j_dat)), replace = FALSE)
        })
        
        #assign("perm_rows", perm_rows, envir = .GlobalEnv)
        #calculate BiK for each permutation of cells
        kpermed = parallel::mclapply(seq(perm_rows), function(perm_n){
          cat(perm_n)
          #extract vector of rows for permutation run
          perm = perm_rows[[perm_n]]
          #subset the x and y coords for those to use for permutation
          dat = spat[perm,1:2]
          #create label vector of anchor and counted cells of length anchor + counted
          label = c(rep(anchor, nrow(i_dat)), rep(counted, nrow(j_dat)))
          #subset permute rows for anchor and counted
          i_dat = dat[label == anchor,]
          j_dat = dat[label == counted,]
          #prep permtued K table
          permed = data.frame(r = r_range,
                              `Theoretical CSR` = pi * r_range^2,
                              iter = perm_n,
                              check.names = FALSE)
          permed$`Permuted CSR` = calculateK(i_dat = i_dat, 
                                           j_dat = j_dat,
                                           anchor = anchor,
                                           counted = counted, 
                                           area = area, win = win, big = big, 
                                           r_range = r_range,
                                           edge_correction = edge_correction,
                                           cores = 1)
          return(permed)
        }, mc.preschedule = FALSE, mc.allow.recursive = TRUE) %>%
          do.call(dplyr::bind_rows, .)
        kpermed$`Exact CSR` = NA
      } else {
        kpermed = data.frame(r = r_range,
                             `Theoretical CSR` = pi * r_range^2,
                             iter = 1,
                             check.names = FALSE)
        kpermed$`Permuted CSR` = NA
        kpermed$`Exact CSR` = exact_K
      }
      
      
      #join the emperical K and the permuted CSR estimate
      final = dplyr::full_join(K_obs,
                               kpermed, by = c("r", "Theoretical CSR")) %>%
        #add the image label to the data frame
        dplyr::mutate(Label = spatial_name, .before = 1)
      
      return(final)
    }) %>% #, mc.cores = cores, mc.preschedule = F,mc.allow.recursive = T
      do.call(dplyr::bind_rows, .)
    #reorder columns to make more sense
    res = res[,c(1,2,7,5,6,3,4,8,9)]
    return(res)
  }, mc.cores = workers, mc.preschedule = FALSE,mc.allow.recursive = TRUE) %>%
    do.call(dplyr::bind_rows, .)%>% #collapse all samples to single data frame
    dplyr::rename(!!mif$sample_id := Label)
  #if user doesn't want the permutation distribution, get average of the permutation estimate
  if(!keep_permutation_distribution & permute){
    out = out %>%
      #remove iter since this is the permutation number
      dplyr::select(-iter) %>%
      #group by those used for permuting
      dplyr::group_by(dplyr::across(mif$sample_id), r, Anchor, Counted) %>%
      #take mean of theoretical, permuted, observed
      dplyr::summarise_all(~mean(., na.rm=TRUE)) %>%
      #calculate the degree of clustering from both the theoretical and permuted
      dplyr::mutate(`Degree of Clustering Permutation` = `Observed K` - `Permuted CSR`,
                    `Degree of Clustering Theoretical` = `Observed K` - `Theoretical CSR`,
                    `Exact CSR` = NA) %>%
      dplyr::mutate(iter = num_permutations, .before = Anchor)
  }
  #if overwrite is true, replace the bivariate count in the derived slot
  if(overwrite){
    mif$derived$bivariate_Count = out %>%
      #calculate the degree of clustering from both the theoretical and permuted
      dplyr::mutate(`Degree of Clustering Theoretical` = `Observed K` - `Theoretical CSR`,
                    `Degree of Clustering Permutation` = `Observed K` - `Permuted CSR`,
                    `Degree of Clustering Exact` = `Observed K` - `Exact CSR`) %>%
      #add run number to differentiate between bivariate compute runs
      dplyr::mutate(Run = 1)
  }
  #if don't overwrite
  if(!overwrite){
    #bind old and new bivar runs together, incrementing Run
    mif$derived$bivariate_Count = mif$derived$bivariate_Count%>%
      dplyr::bind_rows(out %>%
                         #calculate the degree of clustering from both the theoretical and permuted
                         dplyr::mutate(`Degree of Clustering Theoretical` = `Observed K` - `Theoretical CSR`,
                                       `Degree of Clustering Permutation` = `Observed K` - `Permuted CSR`,
                                       `Degree of Clustering Exact` = `Observed K` - `Exact CSR`) %>%
                         dplyr::mutate(Run = ifelse(exists("bivariate_Count", mif$derived),
                                                    max(mif$derived$bivariate_Count$Run) + 1,
                                                    1)))
  }
  #return the final mif object
  return(mif)
}
