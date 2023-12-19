#' Bivariate Ripley's K
#'
#' @param mif mIF object with spatial data frames, clinical, and per-sample summary information
#' @param mnames vector of column names for phenotypes or data frame of marker combinations
#' @param r_range vector range of radii to calculate co-localization *K*
#' @param edge_correction character edge_correction method, one of "translation", "border", "or none" 
#' @param num_permutations integer number of permutations to estimate CSR
#' @param permute whether or not to use permutations to estimate CSR (TRUE) or to calculate exact CSR (FALSE)
#' @param keep_permutation_distribution boolean as to whether to summarise permutations to mean
#' @param overwrite boolean as to whether to replace existing bivariate_Count if exists
#' @param workers integer number of CPU workers to use
#' @param xloc,yloc the x and y positions that correspond to cells. If left as NULL, XMin, XMax, YMin, and YMax must be present in the spatial files
#' @param force logical whether or not to continue if sample has more than 10,000 cells
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
#' x2 = bi_ripleys_k(mif = x, mnames = mnames_good[1:2], 
#'                    r_range = 0:100, edge_correction = "none", permute = FALSE,
#'                    num_permutations = 50, keep_permutation_distribution = FALSE, 
#'                    workers = 1)
bi_ripleys_k = function(mif,
                         mnames,
                         r_range = 0:100,
                         edge_correction = "translation",
                         num_permutations = 50,
                         permute = FALSE, #redo for permutation or estimate
                         keep_permutation_distribution = FALSE,
                         overwrite = TRUE,
                         workers = 6,
                         xloc = NULL,
                         yloc = NULL, force = FALSE){
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
  #check if user should be using wsi method
  if(any(sapply(mif$spatial, nrow) > 10000)){
    if(!force){
      stop("Some samples have a large number of cells - bi_ripleys_k_WSI may be more appropriate.")
    } else {
      message("Some samples have a large number of cells - bi_ripleys_k_WSI may be more appropriate.\nContinuing\n")
    }
  }
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
    core_pp = spatstat.geom::ppp(x = spat[,'xloc'],
                                 y = spat[,'yloc'],
                                 window = win)
    #calculating exact K works now!
    if(!permute){
      #calculate exact K
      exact_K = spatstat.explore::Kest(core_pp,
                                       r = r_range,
                                       correction = edge_correction) %>%
        data.frame() %>%
        dplyr::rename("Theoretical CSR" = 2,
                      "Exact CSR" = 3)
    }
    
    #for the combinations of markers, do bivark and permutations
    res = parallel::mclapply(1:nrow(m_combos), function(combo){
      #pull anchor and counted marker from combos data frame
      anchor = m_combos[combo, ] %>% dplyr::pull(anchor) %>% as.character()
      counted = m_combos[combo, ] %>% dplyr::pull(counted) %>% as.character()
      cat(spatial_name, "\t", combo, "\t", anchor, "\t", counted, "\n")
      #remove rows that are positive for both counted and anchor
      spat_tmp = get_bi_rows(data.frame(spat, check.names = FALSE), c(anchor, counted))
      tabs = table(spat_tmp$Marker)
      #if the number of positive cells for either counted or anchor is less than 2, return empty K
      if(length(tabs) < 2 | any(tabs < 2)){
        final = data.frame(Label = spatial_name,
                           Anchor = anchor,
                           Counted = counted,
                           r = r_range,
                           `Theoretical CSR` = pi*r_range^2,
                           `Permuted CSR` = NA,
                           `Exact CSR` = NA,
                           `Observed K` = NA,
                           check.names=FALSE) 
        if(permute & keep_permutation_distribution){
          final = final %>%
            dplyr::full_join(expand.grid(r = r_range,
                                         iter = seq(num_permutations)),
                             by = "r")
        } else if(permute & !keep_permutation_distribution){
          final$iter = num_permutations
        } else {
          final$iter = 1
        }
        return(final %>%
                 dplyr::relocate(iter, .after = 3))
      }
      #make empty data frame to begin the final K table
      ps = core_pp[spat_tmp$cell]
      spatstat.geom::marks(ps) = spat_tmp$Marker
      K_obs = spatstat.explore::Kcross(ps,
                                       i = anchor, j = counted,
                                       r = r_range,
                                       correction = edge_correction) %>%
        data.frame() %>%
        dplyr::rename("Theoretical CSR" = 2,
                      "Observed K" = 3)
      #set the anchor and counted in final table
      K_obs$Anchor = anchor
      K_obs$Counted = counted
      
      if(permute){
        #randomly sample the rows of possible cell locations for permuting
        perm_rows = lapply(seq(num_permutations), function(x){
          sample(1:nrow(spat), sum(tabs), replace = FALSE)
        })
        
        #assign("perm_rows", perm_rows, envir = .GlobalEnv)
        #calculate BiK for each permutation of cells
        kpermed = parallel::mclapply(seq(perm_rows), function(perm_n){
          cat(perm_n)
          #extract vector of rows for permutation run
          perm = perm_rows[[perm_n]]
          #subset the spatstat object
          ps = core_pp[perm]
          spatstat.geom::marks(ps) = spat_tmp$Marker
          
          permed = spatstat.explore::Kcross(ps,
                                            i = anchor, j = counted,
                                            r = r_range,
                                            correction = edge_correction) %>%
            data.frame() %>%
            dplyr::rename("Theoretical CSR" = 2,
                          "Permuted CSR" = 3) %>%
            dplyr::mutate(iter = perm_n)
          return(permed)
        }, mc.preschedule = FALSE, mc.allow.recursive = TRUE) %>%
          do.call(dplyr::bind_rows, .)
        # #simplify if not keeping perms
        # if(!keep_permutation_distribution){
        #   kpermed = kpermed %>%
        #     dplyr::group_by(r) %>%
        #     dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE))) %>%
        #     dplyr::mutate(iter = 1)
        # }
        kpermed$`Exact CSR` = NA
      } else {
        kpermed = data.frame(r = r_range,
                             `Theoretical CSR` = pi * r_range^2,
                             iter = 1,
                             check.names = FALSE)
        kpermed$`Permuted CSR` = NA
        kpermed$`Exact CSR` = exact_K$`Exact CSR`
      }
      
      
      #join the emperical K and the permuted CSR estimate
      final = dplyr::full_join(K_obs,
                               kpermed, by = c("r", "Theoretical CSR")) %>%
        #add the image label to the data frame
        dplyr::mutate(Label = spatial_name, .before = 1) %>% 
        dplyr::relocate(Anchor, Counted, iter, r, 
                        `Theoretical CSR`, `Permuted CSR`, `Exact CSR`, .after = 1) %>%
        dplyr::group_by(Label, Anchor, Counted, r) %>%
        dplyr::mutate(`Permutations Larger than Observed` = ifelse(permute,
                                                                   sum(`Permuted CSR` > `Observed K`, na.rm = TRUE),
                                                                   NA))
      if(permute & !keep_permutation_distribution){
        final = final %>%
          dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE))) %>%
          dplyr::mutate(iter = num_permutations) %>%
          dplyr::relocate(iter, .after = 3)
      }
        
      
      return(final)
    }, mc.preschedule = F,mc.allow.recursive = T) %>% #
      do.call(dplyr::bind_rows, .)
    #reorder columns to make more sense
    return(res)
  }, mc.cores = workers, mc.preschedule = FALSE,mc.allow.recursive = TRUE) %>%
    do.call(dplyr::bind_rows, .)%>% #collapse all samples to single data frame
    dplyr::rename(!!mif$sample_id := Label)
  out = out %>%
    #calculate the degree of clustering from both the theoretical and permuted
    dplyr::mutate(`Degree of Clustering Theoretical` = `Observed K` - `Theoretical CSR`,
                  `Degree of Clustering Permutation` = `Observed K` - `Permuted CSR`,
                  `Degree of Clustering Exact` = `Observed K` - `Exact CSR`)
  #if overwrite is true, replace the bivariate count in the derived slot
  if(overwrite){
    mif$derived$bivariate_Count = out %>%
      #add run number to differentiate between bivariate compute runs
      dplyr::mutate(Run = 1)
  }
  #if don't overwrite
  if(!overwrite){
    #bind old and new bivar runs together, incrementing Run
    mif$derived$bivariate_Count = mif$derived$bivariate_Count%>%
      dplyr::bind_rows(out %>%
                         dplyr::mutate(Run = ifelse(exists("bivariate_Count", mif$derived),
                                                    max(mif$derived$bivariate_Count$Run) + 1,
                                                    1)))
  }
  #return the final mif object
  return(mif)
}
