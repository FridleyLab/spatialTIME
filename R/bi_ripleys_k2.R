#' Bivariate Ripley's K
#'
#' @param spatial spatial data frame with columns of XMin, XMax, YMin, YMax
#' @param mnames vector of column names for phenotypes or data frame of marker combinations
#' @param r_range vector range of radii to calculate co-localization *K*
#' @param edge_correction character edge_correction method, one of "translation", "border", "or none" 
#' @param num_permutations integer number of permutations to estimate CSR
#' @param keep_permutation_distribution boolean as to whether to summarise permutations to mean
#' @param overwrite boolean as to whether to replace existing bivariate_Count if exists
#' @param workers integer number of CPU workers to use when number of `anchor` or `counted` is greater than 10,000
#' @param big integer used as the threshold for subsetting large samples, default is 1000 either *i* or *j*
#' @param nlarge number of cells in either *i* or *j* to flip to no edge correction - at small (relative to whole spatial region) *r* values differences in results between correction methods is negligible so running a few samples is recommended. Perhaps compute outweighs small differences in correction methods.
#'
#' @return mif object with bivariate ripley's K calculated
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
#' x2 = bi_ripleys_k2(mif = x, mnames = mnames_good, 
#'                    r_range = 0:100, edge_correction = "translation", 
#'                    num_permutations = 50, keep_permutation_distribution = FALSE, 
#'                    workers = 6, big = 1000)
bi_ripleys_k2 = function(mif,
                         mnames,
                         r_range = 0:100,
                         edge_correction = "translation",
                         num_permutations = 50,
                         keep_permutation_distribution = FALSE,
                         overwrite = TRUE,
                         workers = 6,
                         big = 1000,
                         nlarge = 1000){
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
  #split mif into jobs for spatial files
  #split mif into jobs for spatial files
  out = parallel::mclapply(names(mif$spatial), function(spatial_name){
    #prepare spatial data with x and y location (cell centers)
    spat = mif$spatial[[spatial_name]] %>%
      dplyr::mutate(xloc = (XMin + XMax)/2,
                    yloc = (YMin + YMax)/2)
    #find the window of the point process
    win = spatstat.geom::convexhull.xy(spat$xloc, spat$yloc)
    #calculate area of the window
    area = spatstat.geom::area(win)
    #matrix operations are WAY faster than data frame
    #since now all numeric, easy enough to use matrix
    spat = as.matrix(spat[,c("xloc", "yloc", mnames)])
    #get the combinations data frame
    if(inherits(mnames, "data.frame")){
      m_combos = mnames
    }
    if(inherits(mnames, "character")){
      m_combos = expand.grid(anchor = mnames,
                             counted = mnames)
    }
    #remove combinations that have 
    m_combos = m_combos %>%
      rowwise() %>%
      filter(!grepl(gsub(" \\(.*", "+", anchor), gsub(" \\(.*", "+", counted)) &
               !grepl(gsub(" \\(.*", "+", counted), gsub(" \\(.*", "+", anchor)) &
               !grepl(gsub("\\+", "\\\\+", anchor), gsub(" \\(.*", "+", counted)) &
               !grepl(gsub("\\+", "\\\\+", counted), gsub(" \\(.*", "+", anchor))) %>%
      ungroup()
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
                           `Theoretical K` = pi*r_range^2,
                           `Observed K` = NA,
                           `Permuted K` = NA,
                           check.names=FALSE) %>%
          dplyr::full_join(expand.grid(r = r_range,
                                       iter = seq(num_permutations)),
                           by = "r")
        return(final)
      }
      
      K_obs = data.frame(r = r_range,
                         `Theoretical K` = pi * r_range^2,
                         check.names = FALSE)
      #point pattern for anchor for distance matrix
      ppi = spatstat.geom::ppp(x = i_dat[,1],
                               y = i_dat[,2],
                               window = win)
      #point pattern for counted for distance matrix
      ppj = spatstat.geom::ppp(x = j_dat[,1],
                               y = j_dat[,2],
                               window = win)
      #ns
      li = nrow(i_dat)
      lj = nrow(j_dat)
      #find intensity of each marker
      lambdai = li/area
      lambdaj = lj/area
      
      #if the data are too large, distance matrix calculating will fail
      #split the matrix into a tiles to count cells within range
      if(li > big | lj > big){
        #find the number of tiles of distance matrix in columns and rows
        i_slides = floor(li / big)
        j_slides = floor(lj / big)
        #produce vectors of length ns, create T/F vector of which to subset
        i_ranges = getTile(slide = i_slides, l = li, size = big)
        j_ranges = getTile(slide = j_slides, l = lj, size = big)
        #count up those within the specified distance
        counts = parallel::mclapply(i_ranges, function(i_section){
          j_out = parallel::mclapply(j_ranges, function(j_section){
            #subset to tile
            i_tmp = ppi[i_section]
            j_tmp = ppj[j_section]
            
            #if the correction method is set to border, then run the num and den from spatstat Kmulti
            if(correction %in% c("border")){
              bI = spatstat.geom::bdist.points(i_tmp)
              bcloseI = bI[spatstat.geom::crosspairs(i_tmp, j_tmp, max(r_range), what = "ijd")$i]
              RS = spatstat.core::Kount(spatstat.geom::crosspairs(i_tmp, j_tmp, max(r_range), what = "ijd")$d,
                                        spatstat.geom::crosspairs(i_tmp, j_tmp, max(r_range), what = "ijd")$i,
                                        bI, spatstat.geom::handle.r.b.args(r_range, breaks=NULL, win, rmaxdefault = max(r_range)))
              return(RS)
            }
            #calculate distances
            dists = spatstat.geom::crossdist(i_tmp, j_tmp)
            rmv_i = rowSums(dists < max(r_range)) != 0
            rmv_j = colSums(dists < max(r_range)) != 0
            i_tmp = i_tmp[rmv_i]
            j_tmp = j_tmp[rmv_j]
            dists = spatstat.geom::crossdist(i_tmp, j_tmp)
            if(0 %in% dim(dists)){
              return(rep(0, length(r_range)))
            }
            #calculate edge correcion
            if(correction %in% c("trans", "translation")){
              edge = spatstat.core::edge.Trans(i_tmp, j_tmp)
              #count edge correction matrix for cells within range r in distance matrix
              counts = sapply(r_range, function(r){sum(edge[which(dists < r)])})
              #remove large distance and edge correction matrix to keep ram usage down
              rm(dists, edge)
              #return counts for tile
              return(counts)
            }
            if(correction %in% c("none")){
              counts = cumsum(spatstat.geom::whist(spatstat.geom::crosspairs(i_tmp, j_tmp, max(r_range), what = "ijd")$d,
                                                   spatstat.geom::handle.r.b.args(r_range, breaks=NULL, win, rmaxdefault = max(r_range))$val))
              return(counts)
            }
          }) 
          
          if(correction == "border"){
            num = lapply(j_out, function(j_big){
              j_big[[1]]
            }) %>%
              do.call(rbind.data.frame, .) %>%
              colSums() %>%
              unname()
            den = j_out[[1]][[2]]
            return(list(num = num, den = den))
          }
          j_out  %>% #use 1 core per tile
            #bind all j tiles to data frame
            do.call(rbind.data.frame, .) %>%
            #take column sums and return
            colSums()
        }, mc.preschedule = F, mc.allow.recursive = T)
        
        if(correction == "border"){
          num = lapply(counts, function(j_big){
            j_big[[1]]
          }) %>%
            do.call(rbind.data.frame, .) %>%
            colSums() %>%
            unname()
          den = ppi$n
          K_obs$`Observed K` = num / (lambdaj * den)
        } else {
          counts = counts %>%
            #bind i tile counts and sum
            do.call(rbind.data.frame, .) %>%
            colSums() %>%
            unname()
          #calulate clustering from counted edge corrections with intensities of i and j
          K_obs$`Observed K` = (1/(lambdai * lambdaj * area)) * counts
        }
        
        
      }
      #if there are less than 10,000 j or i just compute matrix
      if(!(li > big | lj > big)){
        #calculate distance matrix
        sp_tmp2 = spat_tmp %>%
          data.frame(check.names=F)%>%
          tidyr::gather(Marker, Positive, -xloc, -yloc) %>%
          dplyr::filter(Positive == 1)
        pp_obj = spatstat.geom::ppp(x = sp_tmp2$xloc, y = sp_tmp2$yloc, window = win, marks = factor(sp_tmp2$Marker))
        K_obs = spatstat.core::Kcross(pp_obj,
                                      i = unique(sp_tmp2$Marker)[1],
                                      j = unique(sp_tmp2$Marker)[2],
                                      r = r_range,
                                      correction = correction) %>%
          data.frame() %>%
          dplyr::rename("Theoretical K" = 2,
                        "Observed K" = 3) %>%
          dplyr::full_join(K_obs, ., by = c("r", "Theoretical K"))
      }
      
      K_obs$Anchor = anchor
      K_obs$Counted = counted
      #randomly sample the rows of possible cell locations for permuting
      perm_rows = lapply(seq(num_permutations), function(x){
        sample(1:nrow(spat), sum(nrow(i_dat), nrow(j_dat)), replace = FALSE)
      })
      assign("perm_rows", perm_rows, envir = .GlobalEnv)
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
                            `Theoretical K` = pi * r_range^2,
                            iter = perm_n,
                            check.names = FALSE)
        #point pattern for anchor for distance matrix
        ppi = spatstat.geom::ppp(x = i_dat[,1],
                                 y = i_dat[,2],
                                 window = win)
        #point pattern for counted for distance matrix
        ppj = spatstat.geom::ppp(x = j_dat[,1],
                                 y = j_dat[,2],
                                 window = win)
        #ns
        li = nrow(i_dat)
        lj = nrow(j_dat)
        #find intensity of each marker
        lambdai = li/area
        lambdaj = lj/area
        
        #if the data are too large, distance matrix calculating will fail
        #split the matrix into a tiles to count cells within range
        if(li > big | lj > big){
          #find the number of tiles of distance matrix in columns and rows
          i_slides = floor(li / big)
          j_slides = floor(lj / big)
          #produce vectors of length ns, create T/F vector of which to subset
          i_ranges = getTile(slide = i_slides, l = li, size = big)
          j_ranges = getTile(slide = j_slides, l = lj, size = big)
          #count up those within the specified distance
          counts = parallel::mclapply(i_ranges, function(i_section){
            j_out = parallel::mclapply(j_ranges, function(j_section){
              #subset to tile
              i_tmp = ppi[i_section]
              j_tmp = ppj[j_section]
              
              if(correction %in% c("border")){
                bI = spatstat.geom::bdist.points(i_tmp)
                bcloseI = bI[spatstat.geom::crosspairs(i_tmp, j_tmp, max(r_range), what = "ijd")$i]
                RS = spatstat.core::Kount(spatstat.geom::crosspairs(i_tmp, j_tmp, max(r_range), what = "ijd")$d,
                                          spatstat.geom::crosspairs(i_tmp, j_tmp, max(r_range), what = "ijd")$i,
                                          bI, spatstat.geom::handle.r.b.args(r_range, breaks=NULL, win, rmaxdefault = max(r_range)))
                return(RS)
              }
              #calculate distances
              dists = spatstat.geom::crossdist(i_tmp, j_tmp)
              rmv_i = rowSums(dists < max(r_range)) != 0
              rmv_j = colSums(dists < max(r_range)) != 0
              i_tmp = i_tmp[rmv_i]
              j_tmp = j_tmp[rmv_j]
              dists = spatstat.geom::crossdist(i_tmp, j_tmp)
              if(0 %in% dim(dists)){
                return(rep(0, length(r_range)))
              }
              #calculate edge correcion
              if(correction %in% c("trans", "translation")){
                edge = spatstat.core::edge.Trans(i_tmp, j_tmp)
                #count edge correction matrix for cells within range r in distance matrix
                counts = sapply(r_range, function(r){sum(edge[which(dists < r)])})
                #remove large distance and edge correction matrix to keep ram usage down
                rm(dists, edge)
                #return counts for tile
                return(counts)
              }
              if(correction %in% c("none")){
                counts = cumsum(spatstat.geom::whist(spatstat.geom::crosspairs(i_tmp, j_tmp, max(r_range), what = "ijd")$d,
                                                     spatstat.geom::handle.r.b.args(r_range, breaks=NULL, win, rmaxdefault = max(r_range))$val))
                return(counts)
              }
            }) 
            
            if(correction == "border"){
              num = lapply(j_out, function(j_big){
                j_big[[1]]
              }) %>%
                do.call(rbind.data.frame, .) %>%
                colSums() %>%
                unname()
              den = j_out[[1]][[2]]
              return(list(num = num, den = den))
            }
            j_out  %>% #use 1 core per tile
              #bind all j tiles to data frame
              do.call(rbind.data.frame, .) %>%
              #take column sums and return
              colSums()
          }, mc.preschedule = F, mc.allow.recursive = T)
          
          if(correction == "border"){
            num = lapply(counts, function(j_big){
              j_big[[1]]
            }) %>%
              do.call(rbind.data.frame, .) %>%
              colSums() %>%
              unname()
            den = ppi$n
            K_obs$`Observed K` = num / (lambdaj * den)
          } else {
            counts = counts %>%
              #bind i tile counts and sum
              do.call(rbind.data.frame, .) %>%
              colSums() %>%
              unname()
            #calulate clustering from counted edge corrections with intensities of i and j
            permed$`Permuted K` = (1/(lambdai * lambdaj * area)) * counts
          }
        }
        if(!(li > big | lj > big)){
          #calculate distance matrix
          dat2 = data.frame(dat)
          dat2$Marker = label
          pp_obj = spatstat.geom::ppp(x = dat2$xloc, y = dat2$yloc, window = win, marks = factor(dat2$Marker))
          permed = spatstat.core::Kcross(pp_obj,
                                         i = unique(dat2$Marker)[1],
                                         j = unique(dat2$Marker)[2],
                                         r = r_range,
                                         correction = correction) %>%
            data.frame() %>%
            dplyr::rename("Theoretical K" = 2,
                          "Permuted K" = 3) %>%
            dplyr::full_join(permed, ., by = c("r", "Theoretical K"))
        }
        return(permed)
      }, mc.preschedule = F, mc.allow.recursive = T) %>%
        do.call(dplyr::bind_rows, .)
      #join the emperical K and the permuted CSR estimate
      final = dplyr::full_join(K_obs,
                               kpermed, by = c("r", "Theoretical K")) %>%
        #add the image label to the data frame
        dplyr::mutate(Label = spatial_name, .before = 1)
      
      return(final)
    }) %>% #, mc.cores = cores, mc.preschedule = F,mc.allow.recursive = T
      do.call(dplyr::bind_rows, .)
    #reorder columns to make more sense
    res = res[,c(1,2,7,5,6,3,4,8)]
    return(res)
  }, mc.cores = workers, mc.preschedule = F,mc.allow.recursive = T) %>%
    do.call(dplyr::bind_rows, .)%>% #collapse all samples to single data frame
    rename(!!mif$sample_id := Label)
  #if user doesn't want the permutation distribution, get average of the permutation estimate
  if(!keep_permutation_distribution){
    out = out %>%
      #remove iter since this is the permutation number
      dplyr::select(-iter) %>%
      #group by those used for permuting
      dplyr::group_by(across(mif$sample_id), r, Anchor, Counted) %>%
      #take mean of theoretical, permuted, observed
      dplyr::summarise_all(~mean(., na.rm=TRUE)) %>%
      #calculate the degree of clustering from both the theoretical and permuted
      dplyr::mutate(`Degree of Clustering Permutation` = `Observed K` - `Permuted K`,
                    `Degree of Clustering Theoretical` = `Observed K` - `Theoretical K`)
  }
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

#get the tile vectors of T/F for subsetting the ppp
getTile = function(slide, l, size){
  #apply seq slide count
  lapply(0:slide, function(s){
    #for first tile, get values from 1 to either size or length of positives
    if(s == 0){
      w = 1:ifelse(size<l, size, l)
    }
    #between first and last tile, increment in size increments
    if(s > 0 & s != slide){
      w = (s*size+1):((s+1)*size)
    }
    #for last tile, go from previous max to length of positives
    if(s == slide){
      w = (s*size+1):(l)
    }
    #convert the tile to T/F vector for subsetting
    v = rep(FALSE, l)
    v[w] = TRUE
    v
  })
}
