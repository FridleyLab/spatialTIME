#' Title
#'
#' @param spatial spatial data frame with columns of XMin, XMax, YMin, YMax
#' @param mnames vector of column names for phenotypes or data frame of marker combinations
#' @param r_range range of radii to calculate co-localization K
#' @param correction correction method, either "translation" or "none" currently
#' @param cores number of CPU cores to use when number of `anchor` or `counted` is greater than 10,000
#'
#' @return
#' @export
#'
#' @examples
#' 
bi_ripleys_k2 = function(mif,
                         mnames,
                         r_range = 0:100,
                         correction = "none",
                         num_permutations = 100,
                         keep_permutation_distribution = FALSE,
                         overwrite = TRUE,
                         cores = 6){
  
  if(!inherits(mif, "mif")){
    stop("Please use a mIF object for mif")
  }
  if(!inherits(mnames, "character") & !inherits(mnames, "data.frame")){
    stop("Please use either a character vector or data frame of marker combinations for mnames")
  }
  if(!(0 %in% r_range)){
    r_range = c(0, r_range)
  }
  
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
    m_combos = m_combos[m_combos$anchor != m_combos$counted,]
    #for the combinations of markers, do bivark and permutations
    res = parallel::mclapply(1:nrow(m_combos), function(combo){
      #pull anchor and counted marker from combos data frame
      anchor = as.character(m_combos[combo, 1])
      counted = as.character(m_combos[combo, 2])
      cat(anchor, "\t", counted, "\n")
      #remove rows that are positive for both counted and anchor
      spat_tmp = spat[!(spat[,anchor] == 1 & spat[,counted] == 1),c("xloc", "yloc", anchor, counted)]
      
      #get data for both anchor and counted
      i_dat = spat_tmp[spat_tmp[,3] == 1,]
      j_dat = spat_tmp[spat_tmp[,4] == 1,]
      
      #if the number of positive cells for either counted or anchor is less than 2, return empty K
      if(sum(spat_tmp[,3]) < 2 | sum(spat_tmp[,4]) < 2){
        final = data.frame(r = r_range,
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
      K_obs$`Observed K` = getBiK(i_dat, j_dat, area, r_range, win, correction)
      K_obs$Anchor = anchor
      K_obs$Counted = counted
      
      perm_rows = lapply(seq(num_permutations), function(x){
        sample(1:nrow(spat), sum(nrow(i_dat), nrow(j_dat)), replace = FALSE)
      })
      
      kpermed = parallel::mclapply(seq(perm_rows), function(perm_n){
        perm = perm_rows[[perm_n]]
        dat = spat[perm,1:2]
        label = c(rep(anchor, nrow(i_dat)), rep(counted, nrow(j_dat)))
        i_dat = dat[label == anchor,]
        j_dat = dat[label == counted,]
        permed = data.frame(r = r_range,
                            `Theoretical K` = pi * r_range^2,
                            iter = perm_n,
                            check.names = FALSE)
        permed$`Permuted K` = getBiK(i_dat, j_dat, area, r_range, win, correction)
        return(permed)
      }) %>%
        do.call(dplyr::bind_rows, .)
      
      final = dplyr::full_join(K_obs,
                               kpermed, by = c("r", "Theoretical K"))
      
      return(final)
    }) %>%
      do.call(dplyr::bind_rows, .) %>%
      dplyr::mutate(Label = spatial_name, .before = 1)
    res = res[,c(1,2,7,5,6,3,4,8)]
    return(res)
  }, mc.cores = cores, mc.preschedule = FALSE) %>%
    do.call(dplyr::bind_rows, .)
  
  if(!keep_permutation_distribution){
    out = out %>%
      dplyr::select(-iter) %>%
      dplyr::group_by(Label, r, Anchor, Counted) %>%
      dplyr::summarise_all(~mean(., na.rm=TRUE)) %>%
      dplyr::mutate(`Degree of Clustering Permutation` = `Observed K` - `Permuted K`,
                    `Degree of Clustering Theoretical` = `Observed K` - `Theoretical K`)
  }
  if(overwrite){
    mif$derived$bivariate_Count = out %>%
      dplyr::mutate(Run = 1)
  }
  if(!overwrite){
    mif$derived$bivariate_Count = mif$derived$bivariate_Count%>%
      dplyr::bind_rows(out %>%
                         dplyr::mutate(Run = ifelse(exists("bivariate_Count", mif$derived),
                                                    max(mif$derived$bivariate_Count$Run) + 1,
                                                    1)))
  }
  
  return(mif)
}

#get tile counts
getTile = function(slide, l, size){
  lapply(1:slide, function(s){
    if(s == 1){
      w = 1:size
    }
    if(s != 1){
      w = (s*size):((s+1)*size-1)
    }
    if(s == i_slides){
      w = (s*size):(l)
    }
    v = rep(FALSE, l)
    v[w] = TRUE
    v
  })
}

getBiK = function(i_dat, j_dat, area, r_range, win, correction){
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
  if(li > 10000 | lj > 10000){
    #find the number of tiles of distance matrix in columns and rows
    i_slides = floor(li / 10000)
    j_slides = floor(lj / 10000)
    #produce vectors of length ns, create T/F vector of which to subset
    i_ranges = getTile(slide = i_slides, l = li, size = 1000)
    j_ranges = getTile(slide = j_slides, l = lj, size = 1000)
    #count up those within the specified distance
    is = parallel::mclapply(i_ranges, function(i_section){
      parallel::mclapply(j_ranges, function(j_section){
        #subset to tile
        i_tmp = ppi[i_section]
        j_tmp = ppj[j_section]
        #calculate distances
        dists = spatstat.geom::crossdist(i_tmp, j_tmp)
        #calculate edge correcion
        if(correction %in% c("trans", "translation")){
          edge = spatstat.core::edge.Trans(i_tmp, j_tmp)
        }
        if(correction %in% c("none")){
          edge = dists
        }
        #count edge correction matrix for cells within range r in distance matrix
        counts = sapply(r_range, function(r){sum(edge[which(dists < r)])})
        #remove large distance and edge correction matrix to keep ram usage down
        rm(dists, edge)
        #return counts for tile
        counts
      }) %>% #use 1 core per tile
        #bind all j tiles to data frame
        do.call(rbind.data.frame, .) %>%
        #take column sums and return
        colSums()
    }, mc.cores = cores, mc.preschedule = FALSE)
    #bind i tile counts and sum
    counts = is %>%
      do.call(rbind.data.frame, .) %>%
      colSums() %>%
      unname()
  }
  #if there are less than 10,000 j or i just compute matrix
  if(!(li > 10000 | lj > 10000)){
    #calculate distance matrix
    dists = spatstat.geom::crossdist(ppi, ppj)
    #calculate edge correction
    if(correction %in% c("trans", "translation")){
      edge = spatstat.core::edge.Trans(ppi, ppj)
    }
    if(correction %in% c("none")){
      edge = dists
    }
    #count edge correction matrix for cells within r range in distance matrix
    counts = sapply(r_range, function(r){sum(edge[which(dists < r)])})
  }
  
  #calulate clustering from counted edge corrections with intensities of i and j
  K = (1/(lambdai * lambdaj * area)) * counts
  K
}