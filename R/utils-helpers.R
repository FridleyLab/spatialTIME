# Internal helpers for spatialTIME.
#
# Only helpers reachable from an exported function belong here. The gen-1
# permutation engine that used to live in this file (uni_Rip_K, bi_Rip_K,
# uni_NN_G, bi_NN_G_sample and their K_out/G_out/uni_K/bi_K/uni_G/bi_G/perm_data
# workers) was reachable only from compute_metrics(), which was never exported;
# all of it was removed in 2.0.0 along with get_exactK(), dix_s_z() and
# dix_s_c(). See NEWS.md.

get_bi_rows = function(data, markers){
  data %>%
    dplyr::mutate(cell = 1:dplyr::n()) %>%
    dplyr::select(cell, xloc, yloc, !!markers) %>%
    dplyr::filter(!(get(markers[1]) == 1 & get(markers[2]) == 1)) %>%
    tidyr::gather("Marker", "Positive", -cell, -xloc, -yloc) %>%
    dplyr::filter(Positive == 1) %>%
    dplyr::mutate(Marker = factor(Marker, levels = markers))
}

list.append = function(list, new){
  new_list = list
  new_list[[length(new_list) + 1]] = new
  return(new_list)
}

#get tile counts
getTile = function(slide, l, size){
  if(slide == 1){
    v = rep(TRUE, l)
    return(list(v))
  }
  lapply(1:slide, function(s){
    if(s == 1){
      w = 1:size
    }
    if(s != 1 & s != slide){
      w = (((s-1)*size)+1):((s)*size)
    }
    if(s == slide){
      w = (((s-1)*size)+1):(l)
    }
    v = rep(FALSE, l)
    v[w] = TRUE
    v
  })
}

calculateK = function(i_dat, j_dat, anchor, counted, area, win, big, r_range, edge_correction, cores){
  li = nrow(i_dat)
  lj = nrow(j_dat)
  #find intensity of each marker
  lambdai = li/area
  lambdaj = lj/area
  #point pattern for anchor for distance matrix
  ppi = spatstat.geom::ppp(x = i_dat[,1],
                           y = i_dat[,2],
                           window = win)
  #point pattern for counted for distance matrix
  ppj = spatstat.geom::ppp(x = j_dat[,1],
                           y = j_dat[,2],
                           window = win)
  if(li > big | lj > big){
    #find the number of tiles of distance matrix in columns and rows
    i_slides = ceiling(li / big)
    j_slides = ceiling(lj / big)
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
        if(edge_correction %in% c("border")){
          bI = spatstat.geom::bdist.points(i_tmp)
          bcloseI = bI[spatstat.geom::crosspairs(i_tmp, j_tmp, max(r_range), what = "ijd")$i]
          RS = spatstat.explore::Kount(spatstat.geom::crosspairs(i_tmp, j_tmp, max(r_range), what = "ijd")$d,
                                    spatstat.geom::crosspairs(i_tmp, j_tmp, max(r_range), what = "ijd")$i,
                                    bI, spatstat.geom::handle.r.b.args(r_range, breaks=NULL, win, rmaxdefault = max(r_range)))
          return(RS)
        }
        #calculate distances
        dists = spatstat.geom::crossdist(i_tmp, j_tmp)
        rmv_i = rowSums(dists < max(r_range)) != 0
        rmv_j = colSums(dists < max(r_range)) != 0
        i_tmp = i_tmp[rmv_i] #I
        j_tmp = j_tmp[rmv_j] #J
        dists = spatstat.geom::crossdist(i_tmp, j_tmp)
        if(0 %in% dim(dists)){
          return(rep(0, length(r_range)))
        }
        #calculate edge correcion
        if(edge_correction %in% c("trans", "translation")){
          edge = spatstat.explore::edge.Trans(i_tmp, j_tmp)
          #count edge correction matrix for cells within range r in distance matrix
          counts = sapply(r_range, function(r){sum(edge[which(dists < r)])})
          #remove large distance and edge correction matrix to keep ram usage down
          rm(dists, edge)
          #return counts for tile
          return(counts)
        }
        if(edge_correction %in% c("none")){
          counts = cumsum(spatstat.univar::whist(spatstat.geom::crosspairs(i_tmp, j_tmp, max(r_range), what = "ijd")$d,
                                               spatstat.geom::handle.r.b.args(r_range, breaks=NULL, win, rmaxdefault = max(r_range))$val))
          return(counts)
        }
      }) 
      
      if(edge_correction == "border"){
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
    },mc.cores = cores, mc.preschedule = FALSE, mc.allow.recursive = TRUE)
    
    if(edge_correction == "border"){
      num = lapply(counts, function(j_big){
        j_big[[1]]
      }) %>%
        do.call(rbind.data.frame, .) %>%
        colSums() %>%
        unname()
      den = ppi$n
      estimatedK = num / (lambdaj * den)
    } else {
      counts = counts %>%
        #bind i tile counts and sum
        do.call(rbind.data.frame, .) %>%
        colSums() %>%
        unname()
      #calulate clustering from counted edge corrections with intensities of i and j
      estimatedK = (1/(lambdai * lambdaj * area)) * counts
    }
  }
  #if there are less than 10,000 j or i just compute matrix
  if(!(li > big | lj > big)){
    #calculate distance matrix
    sp_tmp2 = as.data.frame(rbind(i_dat, j_dat)) %>% 
      dplyr::mutate(Marker = c(rep(anchor, nrow(i_dat)), rep(counted, nrow(j_dat))))
    pp_obj = spatstat.geom::ppp(x = sp_tmp2$xloc, y = sp_tmp2$yloc, window = win, marks = factor(sp_tmp2$Marker))
    estimatedK = spatstat.explore::Kcross(pp_obj,
                                          i = unique(sp_tmp2$Marker)[1],
                                          j = unique(sp_tmp2$Marker)[2],
                                          r = r_range,
                                          correction = edge_correction) %>%
      data.frame() %>%
      dplyr::pull(3)
  }
  
  return(estimatedK)
}
