# Internal helpers for spatialTIME.
#
# Only helpers reachable from an exported function belong here. The gen-1
# permutation engine that used to live in this file (uni_Rip_K, bi_Rip_K,
# uni_NN_G, bi_NN_G_sample and their K_out/G_out/uni_K/bi_K/uni_G/bi_G/perm_data
# workers) was reachable only from compute_metrics(), which was never exported;
# all of it was removed in 2.0.0 along with get_exactK(), dix_s_z() and
# dix_s_c(). getTile() and calculateK() went with the tiled Ripley's K
# implementation they existed to serve -- see R/utils-k-engine.R. See NEWS.md.

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
