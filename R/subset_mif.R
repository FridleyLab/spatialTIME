#' Subset mif object on cellular level
#'
#' @description This function allows to subset the mif object into compartments. 
#' For instance a mif object includes all cells and the desired analysis is based
#' on only the tumor or stroma compartment then this function will subset the 
#' spatial list to just the cells in the desired compartment 
#' @param mif An MIF object
#' @param classifier Column name for spatial dataframe to subset
#' @param level Determines which level of the classifier to keep.
#' 
#' @return mif object
#'    
#' @export
#'

subset_mif = function(mif, classifier, level){
  split_spatial = list()
  
  for(a in 1:length(mif$spatial)){
    tmp = mif$spatial[[a]] %>% filter(get(classifier) == level)
    if(nrow(tmp)>2){
      split_spatial = rlist::list.append(split_spatial, tmp)
      names(split_spatial)[length(split_spatial)] = tmp[[mif$sample_id]][1]
    }
  }
  
  mif_new = create_mif(clinical_data = mif$clinical, sample_data = mif$sample,
                       spatial_list = split_spatial, patient_id =  mif$patient_id, 
                       sample_id =  mif$sample_id, clean_columns = FALSE)
  
  return(mif_new)
}
