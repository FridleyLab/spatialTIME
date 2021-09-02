#' Subset mif object on cellular level
#'
#' @description This function allows to subset the mif object into compartments. 
#' For instance a mif object includes all cells and the desired analysis is based
#' on only the tumor or stroma compartment then this function will subset the 
#' spatial list to just the cells in the desired compartment 
#' @param mif An MIF object
#' @param classifier Column name for spatial dataframe to subset
#' @param level Determines which level of the classifier to keep.
#' @param markers vector of 
#' @return mif object where the spatial list only as the cell that are the specified level.
#'    
#' @export
#' @examples 
#' #' #Create mif object
#' library(dplyr)
#' x <- create_mif(clinical_data = example_clinical %>% 
#' mutate(deidentified_id = as.character(deidentified_id)),
#' sample_data = example_summary %>% 
#' mutate(deidentified_id = as.character(deidentified_id)),
#' spatial_list = example_spatial,
#' patient_id = "deidentified_id", 
#' sample_id = "deidentified_sample")
#' 
#' markers = c("CD3..Opal.570..Positive","CD8..Opal.520..Positive",
#' "FOXP3..Opal.620..Positive","PDL1..Opal.540..Positive",
#' "PD1..Opal.650..Positive","CD3..CD8.","CD3..FOXP3.")
#' 
#' mif_tumor = subset_mif(mif = x, classifier = 'Classifier.Label', 
#' level = 'Tumor', markers = markers)

subset_mif = function(mif, classifier, level, markers){
  split_spatial = list()
  summary = data.frame()
  for(a in 1:length(mif$spatial)){
    tmp = mif$spatial[[a]] %>% dplyr::filter(get(classifier) == level)
    if(nrow(tmp)>2){
      split_spatial = list.append(split_spatial, tmp)
      names(split_spatial)[length(split_spatial)] = tmp[[mif$sample_id]][1]
      percent = tmp %>% 
        dplyr::select(!!markers) %>% 
        dplyr::summarize_all(~sum(.)) %>%
        dplyr::mutate_all(.funs = ~./nrow(tmp))
      colnames(percent) = paste0(level, ': % ', colnames(percent))
      counts = tmp %>% 
        dplyr::select(!!markers) %>% 
        dplyr::summarize_all(~sum(.)) %>%
        dplyr::mutate(`Total Cells` = nrow(tmp))
      colnames(counts) = paste0(level, ': ', colnames(counts))
      out = c(tmp[[mif$patient_id]][1], 
              tmp[[mif$sample_id]][1],
              unlist(counts) , unlist(percent))
      names(out)[1:2] = c(mif$patient_id,mif$sample_id)
    }
    summary = rbind.data.frame(summary, t(out))
  }
  
  mif_new = create_mif(clinical_data = mif$clinical, sample_data = summary,
                       spatial_list = split_spatial, patient_id =  mif$patient_id, 
                       sample_id =  mif$sample_id)
  
  return(mif_new)
}