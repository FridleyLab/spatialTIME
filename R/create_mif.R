#' Create Multiplex Immunoflourescent object 
#'
#' @description Creates an MIF object for use in spatialIF functions
#' @param clinical_data A data frame containing patient level data with one row
#' per participant. 
#' @param sample_data A data frame containing sample level data with one row per 
#' sample. Should at a minimum contain a 2 columns: one for sample names and 
#' one for the corresponding patient name.
#' @param spatial_list A named list of data frames with the spatial data from 
#'  each sample making up each individual data frame
#' @param patient_id A character string indicating the column name for patient id in 
#'  sample and clinical data frames. 
#' @param sample_id A character string indicating the column name for sample id
#'  in the sample data frame
#' @param clean_columns A logical value indicating if names in spatial data frames
#'  should be rewritten in a cleaner format. Default is FALSE. 
#' 
#' @return Returns a custom MIF
#'    \item{clinical}{Data frame of clinical data}
#'    \item{sample}{Data frame of sample data}
#'    \item{spatial}{Named list of spatial data}
#'    \item{derived}{List of data derived using the MIF object}
#'    
#' @export
create_mif <- function(clinical_data, sample_data, spatial_list = NULL,
                       patient_id = "patient_id", sample_id = "image_tag",
                       clean_columns = FALSE){
  
  sample_data_clean <- sample_data %>% 
    dplyr::full_join(clinical_data %>% 
                dplyr::select(!!patient_id), by = patient_id) %>% 
    dplyr::select(!!patient_id, !!sample_id, dplyr::everything())
  
  clinical_data_clean <- clinical_data %>% 
    dplyr::full_join(sample_data %>% 
                dplyr::group_by_at(dplyr::vars(patient_id)) %>% 
                dplyr::mutate(sample_string = paste0(!!(as.name(sample_id)), collapse = "|")) %>% 
                dplyr::select(!!patient_id, .data$sample_string) %>% 
                dplyr::slice(1), by = patient_id) %>%
    dplyr::select(!!patient_id, .data$sample_string, dplyr::everything())
  
  if(clean_columns== TRUE){
    sample_data <- sample_data_clean %>% 
      janitor::clean_names()
    
    clinical_data <- clinical_data_clean %>% 
      janitor::clean_names()
    
    if(!is.null(spatial_list)){
      
      spatial_list <- lapply(spatial_list, 
                             function(x) {janitor::clean_names(x)})
      
    }
    
  }
  
  if(!is.null(spatial_list) & is.null(names(spatial_list))){
    
    spatial_names <- lapply(spatial_list, function(x) {x[[sample_id]][[1]]})
    spatial_names <- unlist(spatial_names)
    
    names(spatial_list) <- spatial_names
    
  }
  
  if(is.null(spatial_list)) {
    spatial_list <- list(NA)
  }

  mif <- list(clinical = clinical_data,
              sample = sample_data,
              spatial = spatial_list,
              derived = list())
  
  structure(mif, class="mif")
  
  # return(mif)

}