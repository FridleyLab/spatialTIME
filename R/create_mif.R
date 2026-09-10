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
#' 
#' @return Returns a custom MIF
#'    \item{clinical}{Data frame of clinical data}
#'    \item{sample}{Data frame of sample data}
#'    \item{spatial}{Named list of spatial data}
#'    \item{derived}{List of data derived using the MIF object}
#'    \item{patient_id}{The column name for sample id
#'  in the sample data frame with the clinical data}
#'    \item{sample_id}{The column name for sample id
#'  in the sample data frame to merge with the spatial data}
#'    
#' @export
#' @examples
#' #Create mif object
#' library(dplyr)
#' x <- create_mif(clinical_data = example_clinical %>% 
#' mutate(deidentified_id = as.character(deidentified_id)),
#' sample_data = example_summary %>% 
#' mutate(deidentified_id = as.character(deidentified_id)),
#' spatial_list = example_spatial,
#' patient_id = "deidentified_id", 
#' sample_id = "deidentified_sample")

create_mif <- function(clinical_data, sample_data, spatial_list = NULL,
                       patient_id = "patient_id", sample_id = "image_tag"){
  #checks from Candace Savonen
  stopifnot(
    "clinical_data must be a data frame"= is.data.frame(clinical_data),
    "sample_data must be a data frame" = is.data.frame(sample_data),
    "spatial_list must be a list of data frames" = is.list(spatial_list),
    "All items in 'spatial_list' must be a data frame" = all(sapply(spatial_list, is.data.frame)),
    "patient_id must be a character indicating a column name" = is.character(patient_id),
    "sample_id must be a character indicating a column name" = is.character(sample_id),
    "The column specified by 'patient_id' could not be found in 'clinical_data'" = 
      patient_id %in% colnames(clinical_data),
    "The column specified by 'patient_id' could not be found in 'sample_data'" = 
      patient_id %in% colnames(sample_data),
    "The column specified by 'sample_id' could not be found in 'sample_data'" = 
      sample_id %in% colnames(sample_data), 
    #`names()` returns a character vector, so is.null() on each element was always
    #FALSE and this check could never fail; when names() was NULL, sapply() over it
    #gave list() and all(logical(0)) is TRUE. Auto-naming below handles the
    #genuinely unnamed case, so only reject partial/blank names here.
    "Each item in spatial_list must be named" =
      is.null(names(spatial_list)) ||
        !any(is.na(names(spatial_list)) | !nzchar(names(spatial_list)))
  )
  
  #Removed in 2.0.0: a `sample_data_clean`/`clinical_data_clean` pair was computed
  #here and then never used -- the returned mif below carries the RAW inputs. So the
  #documented `sample_string` column never existed in the product, the work was
  #wasted on every call, and worst of all the discarded full_join() could FAIL on
  #incompatible id types, making create_mif() refuse to build a mif because of a
  #join whose result it threw away. That is why every example in this package used
  #to carry `mutate(deidentified_id = as.character(deidentified_id))`.

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
              derived = list(),
              patient_id = patient_id,
              sample_id = sample_id)
  
  structure(mif, class="mif")
  
  # return(mif)

}