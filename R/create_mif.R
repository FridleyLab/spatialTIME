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
#' @param spatialexp A SpatialExperiment object to create a MIF from a 
#' \href{https://www.bioconductor.org/packages/release/bioc/vignettes/SpatialExperiment/inst/doc/SpatialExperiment.html}{SpatialExperiment} 
#' Bioconductor object.
#' 
#' @import dplyr
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
#' # Create mif object
#' 
#' x <- create_mif(
#'   clinical_data = example_clinical,
#'   sample_data = example_summary,
#'   spatial_list = example_spatial,
#'   patient_id = "deidentified_id", 
#'   sample_id = "deidentified_sample"
#'   )

create_mif <- function(clinical_data, 
                       sample_data, 
                       spatial_list = NULL,
                       patient_id = "patient_id", 
                       sample_id = "image_tag"){
  
  
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
    "Each item in spatial_list must be named" = 
     all(!sapply(names(spatial_list), is.null))
    )
  
  sample_data_clean <- sample_data %>% 
    dplyr::full_join(clinical_data %>% 
                dplyr::select(!!patient_id), by = patient_id) %>% 
    dplyr::select(dplyr::all_of(c(!!patient_id, !!sample_id)), dplyr::everything()) %>% 
    dplyr::group_by_at(patient_id) %>% 
    dplyr::mutate(sample_string = paste0(!!(as.name(sample_id)), collapse = "|")) %>% 
    dplyr::select(dplyr::all_of(c(!!patient_id, 'sample_string'))) %>% 
    dplyr::slice(1)
  
  clinical_data_clean <- clinical_data %>% 
    dplyr::full_join(sample_data_clean, by = patient_id) %>%
    dplyr::select(dplyr::all_of(c(!!patient_id, 'sample_string')), dplyr::all_of(dplyr::everything()))
  
  
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

#' Create Multiplex Immunoflourescent object from a SpatialExperiment object
#'
#' @description Creates an MIF object for use in spatialIF functions using a 
#' SpatialExperiment object as input
#' @param spatial_exp A SpatialExperiment object to create a MIF from a 
#' \href{https://www.bioconductor.org/packages/release/bioc/vignettes/SpatialExperiment/inst/doc/SpatialExperiment.html}{SpatialExperiment} 
#' Bioconductor object.
#' 
#' @import dplyr
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
#' # Create mif object
#' 
#' BiocManager::install("VectraPolarisData")
#' ovarian <- VectraPolarisData::HumanOvarianCancerVP()
#' 
#' mif <- spatial_exp_to_mif(spatial_exp = ovarian)
#'
#'   
spatial_exp_to_mif <- function(spatial_exp, 
                               markers = c("phenotype_cd68", "phenotype_cd3", "phenotype_cd8"), 
                               x_coord = "cell_x_position",
                               y_coord = "cell_y_position", 
                               patient_id = "patient_id",
                               sample_id = "sample_id"){
  
  stopifnot(
    "spatial_exp needs to be of class 'SpatialExperiment'" = 
      class(spatial_exp) == "SpatialExperiment"
  )
  
  # extract information from spatial experiment object
  colData_df <- colData(spatial_exp)
  
  spatialcoords_df <- spatialCoords(spatial_exp)
  
  clinical <- metadata(spatial_exp)$clinical_data
  
  #collapse spatial data
  cell_level <- data.frame(
    cbind(spatialcoords_df, colData_df),
    check.names = FALSE
  )
  
  # Convert to spatialTIME format -------------------------------------------
  
  #spatial data frame
  spat_df = cell_level %>%
    mutate(across(any_of(markers), #for markers
                  ~ifelse(grepl("\\+", .x,), #check if has "+"
                          1, #convert to 1 if yes
                          0) #make 0 if no
    )
    ) %>%
    filter(sample_id %in% meta$sample_id) %>% # and samples with metadata
    dplyr::rename("x" = x_coord,
                  "y" = y_coord)
  
  # convert df to list by samples
  spat_list <- split(spat_df, spat_df$sample_id)
  
  # create summary by sample
  spat_summ <- spat_df %>%
    group_by(sample_id) %>% # for each sample
    summarise(across(any_of(markers), # across all markers
                     ~ sum(.x), # count the number positive
                     .names = "{col} Cells"),
              `Total Cells` = n()) %>% # total number of cells in sample
    mutate(across(contains(markers), # across all markers
                  ~ .x / `Total Cells` * 100, #calculate percent of total
                  .names = "Percent {col} Cells"),
           .before = `Total Cells`)
  
  #create mif
  mif <- create_mif(
    clinical_data = clinical,
    sample_data = spat_summ,
    spatial_list = spat_list,
    patient_id = "patient_id", #linker from summary to clinical
    sample_id = "sample_id" #linker from spatial to summary
  )
  
  return(mif)
}