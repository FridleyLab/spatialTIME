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
#' mif <- create_mif(
#'   clinical_data = example_clinical,
#'  sample_data = example_summary,
#'   spatial_list = example_spatial,
#'   patient_id = "deidentified_id",
#'   sample_id = "deidentified_sample"
#' )
#'

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

  clinical_data <- clinical_data %>%
    dplyr::mutate(patient_id = as.character(!!(as.name(patient_id))))

  sample_data <- sample_data %>%
    dplyr::mutate(patient_id = as.character(!!(as.name(patient_id))))

  sample_data_clean <- sample_data %>%
    dplyr::full_join(clinical_data %>%
                dplyr::select(!!patient_id), by = patient_id) %>%
    dplyr::select(dplyr::all_of(c(patient_id, !!sample_id)), dplyr::everything()) %>%
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
#' Bioconductor SpatialExperiment object as input.
#' @param spatial_exp A SpatialExperiment object to create a MIF from a
#' \href{https://www.bioconductor.org/packages/release/bioc/vignettes/SpatialExperiment/inst/doc/SpatialExperiment.html}{SpatialExperiment}
#' Bioconductor object.
#' @param markers A character vector that contains the name(s) of the columns that contains +/- for markers. This will be used to recode
#' these columns into logical vectors for presence or absence of the marker. By default looks for a + to indicate presence. But this can be
#' altered using the
#' @param x_coord The name of the column in `spatialCoords` that contains the x coordinates for the cell. Default is "cell_x_position".
#' @param y_coord The name of the column in `spatialCoords` that contains the y coordinates for the cell. Default is "cell_y_position".
#' @param patient_id The name of the column in `colData` and `metadata` that contains the sample IDs. Default is "patient_id",
#' @param sample_id The name of the column in `colData` and `spatialCoords` that contains the sample IDs. Default is "sample_id".
#' @param cols_to_keep (Optional) By default only the patient_id, markers, sample_id, and x and y coordinate
#' columns will be kept. But if there is additional information that needs to be kept, you can add those column names
#' here. Or say "everything" if you want all columns to be kept. This may not be advisable for large datasets.
#' @param marker_pos_regex By default will identify all markers columns as having a marker if a `+` is found in the data. But with this variable
#' You can choose a different character to sort in presence/absence. Anything with the regex character given to this parameter will be labeled as
#' positive and all other cells will be labeled as negative.
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
#' if (!("VectraPolarisData" %in% installed.packages())) {
#'  BiocManager::install("VectraPolarisData")
#' }
#' ovarian <- VectraPolarisData::HumanOvarianCancerVP()
#'
#' ova_mif <- spatial_exp_to_mif(spatial_exp = ovarian,
#'                               patient_id = "sample_id",
#'                               markers = c("phenotype_cd68", "phenotype_cd3", "phenotype_cd8"))
#'
#' spe_lung <- VectraPolarisData::HumanLungCancerV3()
#'
#'
#' spe_mif <- spatial_exp_to_mif(spatial_exp = spe_lung,
#'                               patient_id = "slide_id",
#'                               sample_id = "slide_id",
#'                               markers = c(
#'                                 "phenotype_cd4",
#'                                 "phenotype_cd8",
#'                                 "phenotype_cd14",
#'                                 "phenotype_cd19",
#'                                 "phenotype_ck",
#'                                 "phenotype_other"))
#'
spatial_exp_to_mif <- function(spatial_exp,
                               markers,
                               x_coord = "cell_x_position",
                               y_coord = "cell_y_position",
                               patient_id = "patient_id",
                               sample_id = "sample_id",
                               cols_to_keep = NULL,
                               marker_pos_regex = "\\+"){

  stopifnot(
    "spatial_exp needs to be of class 'SpatialExperiment'" =
      class(spatial_exp) == "SpatialExperiment"
  )

  # extract information from spatial experiment object
  colData_df <- colData(spatial_exp)

  spatialcoords_df <- spatialCoords(spatial_exp)

  clinical <- metadata(spatial_exp)$clinical_data %>%
    dplyr::rename(patient_id = dplyr::all_of(patient_id))

  stopifnot(
  "'markers' needs to be a character vector that indicates column names in spatialCoords" =
    is.vector(markers),
  "'x_coord' needs to be a character that is a column name in spatialCoords data" =
    x_coord %in% colnames(spatialcoords_df),
  "'y_coord' needs to be a character that is a column name in the spatialCoords data" =
    y_coord %in% colnames(spatialcoords_df),
  "'sample_id' needs to be a character that is a column name in colData and spatialCoords " =
    all(sample_id %in% colnames(colData_df)),
  "All characters in the markers vector need to be columns name in colData" =
    all(markers %in% colnames(colData_df)),
  "marker_pos_regex needs to be a string that can be used for regex pattern identification" =
    is.character(marker_pos_regex),
  "cols_to_keep needs to be a chracter vector of column names within colData" =
    all(cols_to_keep %in% colnames(colData_df))
  )

  #collapse spatial data
  cell_level <- data.frame(
    cbind(spatialcoords_df, colData_df),
    check.names = FALSE
  )

  if (is.null(cols_to_keep)) {
    cell_level <- cell_level %>%
      dplyr::select(dplyr::all_of(c(sample_id, x_coord, y_coord, markers)))
  } else {
    cell_level <- cell_level %>%
      dplyr::select(dplyr::all_of(c(sample_id, x_coord, y_coord, markers, cols_to_keep)))
  }


  # Convert to spatialTIME format -------------------------------------------
  sample_ids <- dplyr::pull(clinical, patient_id)

  # Spatial data frame
  spat_df <- cell_level %>%
    dplyr::rename(sample_id = dplyr::all_of(sample_id)) %>%
    dplyr::mutate(across(any_of(markers), #for markers
                  ~ifelse(grepl(marker_pos_regex, .x,), #check if has "+" or whatever is specified
                          1, #convert to 1 if yes
                          0) #make 0 if no
    )
    ) %>%
    dplyr::filter(sample_id %in% sample_ids) %>% # and samples with metadata
    dplyr::rename("x" = x_coord,
                  "y" = y_coord)  %>%
    dplyr::mutate(patient_id = sample_id)

  # convert df to list by samples
  spat_sample_ids <- spat_df %>% dplyr::pull(sample_id)
  spat_list <- split(spat_df, spat_sample_ids)

  # create summary by sample
  spat_summ <- spat_df %>%
    dplyr::group_by(sample_id) %>% # for each sample
    dplyr::summarise(across(any_of(markers), # across all markers
                     ~ sum(.x), # count the number positive
                     .names = "{col} Cells"),
              `Total Cells` = n()) %>% # total number of cells in sample
    dplyr::mutate(across(contains(markers), # across all markers
                  ~ .x / `Total Cells` * 100, #calculate percent of total
                  .names = "Percent {col} Cells"),
           .before = `Total Cells`) %>%
    dplyr::mutate(patient_id = sample_id)

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


utils::globalVariables(c("Total Cells", "colData", "metadata", "spatialCoords"))
