#' Create Multiplex Immunoflourescent object 
#'
#' @description Creates an MIF object for use in spatialIF functions
#' @param clinical_data A data frame containing patient level data with one row
#' per participant. 
#' @param sample_data A data frame containing sample level data with one row per 
#' sample. Should at a minimum contain a 2 columns: one for sample names and 
#' one for the corresponding patient name.
#' @param spatial_list Either a named list of data frames with the spatial data
#'  from each sample making up each individual data frame, or a **named character
#'  vector of paths to parquet files**, one sample per file, which builds a
#'  disk-backed mif that reads each sample only when a computation needs it. The
#'  files are referenced where they are and never copied or modified. Use this for
#'  whole slide images, where holding every sample in memory at once is what stops
#'  a multi-core run from starting -- see [mif_to_disk()] for the sizes involved
#'  and for building a store from data already in memory. If the vector is unnamed,
#'  names are taken from each file's `sample_id` column.
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
  #A character `spatial_list` means "these are parquet files on disk". This is a
  #pure addition: the stopifnot() below has always required a list of data frames,
  #so a character vector was an error in every previous version and no existing
  #call can change behaviour. See ?mif_to_disk for what a disk-backed mif is for.
  if(is.character(spatial_list)){
    return(create_mif_from_paths(clinical_data, sample_data, spatial_list,
                                 patient_id, sample_id))
  }
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


#' Build a disk-backed mif from a vector of parquet paths
#'
#' Reference mode: the files stay where they are and are never copied or modified,
#' which is the point -- for whole slide images they are typically the primary data
#' and may be tens of gigabytes. Only each file's footer is read here, so building
#' the mif costs about a millisecond per sample regardless of how large they are.
#'
#' Sample names: supplied names are trusted. When the vector is unnamed, the
#' `sample_id` column is read from each file and its single unique value used, which
#' matches how the in-memory path derives names (`x[[sample_id]][[1]]`) and catches
#' a file holding more than one sample. That case matters because every metric takes
#' the observation window to be the convex hull of *all* cells in a file, so two
#' samples in one file silently produces the wrong window for both.
#'
#' @keywords internal
#' @noRd
create_mif_from_paths <- function(clinical_data, sample_data, spatial_list,
                                  patient_id, sample_id) {
  stopifnot(
    "clinical_data must be a data frame" = is.data.frame(clinical_data),
    "sample_data must be a data frame" = is.data.frame(sample_data),
    "patient_id must be a character indicating a column name" = is.character(patient_id),
    "sample_id must be a character indicating a column name" = is.character(sample_id),
    "The column specified by 'patient_id' could not be found in 'clinical_data'" =
      patient_id %in% colnames(clinical_data),
    "The column specified by 'patient_id' could not be found in 'sample_data'" =
      patient_id %in% colnames(sample_data),
    "The column specified by 'sample_id' could not be found in 'sample_data'" =
      sample_id %in% colnames(sample_data)
  )
  if (!length(spatial_list)) {
    stop("`spatial_list` is an empty character vector; no spatial files to use.",
         call. = FALSE)
  }
  missing <- spatial_list[!file.exists(spatial_list)]
  if (length(missing)) {
    stop("Spatial file", if (length(missing) > 1) "s" else "", " not found:\n  ",
         paste(utils::head(missing, 5), collapse = "\n  "),
         if (length(missing) > 5) paste0("\n  ... and ", length(missing) - 5, " more"),
         call. = FALSE)
  }
  bad <- spatial_list[!grepl("\\.parquet$", spatial_list, ignore.case = TRUE)]
  if (length(bad)) {
    stop("A character `spatial_list` must name parquet files; these do not look ",
         "like parquet:\n  ", paste(utils::head(bad, 5), collapse = "\n  "), "\n",
         "  To build a store from data already in memory use `mif_to_disk()`.",
         call. = FALSE)
  }

  probes <- lapply(spatial_list, pq_probe)
  for (k in seq_along(probes)) {
    if (!sample_id %in% probes[[k]]$columns) {
      stop("\"", spatial_list[k], "\" has no column named \"", sample_id,
           "\" (the `sample_id`).\n  It has: ",
           abbreviate_names(probes[[k]]$columns), call. = FALSE)
    }
  }

  samples <- names(spatial_list)
  if (is.null(samples) || anyNA(samples) || !all(nzchar(samples))) {
    samples <- vapply(seq_along(spatial_list), function(k) {
      v <- unique(pq_read(spatial_list[k], sample_id)[[sample_id]])
      if (length(v) != 1L) {
        stop("\"", spatial_list[k], "\" contains ", length(v), " distinct values of ",
             "\"", sample_id, "\", but a spatial file must hold exactly one sample ",
             "-- every metric uses the convex hull of all of a file's cells as the ",
             "observation window.\n  Split it by ", sample_id,
             ", or name `spatial_list` yourself to override this check.",
             call. = FALSE)
      }
      as.character(v)
    }, character(1))
  }
  if (anyDuplicated(samples)) {
    dup <- unique(samples[duplicated(samples)])
    stop("Duplicate sample name", if (length(dup) > 1) "s" else "", ": ",
         paste0("\"", dup, "\"", collapse = ", "), ".", call. = FALSE)
  }

  store <- new_mif_store(
    files   = unname(spatial_list),
    samples = samples,
    nrow    = vapply(probes, function(p) p$nrow, numeric(1)),
    columns = collapse_schema(lapply(probes, function(p) p$columns)),
    types   = collapse_schema(lapply(probes, function(p) p$types)),
    overlay = NA_character_,
    root    = NA_character_
  )
  attr(store, "overlay_columns") <- rep(list(character(0)), length(spatial_list))

  structure(list(clinical = clinical_data,
                 sample = sample_data,
                 spatial = store,
                 derived = list(),
                 patient_id = patient_id,
                 sample_id = sample_id),
            class = "mif")

}