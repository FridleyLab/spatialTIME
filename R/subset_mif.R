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
  if(!inherits(mif, "mif")){
    stop("mIF should be of class `mif` created with function `create_mif()`")
  }
  split_spatial = list()
  #Collect one row per RETAINED sample and bind at the end.
  #
  #Before 2.0.0 this loop assigned `out` only inside `if(nrow(tmp) > 2)` but ran
  #`summary = rbind.data.frame(summary, t(out))` unconditionally, which failed two
  #different ways: if the first sample had <=2 cells at the requested level the
  #function died with "object 'out' not found", and if a LATER sample did, the
  #previous sample's row was silently appended a second time -- leaving a summary
  #row that corresponded to no retained spatial frame. The silent case was the
  #dangerous one.
  summary_rows = list()

  for(a in seq_along(mif$spatial)){
    tmp = mif$spatial[[a]] %>% dplyr::filter(.data[[classifier]] == level)
    if(nrow(tmp) <= 2) next

    sample_name = tmp[[mif$sample_id]][1]
    patient = mif$sample[[mif$patient_id]][mif$sample[[mif$sample_id]] == sample_name]
    #Guard both directions: no match, and more than one row for this sample id.
    #Length > 1 used to shift every subsequent value in the row by one, silently.
    patient = if(length(patient) < 1) NA else patient[1]

    split_spatial = list.append(split_spatial, tmp)
    names(split_spatial)[length(split_spatial)] = sample_name

    pos = tmp %>%
      dplyr::select(dplyr::all_of(markers)) %>%
      dplyr::summarize(dplyr::across(dplyr::everything(), ~ sum(.x)))

    counts = pos
    counts$`Total Cells` = nrow(tmp)
    colnames(counts) = paste0(level, ': ', colnames(counts))

    #x100. These columns are labelled "%", and marker_freq_diff() puts true
    #percentages in its "%" columns; before 2.0.0 this one stored a proportion, so
    #the two exported functions disagreed on what "%" meant.
    percent = pos %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), ~ .x / nrow(tmp) * 100))
    colnames(percent) = paste0(level, ': % ', colnames(percent))

    #A typed one-row data frame. The old `c(patient, id, unlist(counts),
    #unlist(percent))` coerced everything to character at the c(), so the whole
    #summary table came back as strings like "3611" and "0.00249238438105788".
    summary_rows[[length(summary_rows) + 1]] = dplyr::bind_cols(
      stats::setNames(
        data.frame(patient, sample_name, stringsAsFactors = FALSE),
        c(mif$patient_id, mif$sample_id)
      ),
      counts, percent
    )
  }

  if(length(summary_rows)){
    summary = do.call(dplyr::bind_rows, summary_rows)
  } else {
    #Nothing survived the filter. A bare data.frame() has no columns, so create_mif
    #would reject it for a missing patient_id; emit a zero-row frame with the right
    #shape instead, and say so rather than handing back a silently empty mif.
    warning("No sample had more than 2 cells at level \"", level, "\" of \"",
            classifier, "\"; returning a mif with no spatial data.", call. = FALSE)
    summary = stats::setNames(
      data.frame(character(0), character(0), stringsAsFactors = FALSE),
      c(mif$patient_id, mif$sample_id)
    )
  }

  mif_new = create_mif(clinical_data = mif$clinical, sample_data = summary,
                       spatial_list = split_spatial, patient_id =  mif$patient_id,
                       sample_id =  mif$sample_id)

  return(mif_new)
}