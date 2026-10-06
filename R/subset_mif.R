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
#' @param path for a disk-backed `mif` only (see [mif_to_disk()]), the directory to
#'   write the subset store to. Required in that case and ignored otherwise: the
#'   subset of a cohort too large for memory is generally also too large for memory,
#'   so it is written out sample by sample rather than returned in one piece.
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

subset_mif = function(mif, classifier, level, markers, path = NULL){
  if(!inherits(mif, "mif")){
    stop("mIF should be of class `mif` created with function `create_mif()`")
  }
  #A disk-backed mif was adopted precisely because the cohort does not fit in
  #memory, so filtering it into an in-memory mif would undo that silently, at the
  #worst possible moment -- after the filtering work is already done. Require a
  #destination instead.
  if(is_disk_mif(mif) && is.null(path)){
    stop("`mif` is disk-backed, so `subset_mif()` needs a `path` to write the ",
         "subset store to:\n",
         "    subset_mif(mif, classifier, level, markers, path = \"subset.mif\")\n",
         "  To get an in-memory subset instead, call `collect_mif(mif)` first -- but ",
         "that materialises every sample.", call. = FALSE)
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

  #When writing a store, each filtered frame is written and released immediately, so
  #peak memory stays at one sample rather than the whole subset.
  to_disk = is_disk_mif(mif)
  if(to_disk){
    root = path.expand(path)
    if(file.exists(root)){
      stop("A store already exists at \"", root, "\"; choose another `path`.",
           call. = FALSE)
    }
    dir.create(file.path(root, "spatial"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(root, "derived"), recursive = TRUE, showWarnings = FALSE)
    written = list()
    #Stems for EVERY sample, computed up front, not sample_file_stem() per retained
    #sample: two distinct names can sanitise to one stem ("A/1" and "A 1" both give
    #"A_1"), and then one sample's parquet would overwrite another's while the
    #manifest listed both pointing at the same file. Disambiguating over the full set
    #also keeps a stem stable regardless of which samples survive the filter.
    all_stems = sample_file_stems(names(mif$spatial))
  }

  for(a in seq_len(n_mif_samples(mif))){
    tmp = mif_spatial(mif, a) %>% dplyr::filter(.data[[classifier]] == level)
    if(nrow(tmp) <= 2) next

    sample_name = tmp[[mif$sample_id]][1]
    patient = mif$sample[[mif$patient_id]][mif$sample[[mif$sample_id]] == sample_name]
    #Guard both directions: no match, and more than one row for this sample id.
    #Length > 1 used to shift every subsequent value in the row by one, silently.
    patient = if(length(patient) < 1) NA else patient[1]

    if(to_disk){
      rel = file.path("spatial", paste0(all_stems[[a]], ".parquet"))
      pq_write(tmp, file.path(root, rel))
      written[[length(written) + 1]] = list(sample = sample_name, rel = rel)
    } else {
      split_spatial = list.append(split_spatial, tmp)
      names(split_spatial)[length(split_spatial)] = sample_name
    }

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

  if(to_disk){
    if(!length(written)){
      #The store would have no samples at all; do not leave an unopenable directory
      #behind for the user to puzzle over.
      unlink(root, recursive = TRUE)
      stop("No sample had more than 2 cells at level \"", level, "\" of \"",
           classifier, "\", so there is nothing to write to \"", root, "\".",
           call. = FALSE)
    }
    rel     = vapply(written, function(w) w$rel, character(1))
    samples = vapply(written, function(w) w$sample, character(1))
    probes  = lapply(file.path(root, rel), pq_probe)
    store = new_mif_store(
      files   = rel, samples = samples,
      nrow    = vapply(probes, function(p) p$nrow, numeric(1)),
      columns = collapse_schema(lapply(probes, function(p) p$columns)),
      types   = collapse_schema(lapply(probes, function(p) p$types)),
      overlay = NA_character_, root = root)
    attr(store, "overlay_columns") = rep(list(character(0)), length(rel))
    mif_new = structure(list(clinical = mif$clinical, sample = summary,
                             spatial = store, derived = list(),
                             patient_id = mif$patient_id,
                             sample_id = mif$sample_id),
                        class = "mif")
    write_store_tables(root, mif_new)
    write_manifest(root, mif_new, store,
                   vapply(file.path(root, rel), file.size, numeric(1)))
    return(mif_new)
  }

  mif_new = create_mif(clinical_data = mif$clinical, sample_data = summary,
                       spatial_list = split_spatial, patient_id =  mif$patient_id,
                       sample_id =  mif$sample_id)

  return(mif_new)
}