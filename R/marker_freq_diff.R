#' Marker Frequency Difference
#'
#' marker_freq_diff allows users to easily calculate the difference in marker abundance (frequency) between different tissue compartments.
#' For us this is typically looking at tumor vs stroma marker abundance.
#' 
#' To calculate a p-value associated with the difference in proportion of positive markers identified in different classifier compartments,
#' we implemented the use of the Fisher's Exact test.
#'
#' @param mif object of class `mif` created with [spatialTIME::create_mif()]
#' @param classifier character that specified which column in the spatial data indicates the different tissue compartments profiled
#' @param ref_level character found in the `classifier` column of spatial data indicating the first value when calculating frequency difference
#' @param diff_level character found in the `classifier` column of spatial data indicating the second value when calculating frequency difference
#' @param mnames character vector of marker names that are found in the spatial data for which to calculate abundance differences
#' @param overwrite boolean/logical for whether to overwrite previously calculated marker frequency differences
#'
#' @return an object of class `mif` with the `derived$frequency_difference` slot filled
#' @export
#'
#' @examples
#' mif <- create_mif(clinical_data = example_clinical %>% 
#'                                  dplyr::mutate(deidentified_id = as.character(deidentified_id)),
#'                                sample_data = example_summary %>% 
#'                                  dplyr::mutate(deidentified_id = as.character(deidentified_id)),
#'                                spatial_list = example_spatial,
#'                                patient_id = "deidentified_id", 
#'                                sample_id = "deidentified_sample")
#' mif = marker_freq_diff(mif = mif,
#'                        classifier = "Classifier.Label",
#'                        ref_level = "Tumor", 
#'                        diff_level = "Stroma",
#'                        mnames = c("CD3..FOXP3.", "CD3..CD8.", "CD3..PD1.", "CD3..PD.L1.",
#'                                   "CD3..Opal.570..Positive", "PD1..Opal.650..Positive"),
#'                        overwrite = TRUE)
marker_freq_diff = function(mif, classifier, ref_level, diff_level, mnames, overwrite = FALSE){
  #make sure that the names marker names is a character vector
  if(!inherits(mnames, "character")){
    stop("Provide a vector of marker names in the spatial files")
  }
  #check classifier/level values
  if(length(classifier) != 1){
    stop("The classifier must be of length 1")
  }
  if(length(ref_level) != 1){
    stop("The ref_level must be of length 1")
  }
  if(length(diff_level) != 1){
    stop("The diff_level must be of length 1")
  }
  
  #get the spatial information
  data = mif$spatial
  #check that the classifier column exists and that the levels exist in samples
  tmp = lapply(data, function(spat){
    data.frame(samp = unique(spat[[mif$sample_id]]), #sample name
               classifier_exists = classifier %in% colnames(spat), 
               ref_exists = ref_level %in% unique(spat[[classifier]]),
               diff_exists = diff_level %in% unique(spat[[classifier]]))
  }) %>%
    do.call(dplyr::bind_rows, .)
  #if missing in all, abort
  if(all(tmp$classifier_exists == FALSE)){
    stop("classifier column not found in any spatial data")
  }
  if(all(tmp$ref_exists == FALSE)){
    stop("ref_level value not found in any spatial data")
  }
  if(all(tmp$diff_exists == FALSE)){
    stop("diff_level value not found in any spatial data")
  }
  #No gc(full=TRUE) here: it forced a full garbage collection on every call, which
  #costs seconds on a large workspace, to reclaim one small audit frame.
  rm(tmp)
  
  #for each of the spatial samples
  out = lapply(data, function(spat){
    #sort marker names to make output later look nicer
    mnames = sort(mnames)
    
    #Per-classifier-level positive counts and totals. Kept as its own object rather
    #than being pivoted straight to one row, because the Fisher tests below need
    #the counts addressable by name -- see the note on that loop.
    agg = spat %>%
      #with the spatial data, select out the markers and the classifier label column
      dplyr::select(dplyr::any_of(c(mnames, classifier))) %>%
      #convert the classifier column to be a factor - possible not needed
      dplyr::mutate(!!classifier := as.factor(get(classifier))) %>%
      #group by the classifier levels and find total number of cells and number of positive cells
      dplyr::group_by(dplyr::across(dplyr::all_of(classifier))) %>%
      dplyr::summarise(Total_Cells = dplyr::n(),
                       dplyr::across(dplyr::any_of(mnames), ~ sum(.x)),
                       .groups = "drop")

    tmp = agg %>%
      #calculate percent of compartment that is positive
      dplyr::mutate(dplyr::across(dplyr::any_of(mnames), ~ .x / Total_Cells * 100, .names = "{col}%")) %>%
      #make data long and sort by column names, currently is 2 rows - want to make 1
      tidyr::pivot_longer(-dplyr::any_of(classifier), names_to = "col", values_to = "val") %>%
      dplyr::arrange(col) %>%
      #make wide to 1 row
      tidyr::pivot_wider(names_from = c(!!classifier, col), values_from = val, names_sep = " ")

    #Fisher's exact test for whether marker positivity differs between the two
    #classifier levels.
    #
    #Two bugs lived here before 2.0.0, both fixed by reading counts out of `agg` by
    #name instead of string-matching the pivoted wide frame:
    #
    #  1. The second row of the table was the compartment TOTAL rather than the
    #     count of marker-NEGATIVE cells, so the margin double-counted the
    #     positives. Every p-value the function produced was wrong -- on
    #     example_spatial[[1]] / CD3..Opal.570..Positive it returned 7.949897e-07
    #     where the correct answer is 6.001431e-07.
    #  2. Selection used dplyr::contains(marker), a substring match. Given two
    #     markers where one name is a prefix of the other -- e.g. "CD3..CD8." and
    #     "CD3..CD8..FOXP3.", both real columns of the shipped data -- the shorter
    #     marker's table absorbed the longer marker's counts and fisher.test() was
    #     handed a 3x2 table, whose r x c p-value was then stored as if it were the
    #     2x2 result for that marker.
    level_row = function(lv) which(as.character(agg[[classifier]]) == lv)
    ft_res = lapply(mnames, function(marker){
      #rows = positive / negative, cols = ref level / diff level
      counts = vapply(c(ref_level, diff_level), function(lv){
        i   = level_row(lv)
        pos = as.numeric(agg[[marker]][i])
        c(positive = pos, negative = as.numeric(agg$Total_Cells[i]) - pos)
      }, numeric(2))
      res = stats::fisher.test(counts)$p.value
      names(res) = paste0(marker, "_p.value")
      return(res)
    }) %>%
      unlist() %>%
      t() %>%
      data.frame()
    
    #calulate the difference in frequency of compartments
    calculated_columns = tmp[,paste0(ref_level, " ",mnames, "%")] - tmp[,paste0(diff_level, " ",mnames, "%")]
    colnames(calculated_columns) = paste0(paste0(ref_level, "_M_", diff_level, " "), mnames) #reference level minus difference level
    #add back to single row above and return
    tmp = dplyr::bind_cols(tmp, calculated_columns) %>%
      dplyr::bind_cols(ft_res) %>%
      dplyr::bind_cols(data.frame(samp = unique(spat[[mif$sample_id]])) %>%
                         dplyr::rename(!!mif$sample_id := 1), .)
    return(tmp)
  }) %>%
    #bind all sample results to single table from list
    do.call(dplyr::bind_rows, .)
  
  #Before 2.0.0 the two inner branches of this were swapped, so overwrite = FALSE
  #(the default) DESTROYED the previous run on a populated mif, and on a fresh mif
  #took the append branch and shipped Run = -Inf -- max() of an empty slot, and of
  #the whole data frame rather than of $Run. write_derived() is the same helper the
  #seven metric functions use, so all of them now share one implementation.
  write_derived(mif, "frequency_difference", out, overwrite)
}