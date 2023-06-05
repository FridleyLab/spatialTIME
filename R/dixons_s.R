#' Dixon's S Segregation Statistic
#'
#' @description This function processes the spatial files in the mif object,
#' requiring a column that distinguishes between different groups i.e. tumor and 
#' stroma
#' @param mif An MIF object
#' @param mnames vector of markers corresponding to spatial columns to check Dixon's S between
#' @param num_permutations Numeric value indicating the number of permutations used. 
#'  Default is 1000.
#' @param type a character string for the type that is wanted in the output which can
#' be "Z" for z-statistic results or "C" for Chi-squared statistic results
#' @param workers Integer value for the number of workers to spawn
#' @param overwrite Logical value determining if you want the results to replace the 
#' current output (TRUE) or be to be appended (FALSE).
#' @param xloc a string corresponding to the x coordinates. If null the average of 
#' XMin and XMax will be used 
#' @param yloc a string corresponding to the y coordinates. If null the average of 
#' YMin and YMax will be used 
#' @importFrom magrittr %>%
#' 
#' @return Returns a data frame for Z-statistic
#'    \item{From}{}
#'    \item{To}{}
#'    \item{Obs.Count}{}
#'    \item{Exp. Count}{}
#'    \item{S}{}
#'    \item{Z}{}
#'    \item{p-val.Z}{}
#'    \item{p-val.Nobs}{}
#'    \item{Marker}{}
#'    \item{Classifier Labeled Column Counts}{}
#'    \item{Image.Tag}{}
#' @return Returns a data frame for C-statistic
#'    \item{Segregation}{}
#'    \item{df}{}
#'    \item{Chi-sq}{}
#'    \item{P.asymp}{}
#'    \item{P.rand}{}
#'    \item{Marker}{}
#'    \item{Classifier Labeled Column Counts}{}
#'    \item{Image.Tag}{}
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
#' @export

dixons_s = function(mif, mnames, num_permutations = 1000, type = c("Z", "C"),
                    workers = 1, overwrite = FALSE, xloc = NULL, yloc = NULL){
  
  if(!is(mnames, "character")){
    stop("Provide a vector of marker names in the spatial files")
  }
  data = mif$spatial
  #filter names of markers because order doesn't matter, still computes both ways
  mnames = mnames %>%
    expand.grid(., .) %>%
    dplyr::filter(Var1 != Var2) %>% 
    dplyr::rowwise() %>%
    dplyr::mutate(Var3 = paste0(sort(c(Var1, Var2)), collapse = ",")) %>%
    dplyr::distinct(Var3, .keep_all = TRUE) %>%
    dplyr::select(1, 2) %>%
    dplyr::ungroup()
  
  #per spatial file
  out = parallel::mclapply(data, function(spat){
    #if locations not provided
    if(is.null(xloc)){
      spat$xloc = (spat$XMin+spat$XMax)/2
    } else {
      spat$xloc = spat[[xloc]]
    }
    if(is.null(yloc)){
      spat$yloc = (spat$YMin+spat$YMax)/2
    } else {
      spat$yloc = spat[[yloc]]
    }
    
    res = parallel::mclapply(1:nrow(mnames), function(r){
      markers = as.character(unlist(mnames[r,]))
      df = spat %>%
        dplyr::select(xloc, yloc, dplyr::any_of(markers)) %>%
        dplyr::filter(!(get(markers[1]) == 1 & get(markers[2]) == 1)) %>% #get rid of dual positives
        tidyr::gather("Marker", "Positive", -xloc, -yloc) %>%
        dplyr::mutate(Marker = factor(Marker, levels = markers)) %>%
        dplyr::filter(Positive == 1)
      
      df_tab = table(df$Marker)
      if(nrow(df_tab) == 0 | nrow(df_tab) == 1 | TRUE %in% (df_tab < 3)){ #set minimum number of cells in each group to 3
        final_df_z = expand.grid(From = names(df_tab),
                               To = names(df_tab)) %>%
          dplyr::mutate(`Image Location` = spat$`Image Location`[1], #
                 !!markers[1] := df_tab[1],
                 !!markers[2] := df_tab[2])
        final_df_c = data.frame(df = rep(NA, 3),
                                `Chi-sq` = rep(NA, 3),
                                P.asymp = rep(NA, 3),
                                P.rand = rep(NA, 3),
                                check.names = FALSE)
        rownames(final_df_c) = c("Overall segregation", paste("From", markers))
        colnames(final_df_c)[1] = mif$sample_id
        return(list(tablaZ = final_df_z,
                    tablaC = final_df_c))
      }
      dixon_val = dixon::dixon(df, nsim = num_permutations)
      return(list(tablaZ = dixon_val$tablaZ %>%
                    dplyr::rename_with(~ .x %>% gsub(" ", "", .)) %>%
                    dplyr::mutate(!!markers[1] := df_tab[1],
                                  !!markers[2] := df_tab[2]),
                  tablaC = dixon_val$tablaC))
    })
    #bring together
    Dixon_Z = lapply(res, function(marks){
      marks$tablaZ %>%
        dplyr::mutate(!!mif$sample_id := spat[[mif$sample_id]][1], .before = 1)
    }) %>%
      do.call(dplyr::bind_rows, .) %>%
      dplyr::mutate(Simulations = num_permutations)
    Dixon_C = lapply(res, function(marks){
      marks$tablaC%>%
        tibble::rownames_to_column("Direction") %>%
        dplyr::mutate(Direction = stringr::str_squish(Direction)) %>%
        dplyr::mutate(!!mif$sample_id := spat[[mif$sample_id]][1], .before = 1)
    }) %>%
      do.call(dplyr::bind_rows, .) %>%
      dplyr::mutate(Simulations = num_permutations)
    return(list(Dixon_Z=Dixon_Z, 
                Dixon_C=Dixon_C))
  }, mc.cores = workers, mc.allow.recursive = T, mc.preschedule = F)
  
  if(overwrite){
    if("Z" %in% type){
      mif$derived$Dixon_Z = lapply(out, function(spat){
        spat$Dixon_Z %>%
          dplyr::mutate(Run = 1)
      }) %>%
        do.call(dplyr::bind_rows, .)
    }
    if("C" %in% type){
      mif$derived$Dixon_C = lapply(out, function(spat){
        spat$Dixon_C %>%
          dplyr::mutate(Run = 1)
      }) %>%
        do.call(dplyr::bind_rows, .)
    }
  } else {
    if("Z" %in% type){
      mif$derived$Dixon_Z = dplyr::bind_rows(mif$derived$Dixon_Z, 
                                             lapply(out, function(spat){
                                               spat$Dixon_Z %>%
                                                 dplyr::mutate(Run = max(Run) + 1)
                                               }) %>% 
                                               do.call(dplyr::bind_rows, .)
      )
    }
    if("C" %in% type){
      mif$derived$Dixon_C = dplyr::bind_rows(mif$derived$Dixon_C, 
                                             lapply(out, function(spat){
                                               spat$Dixon_C %>%
                                                 dplyr::mutate(Run = max(Run) + 1)
                                             }) %>% 
                                               do.call(dplyr::bind_rows, .)
      )
    }
  }
  
  structure(mif, class="mif")
}

