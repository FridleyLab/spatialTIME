#' @export
print.mif <- function(x, ...){
  
  # create color schemes for output text 
  emphesis <- crayon::make_style("deepskyblue")
  
  # x$sample <- x$sample %>% 
  #   janitor::clean_names()
  
  cat(emphesis(length(unique(x$clinical[[x$patient_id]]))), "patients spanning",
      emphesis(length(unique(x$sample[[x$sample_id]]))), "samples and",
      emphesis(length(x$spatial)), "spatial data frames were found \n")
  # add number/names of genes in a panel
}