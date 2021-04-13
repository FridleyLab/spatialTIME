#' @export
print.mif <- function(x, ...){
  
  # create color schemes for output text 
  emphesis <- crayon::make_style("deepskyblue")
  
  # x$sample <- x$sample %>% 
  #   janitor::clean_names()
  
  cat(emphesis(length(unique(x$clinical[[1]]))), "patients spanning",
      emphesis(length(unique(x$sample[[2]]))), "samples and",
      emphesis(length(x$spatial)), "spatial data frames were found \n")
  # add number/names of genes in a panel
}