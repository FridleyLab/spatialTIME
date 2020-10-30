#' @export
print.mif <- function(x, ...){
  
  # create color schemes for output text 
  emphesis <- crayon::make_style("deepskyblue")
  
  cat(emphesis(length(unique(x$clinical$patient_id))), "patients spanning",
      emphesis(length(unique(x$sample$image.tag))), "samples and",
      emphesis(length(x$spatial)), "spatial data frames were found \n")
  # add number/names of genes in a panel
}