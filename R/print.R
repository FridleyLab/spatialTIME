#' @export
print.mif <- function(x, ...){
  
  # create color schemes for output text 
  emphesis <- crayon::make_style("deepskyblue")
  
  # x$sample <- x$sample %>% 
  #   janitor::clean_names()
  
  #NULL[["id"]] is an ERROR in R, not NULL, so a partially built mif used to be
  #impossible even to echo at the console -- auto-printing threw
  #"subscript out of bounds".
  n_unique <- function(tbl, col) {
    if (is.null(tbl) || is.null(col) || !col %in% colnames(tbl)) return(NA_integer_)
    length(unique(tbl[[col]]))
  }
  fmt <- function(n) if (is.na(n)) emphesis("?") else emphesis(n)

  cat(fmt(n_unique(x$clinical, x$patient_id)), "patients spanning",
      fmt(n_unique(x$sample, x$sample_id)), "samples and",
      emphesis(length(x$spatial)), "spatial data frames were found \n")
  if (length(x$derived)) {
    cat("  derived:", paste(names(x$derived), collapse = ", "), "\n")
  }
  #Print methods must return their argument invisibly, or `y <- print(mif)` yields
  #NULL and print() cannot sit in a pipe.
  invisible(x)
}