#' Check for presence of named columns
#' 
#' @param check.vector String(s) to check, e.g. c("peaks", "gauss.width") 
#' @param input The names() of this will be checked.
#' @param input.name Name of the input, e.g. "start.params"
#' @param function.name Name of function for which this test is carried out.
#'
#' @return Returns an error and stop signal if entries of check.vector are missing in the input data. 
#' @details example usage: check_columns(c(what, color.by, split.by, y.lab, x, y), long.df, "long.df", "plot_profile")
check_columns <- function(check.vector, input, input.name, function.name){
  
  check <- check.vector %in% names(input) 
  
  if( ! all( check )){
    stop(paste("The input (", input.name , ") supplied to ", function.name, " is missing the following:",   paste(check.vector[!check], collapse = ", ") ))
  }
}
