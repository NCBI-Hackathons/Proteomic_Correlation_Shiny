#' Start the app
#' 
#' This will start the shiny app.
#' 
#' @export
runPCP <- function() {
  aDir <- system.file("shiny", package = "DepLab")
  if (aDir == "") {
    stop("Could not find shiny directory. Try re-installing `DepLab`.", call. = FALSE)
  }
  
  shiny::runApp(appDir = aDir, display.mode = "normal")
}
