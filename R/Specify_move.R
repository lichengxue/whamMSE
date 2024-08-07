#' Launch the Shiny App
#'
#' This function launches the Shiny app included in the package.
#' @export
Specify_move <- function() {
  appDir <- system.file("shinyApp/app.R", package = "whamMSE")
  if (appDir == "") {
    stop("Could not find the Shiny app directory. Try re-installing `whamMSE`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}

