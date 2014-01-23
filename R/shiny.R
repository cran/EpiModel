
#' @title GUI for epiDCM Models
#'
#' @description This function runs a web browser-based GUI of the epiDCM
#'   application using the \code{shiny} package.
#' 
#' @details
#' Using the \code{shiny} interface, this function runs a web-based GUI of
#'   a one-group \code{epiDCM} model with user input on model type, state
#'   sizes, and parameters. Model output may be plotted, summarized, and
#'   saved as raw data using the core \code{EpiModel} functionality for 
#'   \code{epiDCM} models. 
#' 
#' @seealso \code{\link{epiDCM}}, \code{\link{plot.epiDCM}}, 
#'   \code{\link{summary.epiDCM}}, \code{\link{as.data.frame.epiDCM}}
#' 
#' @keywords GUI
#' @export
#' 
#' @examples
#' \dontrun{
#'  gui.epiDCM()
#' }
#' 
gui.epiDCM <- function() {

  shiny::runApp(system.file('shiny', 'epiDCM', package = 'EpiModel'))

}
