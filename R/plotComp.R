#' 
#' @title Plot Compartment Diagram
#'
#' @description This function plots the compartment diagram for standard models included in EpiModel
#'
#' @param type Type of deterministic model to plot. 
#' 
#' @details Current options only are: 
#' \itemize{
#'   \item \code{'SIR'} simple SIR models
#'   \item \code{'SIRvd'} SIR models with vital dynamics
#'  }
#' 
#' @author Samuel M. Jenness <sjenness@@uw.edu>
#' @keywords model
#' @seealso \code{\link{epiDCM}}
#' @export
#' 
#' @examples
#' plotComp(type='SIR')
#' plotComp(type='SIRvd')
plotComp <- function(type) {

  if (type %in% c('SIR', 'SIRvd')) {
    par(mar=c(0,0,0,0))
    plot(0:100, 0:100, type='n', axes=F)
    text(15, 50, 'Susceptible')
      polygon(c(5, 25, 25, 5), c(40, 40, 60, 60))
    text(50, 50, 'Infected')
      polygon(c(40, 60, 60, 40), c(40, 40, 60, 60))
    text(85, 50, 'Recovered')
      polygon(c(75, 95, 95, 75), c(40, 40, 60, 60))
    text(32.5, 55, expression(lambda), cex=1.2)
      arrows(25, 50, 40, lwd=2, length=0.15)
    text(67.5, 55, expression(nu), cex=1.2)
      arrows(60, 50, 75, lwd=2, length=0.15)
    
    if (type=='SIRvd') {
      text(17.5, 27.5, expression(italic(m[s])))
        arrows(15, 40, 15, 15, lwd=2, length=0.15)
      text(52.5, 27.5, expression(italic(m[i])))
        arrows(50, 40, 50, 15, lwd=2, length=0.15)  
      text(87.5, 27.5, expression(italic(m[r])))
        arrows(85, 40, 85, 15, lwd=2, length=0.15)
      text(17.5, 72.5, expression(italic(b)))
        arrows(15, 85, 15, 60, lwd=2, length=0.15)
    }
  }
  
}
