#' 
#' @title Plot Simulations from a Stochastic Epidemic Model
#'
#' @description This function plots values of compartments in an agent-based
#'  stochastic epidemic model.
#'
#' @param out list of data frames containing runs from an epidemic model
#' @param compart compartments to plot
#' @param simlines logical operator to plot individual simulation lines
#' @param simcol color for simulation lines (default is black)
#' @param simalpha transparency level for sim lines, calculated automatically if not set
#' @param medline logical operator to plot median of individual simulations at each timestep
#' @param medcol color for median line
#' @param medlwd line width for median line
#' @param qnts plots polygon of simulation quantiles, with argument taking the quantile range. 
#' For 95\% interval, set to 0.95.
#' @param qntscol color for polygon
#' @param xlim x-axis limits for plot, set to max of time if not set
#' @param ylim y-axis limits for plot, calculated automatically if not set
#' @param axs axis type (see \code{\link{par}} for details), default at 'i'
#' @param add if \code{TRUE} then adds lines to existing plot
#' @param ... additional arguments to pass to plot
#' 
#' @author Samuel M. Jenness <sjenness@@uw.edu>
#' @keywords model
#' @seealso \code{\link{epiABM}}
#' @export
#' 
#' @examples
#' out <- epiABM(type='SI', s.num=500, i.num=1, beta=0.2, cont=0.2, nsteps=500, nsims=20)
#' 
#' par(mar=c(3,3,1,1), mgp=c(2,1,0))
#' plotABM(out, 'i.num', simlines=TRUE, medline=TRUE, qnts=0.5)
#' 
#' # Plotting multiple compartment values
#' plotABM(out, 'i.num', simlines=FALSE, medline=TRUE, medcol=1, qnts=0.5, axs='r')
#' plotABM(out, 's.num', simlines=FALSE, medline=TRUE, medcol=2, qnts=0.5, qntscol=2, add=TRUE)
#' 
plotABM <- function(out, compart, 
                    simlines=TRUE, simcol, simalpha,
                    medline=TRUE, medcol, medlwd,
                    qnts=NULL, qntscol, 
                    xlim, ylim, axs, add=FALSE, ...) {
    
  # model dimensions
  nts <- max(out$time)
  nsims <- ncol(out$s.num)
 
  if (length(compart)>1) stop('Multiple compartments possible with second call with add=TRUE')
  
  # Palettes
  require(RColorBrewer)
  if (missing(simcol)) simcol <- 'black'
  if (missing(simalpha)) simalpha <- max(c(0.05, 1-log10(nsims)/3))
    simpal <- transco(simcol, simalpha)
  bpal <- brewer.pal(3, 'Set1')
    if (missing(medcol)) medcol <- bpal[1] 
    medpal <- transco(medcol, 0.9)
  if (missing(qntscol)) qntscol <- bpal[2]
    qntspal <- transco(qntscol, 0.4)
  
  # Compartment summaries
  max.prev <- max(out[[compart]])
  min.prev <- min(out[[compart]])
  
  # Missing args
  if (missing(xlim)) xlim <- c(0,nts)
  if (missing(ylim)) ylim <- c(0, max.prev)
  
  # Main plot
  if (missing(axs)) axs <- 'i'
  Time <- Out <- 1
  if (add==FALSE) {
    plot(Time, Out, type='n', xlim=xlim, ylim=ylim,
         bty='n', xaxs=axs, yaxs=axs, ...)
  }
  
  # Quantiles
  if (!is.null(qnts)) {
    quants <- c((1-qnts)/2, 1-((1-qnts)/2))
    qnt.prev <- apply(out[[compart]], 1, function(x) quantile(x, c(quants[1], quants[2])))
    xx <- c(1:nts, nts:1)
    yy <- c(qnt.prev[1,], rev(qnt.prev[2,]))
    polygon(xx, yy, col=qntspal, border=NA)
  }
  
  # Individual sim lines
  if (simlines) {
    for (i in 1:nsims) {
      lines(1:nts, out[[compart]][,i], lwd=0.8, col=simpal)
    }
  }
  
  # Median of sim lines
  if (missing(medlwd)) medlwd <- 2.5
  if (medline) {
    med.prev <- apply(out[[compart]], 1, median)
    lines(1:nts, med.prev, lwd=medlwd, col=medpal)
  }
  
}
