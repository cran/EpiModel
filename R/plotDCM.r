#' 
#' @title Plot Values from a Compartmental Epidemic Model
#'
#' @description This function plots values from an deterministic compartment epidemic model.
#'
#' @param out list of data frames containing runs from an epidemic model
#' @param compart compartments to plot
#' @param run run number for multiple compartment plots
#' @param plt palette to use from \code{\link{RColorBrewer}} palettes. The default is 'Spectral'
#' for single compartment plots and 'Set1' for multiple compartment plots.
#' @param comp.lwd line width for compartments
#' @param alpha transparency level for colors, where 0 = transparent and 1 = opaque
#' @param leg type of legend to plot. Values are 'n' for none, 'full' for full legend, 
#' and 'lim' for limited legend.
#' @param leg.name character string to use for legend, calculated automatically if not supplied
#' @param leg.cex legend size
#' @param xlim x-axis limits for plot, set to max of time if not set
#' @param ylim y-axis limits for plot, calculated automatically if not set
#' @param axs axis type (see \code{\link{par}} for details), default at 'i'
#' @param add if \code{TRUE} then adds lines to existing plot
#' @param ... additional arguments to pass to main plot
#' 
#' @author Samuel M. Jenness <sjenness@@uw.edu>
#' @keywords model
#' @seealso \code{\link{epiDCM}}, \code{\link{display.brewer.all}}
#' @export
#' 
#' @examples
#' out <- epiDCM(type='SIR', s.num=1000, i.num=1, r.num=0, 
#'               beta=0.2, cont=1:10, nu=1/3,
#'               b=0.011, ms=0.01, mi=0.03, mr=0.01, 
#'               dt=500, verbose=TRUE)
#'                
#' par(mar=c(3.5,3,1,1), mgp=c(2,1,0))
#' plotDCM(out, compart='s.prev', plt='Greys')
#' plotDCM(out, compart='i.num', leg='lim', alpha=0.5, xlim=c(0,200), axs='r')
#' 
#' plotDCM(out, compart=c('s.num','r.num'), run=5, leg='full')
#' plotDCM(out, compart=c('s.num','r.num'), run=2, alpha=0.3, add=TRUE)
#' 
plotDCM <- function(out, compart, run, plt, comp.lwd, alpha,
                      leg, leg.name, leg.cex,  
                      xlim, ylim, axs, add=FALSE, ...) {
  
  require(RColorBrewer)
  
  # model dimensions
  nts <- max(out$time)
  nvars <- length(out)
  if (class(out$s.num)=='numeric') nruns <- 1
  if (class(out$s.num)=='data.frame') nruns <- dim(out$s.num)[2]
  lcomp <- length(compart)
  
  # find compartment ymax across runs
  allmax <- vector()
  for (i in 1:length(compart)) allmax[i] <- max(out[[compart[i]]])
  ymax <- ceiling(max(allmax))
  
  # defaults for missing arguments
  if (missing(ylim)) ylim=c(0, ymax)
  if (missing(xlim)) xlim=c(0, nts)
  if (missing(axs)) axs='i'
  if (missing(comp.lwd)) comp.lwd=2.5
  if (missing(leg)) leg='n'
  if (missing(alpha)) alpha=0.9
  if (missing(leg.cex)) leg.cex=0.9
  
  # Plot
  if (add==FALSE) {
    Time <- Out <- 1
    plot(Time, Out, type='n', bty='n', 
         xaxs=axs, yaxs=axs, xlim=xlim, ylim=ylim, ...)
  }
  
  if (lcomp == 1) {
    # Lines
    if (missing(plt)) plt <- 'Spectral'
    pal <- transco(brewer.ramp(nruns, plt), alpha)
  	if (nruns==1) {
  	  lines(out$time, out[[compart]], lwd=comp.lwd, col=pal[i])
  	} else {
      for (i in 1:nruns) {
    		lines(out$time, out[[compart]][,i], lwd=comp.lwd, col=pal[i])
    	}
  	}
    
    # Legend names
    if (missing(leg.name)){
      if (nruns == 1) leg.names <- compart
      if (nruns > 1) leg.names <- names(out[[compart[1]]])
    } else {
      leg.names <- paste(leg.name, 1:nruns, sep=' ')
    }
    
    # Legend
  	if (leg == 'n') {
  	} else if (leg == 'full') {
    	legend('topright', legend=leg.names, 
             bg='white', lty=1, lwd=comp.lwd,
    		     col=pal[1:nruns], cex=leg.cex)
  	} else if (leg == 'lim') {
    	legend('topright', legend=c(leg.names[1], leg.names[nruns]), 
    		     bg='white', lty=1, lwd=comp.lwd,
    		     col=pal[c(1, nruns)], cex=leg.cex) 
  	}
    invisible(pal)
  }
  
  if (lcomp > 1) {
    
    # Lines
    if (missing(plt)) plt <- 'Set1'
    if (lcomp > 3) {
      pal <- transco(brewer.ramp(lcomp, plt), alpha)
    } else {
      pal <- transco(brewer.pal(max(c(lcomp,3)), plt), alpha)
    }
    
    if (nruns==1) {
      for (i in 1:length(compart)) {
        lines(out$time, out[[compart[i]]], 
              lwd=comp.lwd, col=pal[i])
      }
    } 
    if (nruns>1) {
      if (missing(run)) {
        for (j in 1:length(compart)) {
            lines(out$time, out[[compart[j]]][,1], 
                  lwd=comp.lwd, col=pal[j])
        }
      } else {
        for (j in 1:length(compart)) {
            lines(out$time, out[[compart[j]]][,run], 
                  lwd=comp.lwd, col=pal[j])
        } 
      }
    }

    
    # Legend
    if (leg == 'n') {  }
    else if (leg == 'full') {
      legend('topright', legend=compart, bg='white', 
             lty=1, col=pal, cex=leg.cex, lwd=comp.lwd) 
    }
    invisible(pal) 
  }
   
}