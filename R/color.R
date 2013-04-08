#' 
#' @title Obtain Transparent Colors
#'
#' @description This function returns a RGB transparent color from a defined hexidecimal color.
#'
#' @param hexcol hexidecimal color 
#' @param alpha transparency level (where 0 is transparent and 1 is opaque)
#' @param invis supresses printing of the RGB color
#' 
#' @author Samuel M. Jenness <sjenness@@uw.edu>
#' @keywords color
#' @export
#' 
#' @examples
#' n <- 10
#' x <- sort(sample(1:200, n))
#' y <- 10 + 2*x + rnorm(n, 0, 10)
#' z <- rpois(n, 10)
#' 
#' tcol <- transco('steelblue', 0.5)
#' 
#' plot(x,y, cex=z/4, pch=21, col='black', bg=tcol, lwd=1.2, axes=FALSE, 
#'      ylim=c(0,500), xlim=c(0,250), yaxs='r', xaxs='r')
#' axis(2, seq(0,500,100), col='white', las=2, cex.axis=0.9, mgp=c(2,0.5,0))
#' axis(1, seq(0,250,50), cex.axis=0.9, mgp=c(2,0.5,0))
#' abline(h=seq(100,500,100), col=transco('black', 0.35))
#' 
transco <- function(hexcol, alpha=1, invis=TRUE) {
  
  newa <- floor(alpha*255)
  t1 <- col2rgb(hexcol, alpha=F)
  t2 <- rep(NA, length(hexcol))
  
  for (i in 1:length(hexcol)) {
    t2[i] <- rgb(t1[1,i], t1[2,i], t1[3,i], newa, maxColorValue=255)
  }
  
  if (invis==T) {
    invisible(t2)
  } else {
    return(t2)
  }
}

#' 
#' @title RColorBrewer Color Ramp for EpiModel Plots
#'
#' @description This function returns vector of colors consistent with a high-brightness set of colors
#' from an RColorBrewer palette.
#'
#' @param plt RColorBrewer palette from \code{\link{brewer.pal}}
#' @param n number of colors to return
#' @param del.lights delete the lightest colors from the color palette
#' 
#' @details RColorBrewer provides easy access to helpful color palettes, but the built-in palettes
#' are limited to the set of colors in the existing palette. This function expands the palette size
#' to any number of colors by filling in the gaps. Also, colors within the 'div' and 'seq' set of 
#' palettes whose colors are very light (close to white) are deleted by default for better visualization
#' of plots.
#' 
#' @author Samuel M. Jenness <sjenness@@uw.edu>
#' @keywords color
#' @seealso \code{\link{RColorBrewer}}
#' @export
#' 
#' @examples
#' plot(1:100, 1:100, type='n')
#' abline(v=1:100, lwd=6, col=brewer.ramp(100, 'Spectral'))
#' abline(v=1:100, lwd=6, col=brewer.ramp(100, 'Greys'))
#' abline(v=1:100, lwd=6, col=brewer.ramp(100, 'Blues'))
#' abline(v=1:100, lwd=6, col=brewer.ramp(100, 'Set1'))
#' 
brewer.ramp <- function(n, plt, del.lights=TRUE){
  
  require(RColorBrewer)
  
  pltmax <- brewer.pal.info[row.names(brewer.pal.info)==plt,]$maxcolors
  pltcat <- brewer.pal.info[row.names(brewer.pal.info)==plt,]$category
  
  if (pltcat=='div') {
    if (del.lights) {
      colors <- brewer.pal(pltmax, plt)[-c(4:7)]
    } else {
      colors <- brewer.pal(pltmax, plt)
    }
  }
  if (pltcat=='qual') {
    colors <- brewer.pal(pltmax, plt)
  }
  if (pltcat=='seq') {
    if (del.lights) {
      colors <- rev(brewer.pal(pltmax, plt)[-c(1:3)])
    } else {
      colors <- rev(brewer.pal(pltmax, plt))
    }
  }
  
  pal <- colorRampPalette(colors)
  
  return(pal(n))
}