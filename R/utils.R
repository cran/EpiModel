#' 
#' @title Obtain Transparent Colors
#'
#' @description This function returns an RGB transparent color from an 
#'  input R color
#'
#' @param col vector of any of the three kinds of \code{R} color specifications
#'   (named, hexadecimal, or positive integer; see \code{\link{col2rgb}})
#' @param alpha transparency level, where 0 is transparent and 1 is opaque
#' @param invisible supresses printing of the RGB color
#' 
#' @details
#' The purpose of this function is to facilitate color transparency, which is
#' used widely in \code{EpiModel} plots. This is an internal function that is
#' not ordinarily called by the end-user.
#' 
#' @return
#' A vector of length equal to the input \code{col} vector, containing the 
#' transformed color values in hexidemical format.
#' 
#' @seealso \code{\link{rgb}}, \code{\link{col2rgb}}
#' 
#' @export
#' @keywords colorUtils internal
#' 
#' @examples
#' # Three-variable bubble plot to show transparency
#' n <- 25
#' x <- sort(sample(1:200, n))
#' y <- 10 + 2*x + rnorm(n, 0, 10)
#' z <- rpois(n, 10)
#' 
#' # Input named R color and transparency level
#' tcol <- transco("steelblue", 0.5)
#' 
#' plot(x,y, cex=z/4, pch=21, col="black", bg=tcol, lwd=1.2, axes=FALSE, 
#'      ylim=c(0,500), xlim=c(0,250), yaxs="r", xaxs="r")
#' axis(2, seq(0,500,100), col="white", las=2, cex.axis=0.9, mgp=c(2,0.5,0))
#' axis(1, seq(0,250,50), cex.axis=0.9, mgp=c(2,0.5,0))
#' abline(h=seq(100,500,100), col=transco("black", 0.35))
#' 
transco <- function(col, 
                    alpha=1, 
                    invisible=TRUE
                    ) {
  
  if (alpha > 1 || alpha < 0)
    stop("Specify alpha between 0 and 1")
  
  newa <- floor(alpha*255)
  t1 <- col2rgb(col, alpha = FALSE)
  t2 <- rep(NA, length(col))
  
  for (i in 1:length(col)) {
    t2[i] <- rgb(t1[1,i], t1[2,i], t1[3,i], newa, maxColorValue=255)
  }
  
  if (invisible == TRUE) {
    invisible(t2)
  } else {
    return(t2)
  }
}


#' @title RColorBrewer Color Ramp for EpiModel Plots
#'
#' @description This function returns vector of colors consistent with 
#'   a high-brightness set of colors from an \code{RColorBrewer} palette.
#'
#' @param plt \code{RColorBrewer} palette from \code{\link{brewer.pal}}
#' @param n number of colors to return
#' @param delete.lights delete the lightest colors from the color palette,
#'   helps with plotting in many high-contrast palettes
#' 
#' @details 
#' \code{RColorBrewer} provides easy access to helpful color palettes, but the 
#' built-in palettes are limited to the set of colors in the existing palette. 
#' This function expands the palette size to any number of colors by filling 
#' in the gaps. Also, colors within the "div" and "seq" set of palettes whose 
#' colors are very light (close to white) are deleted by default for better 
#' visualization of plots.
#' 
#' @return
#' A vector of length equal to \code{n} with a range of color values consistent
#' with an RColorBrewer color palette.
#' 
#' @seealso \code{\link{RColorBrewer}}
#' @keywords colorUtils internal
#' @export
#' 
#' @examples
#' # Shows a 100-color ramp for 4 RColorBrewer palettes
#' par(mfrow=c(2,2), mar=c(1,1,2,1))
#' pals <- c("Spectral", "Greys", "Blues", "Set1")
#' for (i in 1:length(pals)) {
#'  plot(1:100, 1:100, type="n", axes=FALSE, main=pals[i])
#'  abline(v=1:100, lwd=6, col=brewer.ramp(100, pals[i]))
#' }
#' 
brewer.ramp <- function(n, plt, delete.lights=TRUE){
    
  pltmax <- brewer.pal.info[row.names(brewer.pal.info)==plt,]$maxcolors
  pltcat <- brewer.pal.info[row.names(brewer.pal.info)==plt,]$category
  
  if (pltcat == "div") {
    if (delete.lights == TRUE) {
      colors <- brewer.pal(pltmax, plt)[-c(4:7)]
    } else {
      colors <- brewer.pal(pltmax, plt)
    }
  }
  if (pltcat == "qual") {
    colors <- brewer.pal(pltmax, plt)
  }
  if (pltcat == "seq") {
    if (delete.lights == TRUE) {
      colors <- rev(brewer.pal(pltmax, plt)[-c(1:3)])
    } else {
      colors <- rev(brewer.pal(pltmax, plt))
    }
  }
  
  pal <- colorRampPalette(colors)
  
  return(pal(n))
}


#' @title Creates a TEA Variable for Infection Status for \code{ndtv} Animations
#'
#' @description This function creates a new temporally-extended attribute (TEA) 
#'   variable in a \code{networkDynamic} object with color name output from an 
#'   existing status variable with numeric format.
#' 
#' @param nw an object of class \code{networkDynamic}.
#' @param old.var old TEA variable name.
#' @param old.sus status value for susceptible in old TEA variable.
#' @param old.inf status value for infected in old TEA variable.
#' @param old.rec status value for recovered in old TEA variable.
#' @param new.var new TEA variable name to be stored in nw object.
#' @param new.sus status value for susceptible in new TEA variable.
#' @param new.inf status value for infected in new TEA variable.
#' @param new.rec status value for recovered in new TEA variable.
#' @param verbose print progress for calculations.
#' 
#' @details
#' The \code{ndtv} package allows for animated plots of dynamic network objects, 
#' showing the evolving partnership structure. The \code{EpiModel} package uses 
#' temporally-extended attributes (TEAs) to store longitudinal disease status for 
#' every vertex in a dynamic network. To visualize disease status dynamically in 
#' \code{ndtv}, it is currently necessary to create a TEA containing the colors 
#' to be used in drawing the nodes. 
#'  
#' The convention in \code{\link{plot.epiNet.simTrans}} is to color the susceptible 
#' nodes as blue, infected nodes as red,  and recovered nodes as green. This 
#' function allows the user to take the existing status TEA on the 
#' \code{networkDynamic} object and transform that numeric variable into any colors. 
#' 
#' @seealso \code{\link{epiNet.simTrans}} and the \code{ndtv} package documentation.
#' @keywords colorUtils
#' @export
#'
#' @examples
#' ## See EpiModel Tutorial vignette ##
#' 
colorTEA <- function(nw, 
                     old.var = "status", 
                     old.sus = 0, 
                     old.inf = 1, 
                     old.rec = 2,
                     new.var = "ndtvcol", 
                     new.sus, 
                     new.inf, 
                     new.rec,
                     verbose = TRUE
                     ) {
  
  if (missing(new.inf)) new.inf <- transco("firebrick", 0.75)
  if (missing(new.sus)) new.sus <- transco("steelblue", 0.75)
  if (missing(new.rec)) new.rec <- transco("seagreen", 0.75)
  
  times <- 1:max(get.change.times(nw))
  
  for (ts in times) {
    
    stat <- get.vertex.attribute.active(nw, old.var, at=ts)
    infected <- which(stat == old.inf)
    uninfected <- which(stat == old.sus)
    recovered <- which(stat == old.rec)
    
    activate.vertex.attribute(nw, prefix=new.var, value=new.inf,
                              onset=ts, terminus=Inf, v=infected)
    activate.vertex.attribute(nw, prefix=new.var, value=new.sus, 
                              onset=ts, terminus=Inf, v=uninfected)
    activate.vertex.attribute(nw, prefix=new.var, value=new.rec,
                              onset=ts, terminus=Inf, v=recovered)
    
    if (verbose == TRUE) cat(ts, "/", max(times), "\t", sep="")
  }
  
  return(nw)
}
