\name{plotDCM}
\alias{plotDCM}
\title{Plot Values from a Compartmental Epidemic Model}
\usage{
  plotDCM(out, compart, run, plt, comp.lwd, alpha, leg,
    leg.name, leg.cex, xlim, ylim, axs, add = FALSE, ...)
}
\arguments{
  \item{out}{list of data frames containing runs from an
  epidemic model}

  \item{compart}{compartments to plot}

  \item{run}{run number for multiple compartment plots}

  \item{plt}{palette to use from \code{\link{RColorBrewer}}
  palettes. The default is 'Spectral' for single
  compartment plots and 'Set1' for multiple compartment
  plots.}

  \item{comp.lwd}{line width for compartments}

  \item{alpha}{transparency level for colors, where 0 =
  transparent and 1 = opaque}

  \item{leg}{type of legend to plot. Values are 'n' for
  none, 'full' for full legend, and 'lim' for limited
  legend.}

  \item{leg.name}{character string to use for legend,
  calculated automatically if not supplied}

  \item{leg.cex}{legend size}

  \item{xlim}{x-axis limits for plot, set to max of time if
  not set}

  \item{ylim}{y-axis limits for plot, calculated
  automatically if not set}

  \item{axs}{axis type (see \code{\link{par}} for details),
  default at 'i'}

  \item{add}{if \code{TRUE} then adds lines to existing
  plot}

  \item{...}{additional arguments to pass to main plot}
}
\description{
  This function plots values from an deterministic
  compartment epidemic model.
}
\examples{
out <- epiDCM(type='SIR', s.num=1000, i.num=1, r.num=0,
              beta=0.2, cont=1:10, nu=1/3,
              b=0.011, ms=0.01, mi=0.03, mr=0.01,
              dt=500, verbose=TRUE)

par(mar=c(3.5,3,1,1), mgp=c(2,1,0))
plotDCM(out, compart='s.prev', plt='Greys')
plotDCM(out, compart='i.num', leg='lim', alpha=0.5, xlim=c(0,200), axs='r')

plotDCM(out, compart=c('s.num','r.num'), run=5, leg='full')
plotDCM(out, compart=c('s.num','r.num'), run=2, alpha=0.3, add=TRUE)
}
\author{
  Samuel M. Jenness <sjenness@uw.edu>
}
\seealso{
  \code{\link{epiDCM}}, \code{\link{display.brewer.all}}
}
\keyword{model}

