\name{plotABM}
\alias{plotABM}
\title{Plot Simulations from a Stochastic Epidemic Model}
\usage{
  plotABM(out, compart, simlines = TRUE, simcol, simalpha,
    medline = TRUE, medcol, medlwd, qnts = NULL, qntscol,
    xlim, ylim, axs, add = FALSE, ...)
}
\arguments{
  \item{out}{list of data frames containing runs from an
  epidemic model}

  \item{compart}{compartments to plot}

  \item{simlines}{logical operator to plot individual
  simulation lines}

  \item{simcol}{color for simulation lines (default is
  black)}

  \item{simalpha}{transparency level for sim lines,
  calculated automatically if not set}

  \item{medline}{logical operator to plot median of
  individual simulations at each timestep}

  \item{medcol}{color for median line}

  \item{medlwd}{line width for median line}

  \item{qnts}{plots polygon of simulation quantiles, with
  argument taking the quantile range. For 95\% interval,
  set to 0.95.}

  \item{qntscol}{color for polygon}

  \item{xlim}{x-axis limits for plot, set to max of time if
  not set}

  \item{ylim}{y-axis limits for plot, calculated
  automatically if not set}

  \item{axs}{axis type (see \code{\link{par}} for details),
  default at 'i'}

  \item{add}{if \code{TRUE} then adds lines to existing
  plot}

  \item{...}{additional arguments to pass to plot}
}
\description{
  This function plots values of compartments in an
  agent-based stochastic epidemic model.
}
\examples{
out <- epiABM(type='SI', s.num=500, i.num=1, beta=0.2, cont=0.2, nsteps=500, nsims=20)

par(mar=c(3,3,1,1), mgp=c(2,1,0))
plotABM(out, 'i.num', simlines=TRUE, medline=TRUE, qnts=0.5)

# Plotting multiple compartment values
plotABM(out, 'i.num', simlines=FALSE, medline=TRUE, medcol=1, qnts=0.5, axs='r')
plotABM(out, 's.num', simlines=FALSE, medline=TRUE, medcol=2, qnts=0.5, qntscol=2, add=TRUE)
}
\author{
  Samuel M. Jenness <sjenness@uw.edu>
}
\seealso{
  \code{\link{epiABM}}
}
\keyword{model}

