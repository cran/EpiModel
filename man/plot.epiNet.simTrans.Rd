\name{plot.epiNet.simTrans}
\alias{plot.epiNet.simTrans}
\title{Plot Simulations from a Stochastic Network Epidemic Model}
\usage{
\method{plot}{epiNet.simTrans}(x, type = "sim", sim, at = 1,
  col.inf = FALSE, shp.bip = NULL, zeromarg = TRUE, alpha, lwd,
  plots.joined, ...)
}
\arguments{
  \item{x}{an \code{EpiModel} model object of class
  \code{epiNet.simTrans}.}

  \item{type}{type of plot: \code{type="sim"} for epidemic
  model results, \code{type="network"} for a static network
  plot (\code{plot.network}), or \code{type="formation"}
  for network formation statistics.}

  \item{sim}{if \code{type="network"}, simulation number
  for network graph.}

  \item{at}{if \code{type="network"}, time step for network
  graph.}

  \item{col.inf}{if \code{TRUE} and \code{type="network"},
  automatic disease status colors (blue = susceptible, red
  = infected, , green = recovered).}

  \item{shp.bip}{if \code{type="network"}, specify shapes
  in bipartite network for the second mode vertices, with
  acceptable inputs of "triangle" and "square".}

  \item{zeromarg}{if \code{TRUE} and \code{type="network"},
  automatically sets plot margins to 0 on all sides.}

  \item{alpha}{if \code{type="formation"}, transparency
  level for lines, where 0 = transparent and 1 = opaque.}

  \item{lwd}{if \code{type="formation"}, line width for
  output values.}

  \item{plots.joined}{if \code{TRUE} and
  \code{type="formation"}, combine all target statistics in
  one plot, versus one plot per target statistic if
  \code{FALSE}.}

  \item{...}{additional arguments to pass to either plot
  type.}
}
\description{
This function plots three types of output from a stochastic
network model simulated through
\code{\link{epiNet.simTrans}}.
}
\details{
This generic plot function can produce three types of plots
given a stochastic network model simulated through
\code{\link{epiNet.simTrans}}: \enumerate{ \item If
\code{type="sim"}, epidemic model results (e.g., disease
prevalence and incidence) may be plotted. In this case,
this plotting function wraps the \code{\link{plot.epiICM}}
function, as the stochastic epidemiological results are in
the same data structure. Consult the help page for
\code{plot.epiICM} for all the plotting parameters.  \item
If \code{type="network"}, a static network plot will be
generated. A static network plot of a dynamic network of
the sort simulated in \code{EpiModel} is a cross-sectional
extraction of that dynamic network at a specific time
point. This plotting function wraps the
\code{\link{plot.network}} function in the \code{network}
package. Consult the help page for \code{plot.network} for
all the plotting parameters. In addition, five plotting
parameters specific to \code{EpiModel} plots are available:
\code{sim}, \code{at}, \code{col.inf}, \code{shp.bip}, and
\code{zeromarg}.  \item If \code{type="formation"}, summary
network statistics related to partnership formation will be
plotted. These formation plots are similar to the formation
plots for both \code{\link{epiNet.est}} and
\code{\link{epiNet.simNet}} objects.  When running a
\code{epiNet.simTrans} simulation, one must specify there
that \code{save.stats=TRUE}; the plot here will then show
the network statistics requested explicitly in
\code{stat.formula}, or use the formation formula set in
\code{epiNet.est} otherwise. Note, that these network
statistics are not saved for independent
\code{epiNet.simTrans} simulations, since all of the
network data, including these statistics, is simulated and
saved in the call to \code{epiNet.simNet}; one should plot
the \code{epiNet.simNet} object in that situation. }
}
\examples{
## See EpiModel Tutorial vignette for examples ##
}
\seealso{
\code{\link{plot.epiICM}}, \code{\link{plot.network}}
}
\keyword{plot}

