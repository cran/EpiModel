\name{plot.epiNet.est}
\alias{plot.epiNet.est}
\title{Plot Diagnostics from an epiNet.est Object}
\usage{
\method{plot}{epiNet.est}(x, type = "formation", dx.start, dx.end,
  dx.leg = TRUE, plots.joined, ...)
}
\arguments{
  \item{x}{an \code{EpiModel} object of class
  \code{\link{epiNet.est}}.}

  \item{type}{plot type, with options of
  \code{type="formation"} for partnership formation
  statistics or other network statistics specified in
  \code{epiNet.est}, or \code{type="duration"} to plot the
  mean ages of partnerships over time.}

  \item{dx.start}{start time for diagnostic plots. This
  must be a positive integer.}

  \item{dx.end}{end time for diagnostic plots. This must be
  less than or equal to the number of time steps simulated
  in the \code{epiNet.est} diagnostics.}

  \item{dx.leg}{if \code{TRUE}, show legend (only if
  \code{plots.joined=TRUE})}

  \item{plots.joined}{if \code{TRUE}, combine all target
  statistics in one plot, versus one plot per target
  statistic if \code{FALSE}}

  \item{...}{additional arguments to pass to main plot
  window}
}
\description{
This function plots values from diagnostic simulations in
\code{epiNet.est}.
}
\details{
The plot function for \code{epiNet.est} objects will
generate plots of two types of model diagnostic statistics
that were run as part of the diagnostic tools within that
estimation function. The \code{formation} plot shows the
summary statistics requested in \code{stats.form}, which
defaults to those in the formation formula. The
\code{duration} plot shows the average age of all
partnerships at each time step, up until the maximum time
step requested. This is calculated with the
\code{\link{edgelist.meanage}} function.

The \code{plots.joined} argument will control whether the
target statistics in the \code{formation} plot are joined
in one plot or plotted separately. The default is based on
the number of network statistics requested. The layout of
the separate plots within the larger plot window is also
based on the number of statistics.

Required for these plots is that the estimation diagnostics
are run in \code{epiNet.est}. This happens by default, and
is set with the \code{stats=TRUE} argument in that
function. Since these diagnostics only generate one network
simulation, the plots here will only show one line per
network statistic. If the range of possible statistics is
of interest, then the fitted STERGM should be simulated
multiple times, which is possible with
\code{\link{epiNet.simNet}}, which has its own plotting
function, \code{\link{plot.epiNet.simNet}}.
}
\examples{
\dontrun{
# Initialize bipartite network
nw <- network.initialize(500, bipartite=250, directed=FALSE)

# Fit model with limited degree terms, but monitor wider range
est <- epiNet.est(nw,
                  formation = ~ edges + b1degree(0:1) + b2degree(0:1),
                  dissolution = ~ offset(edges),
                  target.stats = c(165, 100, 137.5, 120, 102.5),
                  coef.diss = dissolution.coefs(~ offset(edges), duration=25),
                  stats.form = ~ edges + b1degree(0:5) + b2degree(0:5))

# Formation plot with plots split given the 11 statistics
plot(est)

# Formation plot with plots joined, for time steps 400 to 500
plot(est, plots.joined=TRUE, dx.start=400, dx.end=500)

# Partnership duration plot
plot(est, type="duration")

# Truncate plot to start after time 200, given the age ramp up period
plot(est, type="duration", dx.start=200)
}
}
\seealso{
\code{\link{epiNet.est}}
}
\keyword{plot}

