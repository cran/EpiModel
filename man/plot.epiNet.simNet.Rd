\name{plot.epiNet.simNet}
\alias{plot.epiNet.simNet}
\title{Plot Diagnostics from an epiNet.simNet Object}
\usage{
\method{plot}{epiNet.simNet}(x, type = "formation", sim, dx.start, dx.end,
  dx.leg = TRUE, plots.joined, alpha, lwd, ...)
}
\arguments{
  \item{x}{an \code{EpiModel} object of class
  \code{\link{epiNet.simNet}}.}

  \item{type}{plot type, with options either of
  \code{type="formation"} for partnership formation
  statistics or other network statistics specified in
  \code{epiNet.simNet}, or \code{type="duration"} to plot
  the mean ages of partnerships over time.}

  \item{sim}{network simulation number to plot, with
  default to plot all simulations in \code{epiNet.simNet}
  object.}

  \item{dx.start}{start time for diagnostic plots. This
  must be a positive integer.}

  \item{dx.end}{end time for diagnostic plots. This must be
  less than or equal to the number of time steps simulated
  in the \code{epiNet.est} diagnostics.}

  \item{dx.leg}{if \code{TRUE}, show legend (only if
  \code{plots.joined=TRUE})}

  \item{plots.joined}{if \code{TRUE}, combine all target
  statistics in one plot, versus one plot per target
  statistic if \code{FALSE}.}

  \item{alpha}{transparency level for lines, where 0 =
  transparent and 1 = opaque.}

  \item{lwd}{line width for output values.}

  \item{...}{additional arguments to pass to main plot.}
}
\description{
This function plots values from diagnostic simulations in
\code{epiNet.simNet}.
}
\details{
The plot function for \code{epiNet.simNet} objects will
generate plots of two types of model diagnostic statistics
that were run as part of the diagnostic tools within that
simulation function. The \code{formation} plot shows the
summary statistics requested in \code{stats.form}, which
defaults to those in the formation formula. The
\code{duration} plot shows the average age of all
partnerships at each time step, up until the maximum time
step requested. This is estimated with the
\code{\link{edgelist.meanage}} function.

The \code{plots.joined} argument will control whether the
target statistics in the \code{formation} plot are joined
in one plot or plotted separately. The default is based on
the number of network statistics requested. The layout of
the separate plots within the larger plot window is also
based on the number of statistics.

Required for these plots is that the estimation diagnostics
are run in \code{epiNet.simNet}. This happens by default,
and is set with the \code{stats=TRUE} argument in that
function.
}
\examples{
\dontrun{
# See epiNet.est help to understand these steps
nw <- network.initialize(n = 500, directed = FALSE)
nw <- set.vertex.attribute(nw, "race", value = rep(0:1, each = 250))
formation <- ~ edges + nodematch("race") + degree(0) + concurrent
target.stats <- c(225, 187, 180, 90)
dissolution <- ~ offset(edges)
coef.diss <- dissolution.coefs(dissolution, durations = 20)

# Model estimation using the direct STERGM fit, then simulation
est <- epiNet.est(
  nw,
  formation,
  dissolution,
  target.stats,
  coef.diss,
  edapprox = FALSE)
nwsims <- epiNet.simNet(est, nsteps = 250, nsims = 10)

# Default will plot all simulations for target statistics formation plot
plot(nwsims)

# If that is too busy, try adjusting the line width or transparency with alpha
plot(nwsims, lwd=0.5)
plot(nwsims, lwd=0.75, alpha=0.5)

# If still too busy, then split plots
plot(nwsims, plots.joined=FALSE)

# Or plot just a subset of simulations
plot(nwsims, sim=1)
plot(nwsims, sim=4:6, lwd=0.6)

# Or plot just a subset of time
plot(nwsims, dx.start=100, dx.end=150, alpha=0.4)

# Duration plot shows the average age of partnerships along each time step
# Since a nodematch term was used in the dissolution, it will plot
#   partnerships matched and unmatched separately
plot(nwsims, type="duration", alpha=0.5)

# Truncate plot to start after time 100, given the age ramp up period
plot(nwsims, type="duration", dx.start=100, lwd=0.5)
}
}
\seealso{
\code{\link{epiNet.simNet}}
}
\keyword{plot}

