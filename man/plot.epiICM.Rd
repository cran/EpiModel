\name{plot.epiICM}
\alias{plot.epiICM}
\title{Plot Simulations from a Individual Contact Epidemic Model}
\usage{
\method{plot}{epiICM}(x, y, popfrac, sim.lines = FALSE, sims, sim.col,
  sim.lwd, sim.alpha, mean.line = TRUE, mean.extinct = TRUE, mean.col,
  mean.lwd, mean.lty, qnts, qnts.col, qnts.alpha, leg, leg.cex, xlim, ylim,
  main, axs, add = FALSE, ...)
}
\arguments{
  \item{x}{an \code{EpiModel} object of class
  \code{epiICM}.}

  \item{y}{output compartments or flows from \code{epiICM}
  object to plot.}

  \item{popfrac}{if \code{TRUE}, plot prevalence of values
  rather than numbers (see details).}

  \item{sim.lines}{if \code{TRUE}, plot individual
  simulation lines.}

  \item{sims}{a vector representing which individual
  simulation lines to plot, with default to plot all
  lines.}

  \item{sim.col}{a vector of any standard R color format
  for simulation lines.}

  \item{sim.lwd}{line width for simulation lines.}

  \item{sim.alpha}{transparency level for simulation lines,
  where 0 = transparent and 1 = opaque (see
  \code{\link{transco}}).}

  \item{mean.line}{if \code{TRUE}, plot mean of simulations
  across time.}

  \item{mean.extinct}{if \code{TRUE}, include extinct
  simulations in mean calculation (see details).}

  \item{mean.col}{a vector of any standard R color format
  for mean lines.}

  \item{mean.lwd}{line width for mean lines.}

  \item{mean.lty}{line type for mean lines.}

  \item{qnts}{if numeric, plot polygon of simulation
  quantiles based on the range implied by the argument (see
  details). If \code{FALSE}, supress polygon plot.}

  \item{qnts.col}{a vector of any standard R color format
  for polygons.}

  \item{qnts.alpha}{transparency level for quantile
  polygons, where 0 = transparent and 1 = opaque (see
  \code{\link{transco}}).}

  \item{leg}{if \code{TRUE}, plot default legend.}

  \item{leg.cex}{a numeric value to scale the legend size.}

  \item{xlim}{x-axis scale limits for plot, with default
  calculated based on model time steps.}

  \item{ylim}{y-axis scale limits for plot, with default
  calculated based on range of data.}

  \item{main}{character string for main plot title.}

  \item{axs}{plot axis type (see \code{\link{par}} for
  details), with default to \code{"r"}.}

  \item{add}{if \code{TRUE}, new plot is not called and
  lines are added to existing plot.}

  \item{...}{additional arguments to pass to main plot (see
  \code{\link{plot.default}}).}
}
\description{
This function plots values from an individual contact
epidemic model simulated with \code{epiICM}.
}
\details{
This plotting function will extract the simulation output
values from an object of class \code{epiICM} and plot the
requested time series data of disease prevalence and other
results. The summary statistics that the function
calculates and plots are individual simulation lines, means
of the individual simulation lines, and quantiles of those
individual simulation lines. The mean line, toggled on with
\code{mean.line=TRUE} is calculated as the row mean across
simulations at each time step.

Compartment prevalences are the size of a compartment over
some denominator. To plot the raw numbers from any
compartment, use \code{popfrac=FALSE}; this is the default
for any plots of flows. The \code{popfrac} parameter
calculates and plots the denominators of all specified
compartments using these rules: 1) for one-group models,
the prevalence of any compartment is the compartment size
divided by the total population size; 2) for two-group
models, the prevalence of any compartment is the
compartment size divided by the group population size.

The quantiles show the range of outcome values within a
certain specified quantile range. By default, the
interquartile range is shown. This is specified by
\code{qnts=0.5}, that is the middle 50% of the data. If one
wanted to show the middle 95% of the date, one would
specify \code{qnts=0.95}. To toggle off the polygons where
they are plotted by default, specify \code{qnts=FALSE}.
}
\examples{
\dontrun{
# Plotting multiple compartment values automatically
mod <- epiICM(type="SIR", s.num=500, i.num=1,
             trans.rate=0.2, act.rate=0.25, rec.rate=1/50,
             nsteps=1000, nsims=10)
plot(mod)

# Plot only infecteds and show individual simulation lines
mod <- epiICM(type="SI", s.num=500, i.num=1,
              trans.rate=0.2, act.rate=0.2,
              nsteps=500, nsims=10)
plot(mod, y="i.num", sim.lines=TRUE)
}
}
\seealso{
\code{\link{epiICM}}
}
\keyword{plot}

