\name{plot.epiDCM}
\alias{plot.epiDCM}
\title{Plot Values from a Deterministic Compartmental Epidemic Model}
\usage{
\method{plot}{epiDCM}(x, y, popfrac, run, col, lwd, lty, alpha, leg, leg.name,
  leg.cex, xlim, ylim, main, axs, add = FALSE, ...)
}
\arguments{
  \item{x}{an \code{EpiModel} object of class
  \code{\link{epiDCM}}.}

  \item{y}{output compartments or flows from \code{epiDCM}
  object to plot.}

  \item{popfrac}{if \code{TRUE}, plot prevalence of values
  rather than numbers (see details).}

  \item{run}{run number to plot for models with multiple
  runs; with a default of run 1.}

  \item{col}{color for lines, either specified as a single
  color in a standard R color format, or alternatively, a
  color palette from \code{\link{RColorBrewer}} (see
  details).}

  \item{lwd}{line width for output values.}

  \item{lty}{line type for output values.}

  \item{alpha}{transparency level for lines, where 0 =
  transparent and 1 = opaque (see \code{\link{transco}}).}

  \item{leg}{type of legend to plot. Values are "n" for no
  legend, "full" for full legend, and "lim" for limited
  legend (see details).}

  \item{leg.name}{character string to use for legend,
  calculated automatically if not supplied.}

  \item{leg.cex}{legend scale size, with default of 0.8.}

  \item{xlim}{x-axis scale limits for plot, with default
  calculated based on model time steps.}

  \item{ylim}{y-axis scale limits for plot, with default
  calculated based on range of data.}

  \item{main}{character string for main plot title.}

  \item{axs}{plot axis type (see \code{\link{par}} for
  details), with default to "r".}

  \item{add}{if \code{TRUE}, new plot is not called and
  lines are added to existing plot.}

  \item{...}{additional arguments to pass to main plot (see
  \code{\link{plot.default}}).}
}
\description{
This function plots values from an deterministic
compartment epidemic model solved with \code{epiDCM}.
}
\details{
The goal of this plotting function is to facilitate
plotting values from a deterministic compartmental model
solved with \code{\link{epiDCM}}. Depending on the number
of model runs (sensitivity analyses) and number of groups,
the default plot without any options set will generally
plot the disease prevalence from the model. The specific
compartments or flows to plot may be set using the \code{y}
parameter, and in multiple run models the specific run may
also be set.

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

Since \code{\link{epiDCM}} supports multiple run
sensitivity models, plotting the results of such models
must use a complex visual scheme for easy representation of
results. This is accomplished mainly using the
\code{\link{RColorBrewer}} color palettes, in which one may
specify a range of linked colors using named palettes. For
\code{plot.epiDCM}, one may either specify a Brewer color
palette listed in \code{\link{display.brewer.all}}, or
alternatively a single or vector of standard R colors
(named, hexidecimal, or positive integers; see
\code{\link{col2rgb}}).

There are three automatic legend types available, and the
legend is automatically added by default for most plots. To
shut off the legend, use \code{leg="n"}. To plot a legend
with values for every line, use \code{leg="full"}. With
sensitivity models with many runs, this is generally not
visually appealing Therefore, in those cases, use
\code{leg="lim"} to plot a legend limited to the highest
and lowest of the varying parameters in the model. In cases
where the default legend names are not helpful, one may
override those names with the \code{leg.names} argument.
}
\examples{
# Deterministic SIR model with varying act rate
mod <- epiDCM(type="SIR", s.num=1000, i.num=1, r.num=0,
              trans.rate=0.2, act.rate=1:10, rec.rate=1/3,
              b.rate=0.011, ds.rate=0.01, di.rate=0.03,
              dr.rate=0.01, nsteps=500, verbose=TRUE)

# Plot disease prevalence by default
plot(mod)

# Plot prevalence of susceptibles
plot(mod, y="s.num", col="Greys")

# Plot number of susceptibles
plot(mod, y="s.num", popfrac=FALSE, col="Greys")

# One way to plot multiple runs of multiple compartments together
plot(mod, y=c("s.num", "i.num"), run=5, xlim=c(0,200))
plot(mod, y=c("s.num", "i.num"), run=10, alpha=0.3, leg="n", add=TRUE)
}
\seealso{
\code{\link{epiDCM}}, \code{\link{display.brewer.all}},
\code{\link{transco}}
}
\keyword{plot}

