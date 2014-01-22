\name{plot.conc.microsim}
\alias{plot.conc.microsim}
\title{Plot Values from a Concurrency Microsimulation Model}
\usage{
\method{plot}{conc.microsim}(x, xlim, ylim, alpha = 1, ...)
}
\arguments{
  \item{x}{an \code{EpiModel} object of class
  \code{conc.microsim}.}

  \item{xlim}{x-axis scale limits for plot, with default
  calculated based on model time steps.}

  \item{ylim}{y-axis scale limits for plot, with default
  calculated based on range of data.}

  \item{alpha}{transparency level for simulation lines,
  where 0 = transparent and 1 = opaque (see
  \code{\link{transco}}).}

  \item{...}{additional arguments to pass to main plot (see
  \code{\link{plot.default}}).}
}
\description{
This function plots values from an the concurrency
microsimulation epidemic model simulated with
\code{conc.microsim}.
}
\details{
This function will extract the simulated values output
values from an \code{conc.microsim} model and plot relevant
time series data of disease prevalence and other results.
Currently, the plot function is limited to individual
simulation lines of disease prevalence; future releases
will standardize the plotting options to those available
with \code{\link{epiICM}} class models.
}

