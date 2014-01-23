\name{as.data.frame.epiDCM}
\alias{as.data.frame.epiDCM}
\title{Extract Model Run for epiDCM Object}
\usage{
\method{as.data.frame}{epiDCM}(x, row.names = NULL, optional = FALSE,
  run = 1, ...)
}
\arguments{
  \item{x}{an \code{EpiModel} object of class
  \code{\link{epiDCM}}.}

  \item{run}{run number for model, for multiple-run
  sensitivity models run in \code{\link{epiDCM}}.}

  \item{row.names}{see
  \code{\link{as.data.frame.default}}.}

  \item{optional}{see \code{\link{as.data.frame.default}}.}

  \item{...}{see \code{\link{as.data.frame.default}}.}
}
\description{
This function extracts a model run from an object of class
\code{epiDCM} into a data frame using the generic
\code{as.data.frame} function.
}
\details{
Model output from an \code{epiDCM} run are available
through a data frame with this helper function. The output
data frame will include columns for time, the size of each
compartment, the overall population size (the sum of
compartment sizes), and the size of each flow.
}
\examples{
# One group SIS model with varying act.rate
mod <- epiDCM(type="SIS", s.num=500, i.num=1, r.num=0,
              trans.rate=0.2, act.rate=seq(0.05,0.5,0.05),
              rec.rate=1/50, nsteps=500, verbose=TRUE)
head(as.data.frame(mod, run=3))

# Two group SIR model with vital dynamics
mod <- epiDCM(type = "SIR", groups = 2,
              s.num = 500, i.num = 1, s.num.g2 = 500, i.num.g2 = 1,
              trans.rate = 0.2, trans.rate.g2 = 0.1, act.rate = 3,
              balance = "g1", rec.rate = 1/50, rec.rate.g2 = 1/50,
              b.rate = 1/100, b.rate.g2 = NA, ds.rate = 1/100,
              ds.rate.g2 = 1/100, di.rate = 1/90, di.rate.g2 = 1/90,
              dr.rate = 1/100, dr.rate.g2 = 1/100, nsteps = 500)
head(as.data.frame(mod))
}
\keyword{extract}

