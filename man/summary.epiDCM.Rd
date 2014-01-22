\name{summary.epiDCM}
\alias{summary.epiDCM}
\title{Summary Model Statistics}
\usage{
\method{summary}{epiDCM}(object, time, run = 1, digits = 3,
  comp.plot = FALSE, ...)
}
\arguments{
  \item{object}{an \code{EpiModel} object of class
  \code{epiDCM}.}

  \item{time}{time step for model statistics.}

  \item{run}{run number for model, if sensitivity analyses
  were conducted through \code{epiDCM}.}

  \item{digits}{number of significant digits to print.}

  \item{comp.plot}{if \code{TRUE}, plot compartments and
  flows in summary. This can be called separately using the
  \code{\link{comp.plot}} function.}

  \item{...}{additional summary function arguments.}
}
\description{
This function extracts and prints model statistics solved
with the \code{epiDCM} function.
}
\details{
Summary statistics for the main epidemiological outcomes
(state and transition size and prevalence) from an
\code{epiDCM} model. Time-specific summary measures are
provided, so it is necessary to input a time of interest.
For multiple-run models (e.g., sensitivity analyses), one
may also input a model run number. See the examples below.
}
\examples{
\dontrun{
# Deterministic SIR model with varying act.rate
mod <- epiDCM(type="SIR", s.num=1000, i.num=1, r.num=0,
              trans.rate=0.2, act.rate=2:8, rec.rate=1/3,
              b.rate=0.011, ds.rate=0.01, di.rate=0.03,
              dr.rate=0.01, nsteps=500, verbose=TRUE)
summary(mod, time=25, run=1, comp.plot=TRUE)
summary(mod, time=25, run=7, comp.plot=TRUE)
summary(mod, time=26, run=7, comp.plot=TRUE)
}
}
\seealso{
\code{\link{epiDCM}}
}
\keyword{extract}

