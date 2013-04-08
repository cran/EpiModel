\name{statsDCM}
\alias{statsDCM}
\title{Extract Model Statistics}
\usage{
  statsDCM(out, time, run = 1, digits = 3, console = TRUE)
}
\arguments{
  \item{out}{array containing model results}

  \item{time}{timestep for model statistics}

  \item{run}{run number for model if sensitivity analyses
  are conducted}

  \item{digits}{number of significant digits to print}

  \item{console}{print to R console}
}
\description{
  This function extracts and prints model statistics solved
  with the \code{mod.} class functions.
}
\examples{
out <- epiDCM(type='SIR', s.num=1000, i.num=1, r.num=0,
                beta=0.2, cont=2:8, nu=1/3,
                b=0.011, ms=0.01, mi=0.03, mr=0.01,
                dt=500, verbose=TRUE)
statsDCM(out, time=25, run=1)
statsDCM(out, time=25, run=7)
statsDCM(out, time=26, run=7)
}
\author{
  Samuel M. Jenness <sjenness@uw.edu>
}
\keyword{model}

