\name{getDCM}
\alias{getDCM}
\title{Extract Model Run}
\usage{
  getDCM(out, run = 1)
}
\arguments{
  \item{out}{array containing model results}

  \item{run}{run number for model if sensitivity analyses
  are conducted}
}
\description{
  This function extracts a specific model run for epiDCM
}
\examples{
out <- epiDCM(type='SIR', s.num=1000, i.num=1, r.num=0,
               beta=0.2, cont=1:4, nu=1/3,
               b=0.011, ms=0.01, mi=0.03, mr=0.01,
               dt=500, verbose=TRUE)
mod.r5 <- getDCM(out, run=3)
head(mod.r5)
}
\author{
  Samuel M. Jenness <sjenness@uw.edu>
}
\keyword{model}

