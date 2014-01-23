\name{recovery}
\alias{recovery}
\title{Recovery: epiNet.simTrans Module}
\usage{
recovery(all, at)
}
\arguments{
  \item{all}{a list object containing a
  \code{networkDynamic} object and other initialization
  information passed from \code{\link{epiNet.simTrans}}.}

  \item{at}{time step to query.}
}
\description{
This function simulates recovery from the infected state
either to an distinct recovered state (SIR model type) or
back to a susceptible state (SIS model type), for use in
\code{\link{epiNet.simTrans}}.
}
\keyword{internal}

