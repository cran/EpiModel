\name{deaths.rec}
\alias{deaths.rec}
\title{Deaths for Recovered: epiNet.simTrans Module}
\usage{
deaths.rec(all, at)
}
\arguments{
  \item{all}{a list object containing a
  \code{networkDynamic} object and other initialization
  information passed from \code{\link{epiNet.simTrans}}.}

  \item{at}{time step for simulation.}
}
\description{
This function simulates death for those in the susceptible
state for use in \code{\link{epiNet.simTrans}} simulations.
}
\seealso{
\code{\link{epiNet.simTrans}}, \code{\link{deaths.sus}},
\code{\link{deaths.inf}}
}
\keyword{epiNetModule}
\keyword{internal}

