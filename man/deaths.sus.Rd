\name{deaths.sus}
\alias{deaths.sus}
\title{Deaths for Suceptibles: epiNet.simTrans Module}
\usage{
deaths.sus(all, at)
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
\code{\link{epiNet.simTrans}}, \code{\link{deaths.inf}},
\code{\link{deaths.rec}}
}
\keyword{epiNetModule}
\keyword{internal}

