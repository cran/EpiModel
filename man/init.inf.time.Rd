\name{init.inf.time}
\alias{init.inf.time}
\title{Infection Time Module for epiNet.simTrans}
\usage{
init.inf.time(all)
}
\arguments{
  \item{all}{a list object containing a
  \code{networkDynamic} object and other initialization
  information passed from \code{\link{epiNet.simTrans}}.}
}
\description{
This function sets the initial time of infection for those
with disease status infected.
}
\details{
This module sets the time of infection for those nodes set
infected at the starting time of \code{epiNet} simulations,
\eqn{t_1}. For vital dynamics models, the infection time
for those nodes is a random draw from an exponential
distribution with the rate parameter defined by the
\code{di.rate} argument. For models without vital dynamics,
the infection time is a random draw from a uniform
distribution of integers with a minimum of 1 and a maximum
of the number of time steps in the model. In both cases, to
set the infection times to be in the past, these times are
multiplied by -1, and 2 is added to allow for possible
infection times up until step 2, when the disease
simulation time loop starts.
}
\seealso{
This is an initialization module for
\code{\link{epiNet.simTrans}}.  It follows the infection
status module in \code{\link{init.status}}.
}
\keyword{epiNetModule}
\keyword{internal}

