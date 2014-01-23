\name{init.status}
\alias{init.status}
\title{Disease Status Initialization Module for epiNet.simTrans}
\usage{
init.status(all, tea = FALSE)
}
\arguments{
  \item{all}{a list object containing a
  \code{networkDynamic} object and other initialization
  information passed from \code{\link{epiNet.simTrans}}.}

  \item{tea}{if \code{TRUE}, set status in a temporally
  extended attribute, otherwise use a static vertex
  attribute.}
}
\description{
This function sets the initial disease status on the
network given the supplied prevalence.
}
\details{
This internal function sets, either randomly or
deterministically, the nodes that are infected at the
starting time of \code{epiNet} simulations, \eqn{t_1}. If
the number to be initially infected is passed, this
function may set the initial number infected based on the
number specified, either as a a set of random draws from a
binomial distribution or as the exact number specified. In
either case, the specific nodes infected are a random
sample from the network. In contrast, a set of specific
nodes may be infected by passing the vector to
\code{\link{epiNet.simTrans}}.
}
\seealso{
This is an initialization module for
\code{\link{epiNet.simTrans}}.  It precedes the infection
time module in \code{\link{init.inf.time}}.
}
\keyword{epiNetModule}
\keyword{internal}

