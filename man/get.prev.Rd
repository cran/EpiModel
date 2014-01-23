\name{get.prev}
\alias{get.prev}
\title{Get State Size from a NetworkDynamic Object}
\usage{
get.prev(all, at, set.all = TRUE)
}
\arguments{
  \item{all}{a list object containing a
  \code{networkDynamic} object and other initialization
  information passed from \code{\link{epiNet.simTrans}}.}

  \item{at}{time step to query.}

  \item{set.all}{if \code{TRUE}, sets output on all object,
  otherwise will output directly.}
}
\description{
This function provides all active model state sizes from
the network at the specified time step, output to a list of
vectors.
}
\details{
This network utility is used during the
\code{\link{epiNet.simTrans}} simulation process to
efficiently query the current size of each state or
compartment in the model at any given timestep. For a
bipartite network, the current state size for each mode,
and overall is provided. A more flexible, but inefficient
status query function is \code{get.stat}.
}
\seealso{
\code{\link{get.stat}}
}
\keyword{epiNetUtils}
\keyword{internal}

