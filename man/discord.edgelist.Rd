\name{discord.edgelist}
\alias{discord.edgelist}
\title{Discordant Edgelist from NetworkDynamic Object}
\usage{
discord.edgelist(all, infecteds, susceptible, at)
}
\arguments{
  \item{all}{a list object containing a
  \code{networkDynamic} object and other initialization
  information passed from \code{\link{epiNet.simTrans}}.}

  \item{infecteds}{vector of IDs for currently infecteds.}

  \item{susceptible}{vector of IDs for currently
  susceptible.}

  \item{at}{current time step.}
}
\value{
This function returns a \code{data.frame} with the
following columns: \itemize{ \item \strong{time:} time step
queried \item \strong{inf:} ID number for the infected
partner \item \strong{sus:} ID number for the susceptible
partner \item \strong{eid:} Edge ID number for that unique
dyad \item \strong{infdeg:} Current degree for the infected
partner \item \strong{susdeg:} Current degree for the
susceptible partner } The output from this function is
added to the transmission \code{data.frame} object that is
requested as output in \code{epiNet.simTrans} simulations
with the \code{save.trans=TRUE} argument.
}
\description{
This function returns a \code{data.frame} with a discordant
edgelist, defined as the set of edges in which the status
of the two partners is one susceptible and one infected.
}
\details{
This internal function works within the parent
\code{\link{infection}} function to pull the current
edgelist from the dynamic network object, look up the
disease status of the head and tails on the edge, and
subset the list to those edges with one susceptible and one
infected node.
}
\seealso{
\code{\link{epiNet.simTrans}}, \code{\link{infection}}
}
\keyword{epiNetModule}
\keyword{internal}

