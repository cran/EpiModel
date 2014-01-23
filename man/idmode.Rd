\name{idmode}
\alias{idmode}
\title{Mode Numbers for Bipartite Network}
\usage{
idmode(nw, ids)
}
\arguments{
  \item{nw}{an object of class \code{network} or
  \code{networkDynamic}.}

  \item{ids}{a vector of ID numbers for which the mode
  number should be returned.}
}
\description{
This function outputs mode numbers give ID numbers for a
bipartite network.
}
\examples{
nw <- network.initialize(10, bipartite = 5)
idmode(nw)
idmode(nw, ids = 3)
}
\seealso{
\code{\link{modeids}} provides the reverse functionality.
}
\keyword{epiNetUtils}
\keyword{internal}

