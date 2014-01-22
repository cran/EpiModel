\name{modeids}
\alias{modeids}
\title{ID Numbers for Bipartite Network}
\usage{
modeids(nw, mode)
}
\arguments{
  \item{nw}{an object of class \code{network} or
  \code{networkDynamic}.}

  \item{mode}{mode number to return ID numbers for.}
}
\description{
This function outputs ID numbers for a mode number for a
bipartite network.
}
\examples{
nw <- network.initialize(10, bipartite=5)
modeids(nw, mode=2)
}
\seealso{
\code{\link{idmode}} provides the reverse functionality.
}
\keyword{epiNetUtils}
\keyword{internal}

