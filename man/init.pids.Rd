\name{init.pids}
\alias{init.pids}
\title{Persistent ID Initialization}
\usage{
init.pids(nw, prefixes = c("F", "M"))
}
\arguments{
  \item{nw}{an object of class \code{networkDynamic}.}

  \item{prefixes}{character string prefix for mode-specific
  ID, with default to c('F', 'M').}
}
\description{
This function initializes the persistent IDs for a
\code{networkDynamic} object.
}
\details{
This function is used for \code{\link{epiNet.simTrans}}
simulations over bipartite networks for populations with
vital dynamics. Persistent IDs are required in this
situation because when new nodes are added to the first
mode in a bipartite network, the IDs for the second mode
shift upward. Persistent IDs allow for an analysis of
disease transmission chains for these simulations. These
IDs are also invoked in the \code{\link{births}} module
when the persistent IDs of incoming nodes must be set.
}
\examples{
# Initialize network with 25 female and 75 male
nw <- network.initialize(100, bipartite=25)

# Set persistent IDs using the default F/M prefix
nw <- init.pids(nw)
nw \%v\% 'vertex.names'

# Use another prefix combination
nw <- init.pids(nw, c('A', 'B'))
nw \%v\% 'vertex.names'
}
\seealso{
\code{\link{initialize.pids}}
}
\keyword{epiNetModule}
\keyword{internal}

