\name{node.active}
\alias{node.active}
\title{Query Active Nodes in NetworkDynamic Object}
\usage{
node.active(nw, at, out, mode, active.default = FALSE)
}
\arguments{
  \item{nw}{an object of class \code{networkDynamic}.}

  \item{at}{time step to query.}

  \item{out}{function output, with options of
  \code{out="vec"} for a T/F vector of whether the node is
  active, \code{out="ids"} for a vector of IDs active,
  \code{out="prev"} for the number of nodes that are
  active, and \code{out="all"} to return a list of the
  prior three elements.}

  \item{mode}{if \code{nw} is bipartite, the mode number
  for status (may be ignored if requesting output for both
  modes).}

  \item{active.default}{if \code{TRUE}, elements without an
  activity attribute will be regarded as active.}
}
\description{
This function outputs information on the active nodes in a
\code{networkDynamic} object.
}
\details{
This is a specialized version of \code{\link{is.active}}
from the \code{networkDynamic} package that allows for key
output to be efficiently generated for use in
\code{\link{epiNet.simTrans}} simulations. For functions
that query active nodes with specific values of the disease
status TEA, use \code{\link{get.stat}} or
\code{\link{get.prev}}.
}
\examples{
# Initialize NW and activate vertices
nw <- network.initialize(20)

# Activate all vertices, then deactive half at time 5
activate.vertices(nw, onset = 1, terminus = 10)
deactivate.vertices(nw, onset = 5, terminus = 10, v = 1:10)

# Output all information for vertices at time 1 and time 5
node.active(nw, at = 1, out = "all")
node.active(nw, at = 5, out = "all")
}
\seealso{
\code{\link{is.active}}, \code{\link{get.stat}},
\code{\link{get.prev}}
}
\keyword{epiNetUtils}
\keyword{internal}

