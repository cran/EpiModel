\name{get.stat}
\alias{get.stat}
\title{Get Active Status from NetworkDynamic Object}
\usage{
get.stat(all, stat = 1, at, out, mode, require.active = TRUE)
}
\arguments{
  \item{all}{a list object containing a
  \code{networkDynamic} object and other initialization
  information passed from \code{\link{epiNet.simTrans}}.}

  \item{stat}{a vector of status values in numeric to query
  (the standard status representation is 0 = negative, 1 =
  infected, 2 = recovered).}

  \item{at}{time step to query.}

  \item{out}{function output, with options of
  \code{out="vec"} for a vector of status, \code{out="ids"}
  for a vector of IDs with status equal to \code{stat}, and
  \code{out="prev"} for the number current with status
  equal to \code{stat}, and \code{out="all"} to return a
  list of the prior three elements.}

  \item{mode}{if \code{nw} is bipartite, the mode number
  for status (may be ignored if requesting output for both
  modes).}

  \item{require.active}{if \code{TRUE}, NA will be returned
  instead of an attribute value if the associated vertex or
  edge is inactive for the query period (see
  \code{\link{get.vertex.attribute.active}}).}
}
\description{
This function outputs information on the disease or outcome
status stored in the networkDyamic status TEA used by
\code{\link{epiNet.simTrans}}.
}
\details{
This network utility is used during the
\code{\link{epiNet.simTrans}} simulation process to query
the ID numbers and size of currently active nodes with a
specific disease status. A more efficient query just for
state sizes of status = 1 is found in \code{get.prev}.
}
\seealso{
\code{\link{get.prev}}
}
\keyword{epiNetUtils}
\keyword{internal}

