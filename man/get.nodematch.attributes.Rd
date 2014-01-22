\name{get.nodematch.attributes}
\alias{get.nodematch.attributes}
\title{Get Nodematch Term Proportions}
\usage{
get.nodematch.attributes(nw, formation)
}
\arguments{
  \item{nw}{an object of class \code{network} or
  \code{networkDynamic}.}

  \item{formation}{a right-hand sided STERGM formation
  formula (see \code{\link{epiNet.est}}).}
}
\description{
This function outputs a list of values and proportions of
static vertex attributes present in the formation formula
as nodematch terms.
}
\details{
Incoming nodes in a dynamic network simulation require
values for any attributes that are used in the formation
formula through nodematch network terms. This function gets
the proportional distribution of those attribute values of
existing nodes in the network. The distributions may then
be fed into \code{set.nodematch.attributes} to update the
attributes of new nodes.
}
\examples{
# Set network with two vertex attributes
nw <- network.initialize(100)
nw <- set.vertex.attribute(nw, "race", value = rep(0:1, each=50))
nw <- set.vertex.attribute(nw, "sex",
                           value = c(rep("M", 10), rep("F", 90)))

# Formation formula uses both attributes
formation <- ~ edges + nodematch("race") + nodematch("sex")

# Get existing proportions of those two terms
get.nodematch.attributes(nw, formation)
}
\seealso{
\code{\link{set.nodematch.attributes}}
}
\keyword{epiNetUtils}
\keyword{internal}

