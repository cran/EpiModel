\name{set.nodematch.attributes}
\alias{set.nodematch.attributes}
\title{Set Nodematch Values from Proportions}
\usage{
set.nodematch.attributes(nw, nm.attr, new.nodes)
}
\arguments{
  \item{nw}{an object of class \code{network} or
  \code{networkDynamic}.}

  \item{nm.attr}{object output from
  \code{\link{get.nodematch.attributes}}.}

  \item{new.nodes}{vector of new nodes with missing (that
  is, \code{NA}) attribute values.}
}
\description{
This function sets the values of static vertex attributes
present in the formation formula as nodematch terms.
}
\details{
Incoming nodes in a dynamic network simulation require
values for any attributes that are used in the formation
formula through nodematch network terms. This function sets
the values of those attributes based on the proportional
distribution of the existing node attributes. The
distribution must first be obtained from
\code{get.nodematch.attributes}.
}
\examples{
# Set network with two vertex attributes and check distribution
nw <- network.initialize(100)
nw \%v\% "race" <- rep(0:1, each=50)
nw \%v\% "sex" <- c(rep("M", 10), rep("F", 90))
prop.table(table(nw \%v\% "race"))
prop.table(table(nw \%v\% "sex"))

# Formation formula uses both attributes
formation <- ~ edges + nodematch("race") + nodematch("sex")

# Need to get existing proportions of those two terms first
nm.attr <- get.nodematch.attributes(nw, formation)

# Now add 50 more nodes, and set the new attributes in proportion
add.vertices(nw, nv=50)
nw <- set.nodematch.attributes(nw, nm.attr, new.nodes=101:150)

# Some random variation due to sampling new node values
prop.table(table(nw \%v\% "race"))
prop.table(table(nw \%v\% "sex"))
}
\seealso{
\code{\link{get.nodematch.attributes}}
}
\keyword{epiNetUtils}
\keyword{internal}

