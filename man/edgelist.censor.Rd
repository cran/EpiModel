\name{edgelist.censor}
\alias{edgelist.censor}
\title{Table of Partnership Censoring}
\usage{
edgelist.censor(el)
}
\arguments{
  \item{el}{a timed edgelist with start and end times from
  extracted from a \code{networkDynamic} object using the
  \code{as.data.frame.networkDynamic} function.}
}
\description{
This function outputs a table of the number and percent of
partnerships that are left, right, both, or uncensored for
a \code{networkDynamic} object.
}
\details{
Given a STERGM simulation over a specified number of time
steps, the edges within that simulation may be
left-censored (started before the first step),
right-censored (continued after the last step), right and
left-censored, or uncensored. The amount of censoring will
tend to increase when the average partnership duration
approaches the length of the simulation.
}
\examples{
nw <- network.initialize(n = 100, directed = FALSE)
nw <- set.vertex.attribute(nw, "race", value = rep(0:1, each = 50))
formation <- ~ edges + nodematch("race")
target.stats <- c(45, 25)
dissolution <- ~ offset(edges)
coef.diss <- dissolution.coefs(dissolution, duration = 20)

# Model estimation, then simulation
est <- epiNet.est(nw,
                  formation,
                  dissolution,
                  target.stats,
                  coef.diss,
                  save.stats = FALSE,
                  verbose = FALSE)
nwsims <- epiNet.simNet(est, nsteps = 100,
                        nsims = 1, verbose = FALSE)

# Extract edgelist from first simulation, then calculate censoring
el.sim1 <- as.data.frame(nwsims$nwsim1)
edgelist.censor(el.sim1)
}
\keyword{epiNetUtils}

