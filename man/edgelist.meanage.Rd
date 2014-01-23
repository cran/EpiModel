\name{edgelist.meanage}
\alias{edgelist.meanage}
\title{Mean Age of Partnerships over Time}
\usage{
edgelist.meanage(x, nw, el, dissolution)
}
\arguments{
  \item{x}{an \code{EpiModel} object of class
  \code{\link{epiNet.est}}.}

  \item{nw}{if not passing \code{x}, a
  \code{networkDynamic} object.}

  \item{el}{if not passing \code{x}, a timed edgelist from
  a \code{networkDynamic} object extracted with the
  \code{as.data.frame.networkDynamic} function.}

  \item{dissolution}{if not passing \code{x}, a right-hand
  sided STERGM dissolution formula (see
  \code{\link{epiNet.est}}).}
}
\description{
This function outputs a vector of mean ages of partnerships
at a series of timesteps
}
\details{
This function calculates the mean partnership age at each
time step over a dynamic network simulation from
\code{\link{epiNet.est}}. These objects contain the
network, edgelist, and dissolution objects needed for the
calculation. Alternatively, one may pass in these objects
separately if \code{epiNet.est} was not used, or statistics
were not run requested after the estimation.

Currently, the calculations are limited to those
dissolution formulas with a single homogenous dissolution
(\code{~offset(edges)}). This functionality will be
expanded in future releases.
}
\examples{
nw <- network.initialize(n = 100, directed = FALSE)
nw <- set.vertex.attribute(nw, "race", value = rep(0:1, each = 50))
formation <- ~ edges + nodematch("race")
target.stats <- c(45, 25)
dissolution <- ~ offset(edges)
coef.diss <- dissolution.coefs(dissolution, duration = 20)

# Model estimation
est <- epiNet.est(
  nw,
  formation,
  dissolution,
  target.stats,
  coef.diss,
  save.stats = FALSE,
  verbose = FALSE)

# Get the formation coefficients, then simulate a dynamic network
# from a simulated static network
coef.form <- est$coef.form
sim <- simulate(
  simulate(est$fit),
  formation = formation,
  dissolution = dissolution,
  coef.form = coef.form,
  coef.diss = coef.diss[[3]],
  time.slices = 200)

# Gather the three objects needed for the age calculation
nw <- est$fit$network
el <- as.data.frame(sim)
dissolution <- est$dissolution
edgelist.meanage(nw=nw, el=el, dissolution=dissolution)

# Alternatively, epiNet.est automatically calculates these with stats = TRUE
est <- epiNet.est(
  nw,
  formation,
  dissolution,
  target.stats,
  coef.diss,
  edapprox = TRUE,
  save.stats = TRUE,
  stats.end = 200)
est$pages
}
\keyword{epiNetUtils}
\keyword{internal}

