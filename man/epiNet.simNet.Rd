\name{epiNet.simNet}
\alias{epiNet.simNet}
\title{Simulate Networks from STERGM Model Fit}
\usage{
epiNet.simNet(x, nsteps, nsims = 1, save.stats = TRUE, stats.formula,
  verbose = TRUE, ...)
}
\arguments{
  \item{x}{an \code{EpiModel} object of class
  \code{\link{epiNet.est}}.}

  \item{nsteps}{number of time steps over which to simulate
  relational dynamics.}

  \item{nsims}{number of network simulations.}

  \item{save.stats}{if \code{TRUE}, run simulation-based
  model diagnostics for the formation and dissolution
  processes.}

  \item{stats.formula}{a right-hand sided ERGM formula that
  includes network statistics of interest, with the default
  to the formation formula terms.}

  \item{verbose}{if \code{TRUE}, print simulation progress
  to the console.}

  \item{...}{additional arguments to pass to
  \code{simulate.stergm}.}
}
\description{
This function simulates a set of dynamic networks of class
\code{networkDynamic} given a STERGM model fit with
\code{\link{epiNet.est}}.
}
\details{


This function is a wrapper function around the
\code{simulate} functions in the \code{tergm} package. The
purpose is to simulate a series of dynamic relational
networks for a set number of time steps, and to provide
access to access to diagnostics on those network
simulations.

The intended use of these simulated networks is in
independent network disease models simulated in
\code{\link{epiNet.simTrans}}. With dependent network
models, one bypasses this network simulation step as that
processes is coincident with the disease simulation. See
Section 4 of \href{../doc/Tutorial.pdf}{EpiModel Tutorial}
for further clarification.
}
\examples{
\dontrun{
# Initialize bipartite network
num.m1 <- 100
num.m2 <- 100
nw <- network.initialize(num.m1 + num.m2,
                         bipartite=num.m1, directed=FALSE)

# Set mode-specific fractional degree distributions
deg.dist.m1 <- c(0.40, 0.55, 0.04, 0.01)
deg.dist.m2 <- c(0.48, 0.41, 0.08, 0.03)

# Set model parameters and control arguments
formation <- ~ edges + b1degree(0:1) + b2degree(0:1)
target.stats <- c(66, 40, 55, 48, 41)
dissolution <- ~ offset(edges)
duration <- 25
coef.diss <- dissolution.coefs(dissolution, duration)

# Estimate the model with all the parameters set above
est <- epiNet.est(
  nw,
  formation,
  dissolution,
  target.stats,
  coef.diss)
est

# Simulate 5 networks over 100 timesteps
nwsims <- epiNet.simNet(est, nsteps = 100, nsims = 5)

# Print the object, plot the formation summary, plot the duration summary
nwsims
plot(nwsims, sim=2)
plot(nwsims, type = "duration", sim = 2)
}
}
\seealso{
\code{\link{plot.epiNet.simNet}}, \code{\link{epiNet.est}},
\code{\link{epiNet.simTrans}}
}
\keyword{model}

