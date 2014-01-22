\name{epiNet.est}
\alias{epiNet.est}
\title{Stochastic Network Estimation with STERGM}
\usage{
epiNet.est(nw, formation, dissolution, target.stats, coef.diss, constraints,
  edapprox = TRUE, save.stats = TRUE, stats.formula, stats.start = 1,
  stats.end = 1000, verbose = TRUE, ...)
}
\arguments{
  \item{nw}{an object of class \code{network}.}

  \item{formation}{a right-hand sided STERGM formation
  formula in the form \code{~ edges + ...}, where
  \code{...} are additional network statistics.}

  \item{dissolution}{a right-hand sided STERGM dissolution
  formula in the form \code{~ offset(edges) + ...}, where
  \code{...} are additional network statistics.}

  \item{target.stats}{a vector of target statistics for
  formation model, with one number for each network
  statistic in the model (see \code{\link{stergm}}).}

  \item{coef.diss}{an object of class
  \code{dissolution.coefs} that is calculated with the
  function of that name containing the target coefficients
  to target.}

  \item{constraints}{a right-hand sided formula specifying
  contraints for the modeled network, in the form
  \code{~...}, where \code{...} are constraint terms
  described in \code{\link{stergm}}. By default, no
  constraints are set.}

  \item{edapprox}{if \code{TRUE}, use the indirect edges
  dissolution approximation method for dynamic modeling
  fitting rather than the more time-intensive full STERGM
  estimation (available with \code{FALSE}). For this to
  function, the \code{coef.diss} terms should be in the
  same order as the formation formula terms.}

  \item{save.stats}{if \code{TRUE}, run simulation-based
  model diagnostics for the formation and dissolution
  processes.}

  \item{stats.formula}{a right-hand sided ERGM formula that
  includes network statistics of interest, with the default
  to the formation formula terms.}

  \item{stats.start}{starting timestep for calculating the
  summary measure diagnostics.}

  \item{stats.end}{ending timestep for calculating the
  summary measure diagnostics.}

  \item{verbose}{if \code{TRUE}, print simulation progress
  to the console.}

  \item{...}{additional arguments to pass to \code{ergm} or
  \code{stergm}.}
}
\description{
Estimates and diagnoses a statistical network model using
the exponential random graph modeling (ERGM) framework with
extensions for dynamic/temporal models (STERGM)
}
\details{


This function is a wrapper function around \code{ergm} and
\code{stergm} to facilitate estimating a dynamic/temporal
ERGM from target statistics. This framework is suggested
for parameterizing a stochastic network epidemic model from
empirical partnership formation and dissolution data
collected in an egocentric format (i.e., study subjects
asked about their recent partners, including questions on
partnership duration).

Additional functionality here includes time-varying
diagnostics for a temporal simulation from the fitted
model, to ensure that target statistics are approximated by
the network simulation. This is available through the
\code{save.stats} argument, with the \code{stats.formula}
argument allowing flexibility in network statistics of
interest outside those in the formation formula. The
diagnostics are available by printing and plotting the
fitted \code{epiNet.est} object.

With a fitted and diagnosed network model, one may proceed
to \code{\link{epiNet.simNet}} for simulating a series of
networks for use in an independent stochastic network
epidemic model, or straight to the epidemic model function
\code{\link{epiNet.simTrans}} if there are dependencies
between the network model and the disease model. See
Section 4 of \href{../doc/Tutorial.pdf}{EpiModel Tutorial}
for further clarification.
}
\examples{
\dontrun{
# Initialize an empty network of 500 nodes divided equally into two races
nw <- network.initialize(n=500, directed=F)
nw <- set.vertex.attribute(nw, "race", value = rep(0:1, each=250))

# Set formation and dissolution formulas
formation <- ~ edges + nodematch("race") + degree(0) + concurrent
dissolution <- ~ offset(edges)

# Set target statistics for formation
target.stats <- c(225, 187, 180, 90)

# Set a vector of durations (here in months, but timestep units are arbitrary)
#  Then use the dissolution.coefs functions to obtain the offset coefficients
durations <- 20
coef.diss <- dissolution.coefs(dissolution, durations)
coef.diss

# Set the stats.formula to include more degree tergms than the formation formula
dx.stats <- ~edges + nodematch("race") + degree(0:5)

# Estimate the STERGM with all the parameters set above with diagnostic
#  simulations through 2000 timesteps, using the edges dissolution approximation
est <- epiNet.est(
  nw,
  formation,
  dissolution,
  target.stats,
  coef.diss,
  edapprox = TRUE,
  stats.formula = dx.stats,
  dx.end = 2000)

# Estimate the STERGM directly and use the default for model statistics
est2 <- epiNet.est(
  nw,
  formation,
  dissolution,
  target.stats,
  coef.diss,
  edapprox = FALSE)
est2

# Print and plot summaries
est
plot(est)
plot(est, type="duration")
}
}
\references{
Krivitsky PN, Handcock MS (2010). A Separable Model for
Dynamic Networks.  \url{http://arxiv.org/abs/1011.1937}.
}
\seealso{
\code{\link{stergm}}, \code{\link{epiNet.simNet}},
\code{\link{epiNet.simTrans}}
}
\keyword{model}

