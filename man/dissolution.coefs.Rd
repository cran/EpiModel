\name{dissolution.coefs}
\alias{dissolution.coefs}
\title{Calcuate Dissolution Coefficients for Stochastic Network Models}
\usage{
dissolution.coefs(dissolution, duration, d.rate = 0)
}
\arguments{
  \item{dissolution}{a right-hand sided STERGM dissolution
  formula (see \code{\link{epiNet.est}}); currently limited
  only to simple \code{~offset(edges)} dissolution models.}

  \item{duration}{an edge duration in arbitrary time
  units.}

  \item{d.rate}{background death rate in the absence of
  disease.}
}
\value{
A list of with the following elements:
}
\description{
This function calculates dissolution coefficients from a
duration to pass as terms to an ERGM/STERGM fit and applies
a coefficient correction to account for deaths.
}
\details{
This function performs two duties prior to the estimation
phase of \code{epiNet} class stochastic network models.
\enumerate{ \item \strong{Transforming duration to
coefficient form:} an average duration of edges in a
population must be mathematically transformed to logit
coefficients. The theory and mathematics of the
transformation are further developed in the
\code{networkDynamic} package vignette.  \item
\strong{Adjusting coefficients to account for death:} In a
dynamic network simulation in an open population (in which
there are births and deaths), it is necessary to adjust the
dissolution coefficient for the STERGM simulations to
account for the death as a competing risk to edge
dissolution. }

Future releases of this software will allow for more
flexibility in the possible dissolution models that may be
calculated here, including models that have heterogenous
dissolution probabilities conditional on nodal or edge
attributes.
}
\examples{
dissolution <- ~offset(edges)
duration <- 25
dissolution.coefs(dissolution, duration)
dissolution.coefs(dissolution, duration, d.rate=0.001)
}
\seealso{
For a detailed explanation of this function, please see the
HTML vignette: \href{../doc/epiNetUtils.html}{EpiModel
Network Utility Functions}.
}
\keyword{epiNetUtils}

