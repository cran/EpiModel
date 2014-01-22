\name{bip.degdist.check}
\alias{bip.degdist.check}
\title{Check Degree Distribution for Bipartite Target Statistics}
\usage{
bip.degdist.check(num.m1, num.m2, deg.dist.m1, deg.dist.m2)
}
\arguments{
  \item{num.m1}{number of vertices in mode 1.}

  \item{num.m2}{number of vertices in mode 2.}

  \item{deg.dist.m1}{vector with fractional degree
  distribution for mode 1.}

  \item{deg.dist.m2}{vector with fractional degree
  distribution for mode 2.}
}
\description{
This function checks for consistency in the degree
distributions between the two modes in a bipartite network,
and if balance is present returns the implied network
statistics for use in a ERGM.
}
\details{
This function outputs the number of nodes of degree 0 to m,
where m is the length of a fractional degree distribution
vector, given that vector and the size of the mode. This
utility is used to check for balance in implied degree
given that fractional distribution within bipartite network
simulations, in which the degree-contrained counts must
match across modes.
}
\examples{
# An imbalanced distribution
bip.degdist.check(num.m1 = 500, num.m2 = 500,
                  deg.dist.m1 = c(0.40, 0.55, 0.03, 0.02),
                  deg.dist.m2 = c(0.48, 0.41, 0.08, 0.03))

# A balanced distribution
targets <- bip.degdist.check(num.m1 = 500, num.m2 = 500,
                  deg.dist.m1 = c(0.40, 0.55, 0.04, 0.01),
                  deg.dist.m2 = c(0.48, 0.41, 0.08, 0.03))
targets
}
\seealso{
For a detailed explanation of this function, see the
included HTML vignette:
\href{../doc/epiNetUtils.html}{EpiModel Network Utility
Functions.}
}
\keyword{epiNetUtils}

