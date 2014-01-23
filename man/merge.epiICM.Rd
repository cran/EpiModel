\name{merge.epiICM}
\alias{merge.epiICM}
\title{Merge Model Simulations Across epiICM Objects}
\usage{
\method{merge}{epiICM}(x, y, ...)
}
\arguments{
  \item{x}{an \code{EpiModel} object of class
  \code{\link{epiICM}}.}

  \item{y}{another \code{EpiModel} object of class
  \code{\link{epiICM}}, with the identical model
  parameterization as \code{x}.}

  \item{...}{additional arguments required by generic merge
  method, but not used.}
}
\description{
This function merges data from two independent
\code{epiSPM} simulations.
}
\details{
The purpose of this generic merge function is to facilitate
analysis of multiple simulations of \code{\link{epiICM}}
class models that may have been simulated under different
function calls, but where the model parameterization
between the two calls is exactly the same. Such a situation
would occur when one runs sets of simulations in parallel
across computing clusters for efficency.

Note that this merge function does not work the same as the
default merge in allowing for a combined object where the
structure differs between the individual elements: instead,
the function checks that objects are identical in model
parameterization in every respect except number of
simulations.
}
\examples{
x <- epiICM(type="SI", s.num=1000, i.num=100,
            trans.rate = 0.2, act.rate = 0.8,
            nsteps = 10, nsims = 3)

y <- epiICM(type="SI", s.num=1000, i.num=100,
            trans.rate = 0.2, act.rate = 0.8,
            nsteps = 10, nsims = 2)

z <- merge(x, y)
x$i.num
y$i.num
z$i.num
}
\keyword{extract}

