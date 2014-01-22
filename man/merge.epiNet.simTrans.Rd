\name{merge.epiNet.simTrans}
\alias{merge.epiNet.simTrans}
\title{Merge Model Simulations Across epiNet.simTrans Objects}
\usage{
\method{merge}{epiNet.simTrans}(x, y, keep.trans = TRUE,
  keep.network = TRUE, keep.stats = TRUE, ...)
}
\arguments{
  \item{x}{an \code{EpiModel} object of class
  \code{\link{epiNet.simTrans}}.}

  \item{y}{another \code{EpiModel} object of class
  \code{\link{epiNet.simTrans}}, with the identical model
  parameterization as \code{x}.}

  \item{keep.trans}{if \code{TRUE}, keep in the merged
  output the transmission data frames in \code{trans} and
  the status matrix in \code{stat.mat} (not only not using
  temporally extended attributes) from the original
  \code{x} and \code{y} elements.}

  \item{keep.network}{if \code{TRUE}, keep in the merged
  output the \code{networkDynamic} objects from the
  original \code{x} and \code{y} elements.}

  \item{keep.stats}{if \code{TRUE}, keep in the merged
  output the network statistics (as set by the
  \code{stats.formula} parameter in \code{epiNet.simTrans})
  from the original \code{x} and \code{y} elements.}

  \item{...}{additional arguments required by generic merge
  method, but not used.}
}
\description{
This function merges data from two independent
\code{epiNet.simTrans} simulations.
}
\details{
The purpose of this generic merge function is to facilitate
analysis of multiple simulations of
\code{\link{epiNet.simTrans}} class models that may have
been simulated under different function calls, but where
the model parameterization between the two calls is exactly
the same. Such a situation would occur when one runs sets
of simulations in parallel across computing clusters for
efficency.

Note that this merge function does not work the same as the
default merge function in generating a combined object
where the structure differs between the individual
elements: instead, the function checks that objects are
identical in model parameterization in every respect except
number of simulations.
}
\examples{
nw <- network.initialize(n = 100, directed = FALSE)
dissolution <- ~offset(edges)
coef.diss <- dissolution.coefs(dissolution, duration=10)
est <- epiNet.est(nw,
                  formation = ~ edges,
                  dissolution = ~offset(edges),
                  target.stats = 25,
                  coef.diss = coef.diss,
                  save.stats = FALSE, verbose = FALSE)
nwsims <- epiNet.simNet(est, nsteps = 10, verbose = FALSE)

x <- epiNet.simTrans(nwsims, type = "SI",
                     i.num = 1, trans.rate = 0.9,
                     sims.per.nw = 2, verbose = FALSE)
y <- epiNet.simTrans(nwsims, type = "SI",
                     i.num = 1, trans.rate = 0.9,
                     sims.per.nw = 3, verbose = FALSE)

z <- merge(x, y)
x$i.num
y$i.num
z$i.num
}
\keyword{extract}

