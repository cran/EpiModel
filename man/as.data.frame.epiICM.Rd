\name{as.data.frame.epiICM}
\alias{as.data.frame.epiICM}
\alias{as.data.frame.epiNet.simTrans}
\title{Extract Model Run for Stochastic Models}
\usage{
\method{as.data.frame}{epiICM}(x, row.names = NULL, optional = FALSE, sim,
  out = "mean", ...)

\method{as.data.frame}{epiNet.simTrans}(x, row.names = NULL,
  optional = FALSE, sim, out = "mean", ...)
}
\arguments{
  \item{x}{an \code{EpiModel} object of class \code{epiICM}
  or \code{epiNet.simTrans}.}

  \item{sim}{simulation number from model, used only if
  \code{out="vals"}.}

  \item{out}{type of output to data frame: use
  \code{"mean"} for mean across simulations, \code{"sd"}
  for standard deviations across all simulations, and
  \code{"vals"} for values from one simulation (with
  simulation set with \code{sim}).}

  \item{row.names}{see
  \code{\link{as.data.frame.default}}.}

  \item{optional}{see \code{\link{as.data.frame.default}}.}

  \item{...}{see \code{\link{as.data.frame.default}}.}
}
\description{
This function extracts model simulations for objects of
classes \code{epiICM} and \code{epiNet.simTrans} into a
data frame using the generic \code{as.data.frame} function.
}
\details{
These generics work from both \code{epiICM} and
\code{epiNet.simTrans} class objects, and output data
frames containing the time-specific means, standard
deviations, or model values (compartment and flow sizes
from individual simulations) from these stochastic model
classes. For summary output types, row mean or standard
deviations are calculated across all simulations. Total
population size for one and two-group models are also
calculated and appended.
}
\examples{
# Individual contact SIS model with 2 simulations
mod <- epiICM(type="SIS", s.num=500, i.num=1,
              trans.rate=0.8, act.rate=2, rec.rate=0.1,
              nsteps=20, nsims=2)

# Default output is mean across simulations
as.data.frame(mod)

# Standard deviations of simulations
as.data.frame(mod, out="sd")

# Individual simulation runs, with default sim=1
as.data.frame(mod, out="vals")
as.data.frame(mod, out="vals", sim=2)
}
\keyword{extract}

