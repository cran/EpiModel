\name{summary.epiICM}
\alias{summary.epiICM}
\title{Summary Model Statistics}
\usage{
\method{summary}{epiICM}(object, time, digits = 3, comp.plot = FALSE, ...)
}
\arguments{
  \item{object}{an \code{EpiModel} object of class
  \code{epiICM}.}

  \item{time}{time step for model statistics.}

  \item{digits}{number of significant digits to print.}

  \item{comp.plot}{if \code{TRUE}, plot compartments and
  flows in summary. This can be called separately using the
  \code{\link{comp.plot}} function.}

  \item{...}{additional summary function arguments.}
}
\description{
This function extracts and prints model statistics solved
with the \code{epiICM} function.
}
\details{
Summary statistics for the main epidemiological outcomes
(state and transition size and prevalence) from an
\code{epiICM} model. Time-specific summary measures are
provided, so it is necessary to input a time of interest.
One may simultaneously obtain console-based statistics and
a compartment plot via the \code{comp.plot} function
through the \code{comp.plot} argument.
}
\examples{
\dontrun{
# Stochastic SI model with 10 simulations
mod <- epiICM(type = "SI", s.num = 500, i.num = 1,
              trans.rate = 0.2, act.rate = 0.25,
              nsteps = 500, nsims = 10)
summary(mod, time = 25)
summary(mod, time = 50, comp.plot = TRUE)
}
}
\seealso{
\code{\link{epiICM}}, \code{\link{comp.plot}}
}
\keyword{extract}

