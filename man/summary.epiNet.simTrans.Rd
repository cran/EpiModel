\name{summary.epiNet.simTrans}
\alias{summary.epiNet.simTrans}
\title{Summary Model Statistics}
\usage{
\method{summary}{epiNet.simTrans}(object, time, digits = 3,
  comp.plot = FALSE, ...)
}
\arguments{
  \item{object}{an \code{EpiModel} object of class
  \code{epiNet.simTrans}.}

  \item{time}{timestep for model statistics.}

  \item{digits}{number of significant digits to print.}

  \item{comp.plot}{if \code{TRUE}, plot compartments and
  flows in summary. This can be called separately using the
  \code{\link{comp.plot}} function.}

  \item{...}{additional summary function arguments.}
}
\description{
This function extracts and prints model statistics solved
with the \code{epiNet.simTrans} function.
}
\details{
Summary statistics for the main epidemiological outcomes
(state and transition size and prevalence) from an
\code{epiNet.simTrans} model. Time-specific summary
measures are provided, so it is necessary to input a time
of interest. One may simultaneously obtain console-based
statistics and a compartment plot via the \code{comp.plot}
function through the \code{comp.plot} argument.
}
\examples{
## See EpiModel Tutorial vignette ##
}
\seealso{
\code{\link{epiNet.simTrans}}, \code{\link{comp.plot}}
}
\keyword{extract}

