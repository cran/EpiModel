\name{plotComp}
\alias{plotComp}
\title{Plot Compartment Diagram}
\usage{
  plotComp(type)
}
\arguments{
  \item{type}{Type of deterministic model to plot.}
}
\description{
  This function plots the compartment diagram for standard
  models included in EpiModel
}
\details{
  Current options only are: \itemize{ \item \code{'SIR'}
  simple SIR models \item \code{'SIRvd'} SIR models with
  vital dynamics }
}
\examples{
plotComp(type='SIR')
plotComp(type='SIRvd')
}
\author{
  Samuel M. Jenness <sjenness@uw.edu>
}
\seealso{
  \code{\link{epiDCM}}
}
\keyword{model}

