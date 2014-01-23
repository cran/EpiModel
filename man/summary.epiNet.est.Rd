\name{summary.epiNet.est}
\alias{summary.epiNet.est}
\title{Summary for Network Model Fit}
\usage{
\method{summary}{epiNet.est}(object, ...)
}
\arguments{
  \item{object}{an \code{EpiModel} object of class
  \code{epiNet.est}.}

  \item{...}{additional summary function arguments.}
}
\description{
This function prints the summary model fit statistics for
an ERGM or STERGM fit.
}
\details{
This function is simply a wrapper function for
\code{summary.ergm} and \code{summary.stergm}.
}
\keyword{extract}

