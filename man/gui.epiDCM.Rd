\name{gui.epiDCM}
\alias{gui.epiDCM}
\title{GUI for epiDCM Models}
\usage{
gui.epiDCM()
}
\description{
This function runs a web browser-based GUI of the epiDCM
application using the \code{shiny} package.
}
\details{
Using the \code{shiny} interface, this function runs a
web-based GUI of a one-group \code{epiDCM} model with user
input on model type, state sizes, and parameters. Model
output may be plotted, summarized, and saved as raw data
using the core \code{EpiModel} functionality for
\code{epiDCM} models.
}
\examples{
\dontrun{
 gui.epiDCM()
}
}
\seealso{
\code{\link{epiDCM}}, \code{\link{plot.epiDCM}},
\code{\link{summary.epiDCM}},
\code{\link{as.data.frame.epiDCM}}
}
\keyword{GUI}

