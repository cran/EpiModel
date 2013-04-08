\name{brewer.ramp}
\alias{brewer.ramp}
\title{RColorBrewer Color Ramp for EpiModel Plots}
\usage{
  brewer.ramp(n, plt, del.lights = TRUE)
}
\arguments{
  \item{plt}{RColorBrewer palette from
  \code{\link{brewer.pal}}}

  \item{n}{number of colors to return}

  \item{del.lights}{delete the lightest colors from the
  color palette}
}
\description{
  This function returns vector of colors consistent with a
  high-brightness set of colors from an RColorBrewer
  palette.
}
\details{
  RColorBrewer provides easy access to helpful color
  palettes, but the built-in palettes are limited to the
  set of colors in the existing palette. This function
  expands the palette size to any number of colors by
  filling in the gaps. Also, colors within the 'div' and
  'seq' set of palettes whose colors are very light (close
  to white) are deleted by default for better visualization
  of plots.
}
\examples{
plot(1:100, 1:100, type='n')
abline(v=1:100, lwd=6, col=brewer.ramp(100, 'Spectral'))
abline(v=1:100, lwd=6, col=brewer.ramp(100, 'Greys'))
abline(v=1:100, lwd=6, col=brewer.ramp(100, 'Blues'))
abline(v=1:100, lwd=6, col=brewer.ramp(100, 'Set1'))
}
\author{
  Samuel M. Jenness <sjenness@uw.edu>
}
\seealso{
  \code{\link{RColorBrewer}}
}
\keyword{color}

