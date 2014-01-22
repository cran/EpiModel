\name{colorTEA}
\alias{colorTEA}
\title{Creates a TEA Variable for Infection Status for \code{ndtv} Animations}
\usage{
colorTEA(nw, old.var = "status", old.sus = 0, old.inf = 1, old.rec = 2,
  new.var = "ndtvcol", new.sus, new.inf, new.rec, verbose = TRUE)
}
\arguments{
  \item{nw}{an object of class \code{networkDynamic}.}

  \item{old.var}{old TEA variable name.}

  \item{old.sus}{status value for susceptible in old TEA
  variable.}

  \item{old.inf}{status value for infected in old TEA
  variable.}

  \item{old.rec}{status value for recovered in old TEA
  variable.}

  \item{new.var}{new TEA variable name to be stored in nw
  object.}

  \item{new.sus}{status value for susceptible in new TEA
  variable.}

  \item{new.inf}{status value for infected in new TEA
  variable.}

  \item{new.rec}{status value for recovered in new TEA
  variable.}

  \item{verbose}{print progress for calculations.}
}
\description{
This function creates a new temporally-extended attribute
(TEA) variable in a \code{networkDynamic} object with color
name output from an existing status variable with numeric
format.
}
\details{
The \code{ndtv} package allows for animated plots of
dynamic network objects, showing the evolving partnership
structure. The \code{EpiModel} package uses
temporally-extended attributes (TEAs) to store longitudinal
disease status for every vertex in a dynamic network. To
visualize disease status dynamically in \code{ndtv}, it is
currently necessary to create a TEA containing the colors
to be used in drawing the nodes.

The convention in \code{\link{plot.epiNet.simTrans}} is to
color the susceptible nodes as blue, infected nodes as red,
and recovered nodes as green. This function allows the user
to take the existing status TEA on the
\code{networkDynamic} object and transform that numeric
variable into any colors.
}
\examples{
## See EpiModel Tutorial vignette ##
}
\seealso{
\code{\link{epiNet.simTrans}} and the \code{ndtv} package
documentation.
}
\keyword{colorUtils}

