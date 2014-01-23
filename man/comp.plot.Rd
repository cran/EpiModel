\name{comp.plot}
\alias{comp.plot}
\alias{comp.plot.epiDCM}
\alias{comp.plot.epiICM}
\alias{comp.plot.epiNet.simTrans}
\title{Plot Compartment Diagram from Deterministic or Stochastic Model}
\usage{
comp.plot(x, time, run, digits, ...)

\method{comp.plot}{epiDCM}(x, time, run = 1, digits = 3, ...)

\method{comp.plot}{epiICM}(x, time = 1, run, digits = 3, ...)

\method{comp.plot}{epiNet.simTrans}(x, time = 1, run, digits = 3, ...)
}
\arguments{
  \item{x}{an \code{EpiModel} object of class
  \code{epiDCM}, \code{epiICM}, or \code{epiNet.simTrans}}

  \item{time}{timestep at which summary data is extracted}

  \item{run}{model run number for plotting, only for
  \code{epiDCM} class models with multiple runs
  (sensitivity analyses)}

  \item{digits}{number of significant digits to print}

  \item{...}{additional arguments passed to plot (not
  currently used)}
}
\description{
This function plots a compartment diagram for all three
classes of epidemic models in EpiModel: \code{epiDCM},
\code{epiICM}, and \code{epiNet.simTrans}.
}
\details{
The \code{comp.plot} function is a visual summary of an
EpiModel at a specific point in time. The information
contained in \code{comp.plot} is essentially the same as in
the \code{summary} functions for a model, but presented in
a way commonly seen in the epidemiological literature.

For \code{epiDCM} class plots, one must specify the
specific run number to be plotted, in case the object
contains multiple runs from a sensitivity analysis. For
\code{epiICM} and \code{epiNet.simTrans} class plots, the
\code{run} argument is not used because the plots show the
means and standard deviations for each state size and flow
at the specified timestep.

The functionality of these plots is currently limited to
one-group or one- mode models in each of the three model
classes, and that functionality will be expanded in future
EpiModel releases.
}
\examples{
# Example for epiDCM: SIR model with varying act.rate
mod <- epiDCM(type="SIR", s.num = 1000, i.num = 1,
              trans.rate = 0.2, act.rate = 5:7, rec.rate = 1/3,
              b.rate = 1/90, ds.rate = 1/100, di.rate = 1/35,
              dr.rate = 1/100, nsteps = 25)
comp.plot(mod, time=25, run=3)

# Example for epiICM: SIR model with 3 simulations
modICM <- epiICM(type = "SIR", s.num = 500, i.num = 1,
                 trans.rate = 0.2, act.rate = 3,
                 rec.rate = 1/50, b.rate = 1/100,
                 ds.rate = 1/100, di.rate = 1/90, dr.rate = 1/100,
                 nsteps = 25, nsims = 3)
comp.plot(modICM, time=25)
}
\keyword{plot}

