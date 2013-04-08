\name{epiABM}
\alias{epiABM}
\title{Stochastic Agent-Based SI-SIS-SIR Models}
\usage{
  epiABM(type, s.num, i.num, r.num = 0, beta, cont, nu = 0,
    nsteps, nsims, verbose = TRUE)
}
\arguments{
  \item{type}{model type (choice of 'SI', 'SIR', or 'SIS')}

  \item{s.num}{number of initial susceptible in population}

  \item{i.num}{number of initial infected in population}

  \item{r.num}{number of initial recovered in poulation}

  \item{beta}{probability of infection per contact}

  \item{cont}{average number of contacts per person per
  unit time}

  \item{nu}{average rate of recovery (1/duration of
  disease)}

  \item{nsteps}{number of time steps to simulate over}

  \item{nsims}{number of simulations to run}

  \item{verbose}{print progress}
}
\description{
  This function simulates an agent-based model of an
  infectious disease with random mixing in the population.
}
\details{
  This function provides flexibility in model choice of a
  susceptible-infected (SI), susceptible-infected-recovered
  (SIR), or a susceptible-infected-susceptible (SIS)
  epidemic model. For an SI model, use the defaults of 0
  for \code{r.num}, \code{nu}, and \code{mu}. For an SIR
  model, specify the recovery rate with the \code{nu}
  parameter, but leave \code{mu} to its 0 default. For an
  SIS model, specify the re-susceptiblity with the
  \code{mu} parameter, but leave the \code{nu} to its 0
  default.

  The model has two sources of stochasticity. First, the
  contact structure is random, with new partnerships
  selected at random from among the actors. The contact
  rate argument \code{cont} is a rate per person per
  timestep, and is transformed to the number of new
  partnerships per timestep by cont*(s.num+i.num)/2.
  Second, the probability of transmission within each
  partnership is based on a standard uniform random deviate
  falling below \code{beta}.
}
\examples{
# SI model
out <- epiABM(type='SI', s.num=500, i.num=1, beta=0.2, cont=0.25, nsteps=250, nsims=10)
par(mar=c(3,3,1,1), mgp=c(2,1,0))
plotABM(out, 'i.num', medline=FALSE)

# SIR model
out <- epiABM(type='SIR', s.num=500, i.num=1, beta=0.2, cont=0.25, nu=1/50, nsteps=250, nsims=10)
plotABM(out, 'i.num', medline=FALSE)

# SIS model
out <- epiABM(type='SIS', s.num=500, i.num=1, beta=0.2, cont=0.25, nu=1/50, nsteps=250, nsims=10)
plotABM(out, 'i.num', medline=FALSE)
}
\author{
  Samuel M. Jenness <sjenness@uw.edu> and Steven M.
  Goodreau <goodreau@uw.edu>
}
\seealso{
  \code{\link{plotABM}}
}
\keyword{model}

