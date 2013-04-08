\name{epiDCM}
\alias{epiDCM}
\title{Deterministic Compartmental SI-SIS-SIR Models}
\usage{
  epiDCM(type, s.num, i.num, r.num = 0, beta, cont, nu = 0,
    b = 0, ms = 0, mi = 0, mr = 0, dt, new.mod = NULL,
    print.mod = FALSE, verbose = FALSE, ...)
}
\arguments{
  \item{type}{model type (choice of 'SI', 'SIR', or 'SIS')}

  \item{s.num}{number of initial susceptible in population}

  \item{i.num}{number of initial infected in population}

  \item{r.num}{number of initial recovered in population}

  \item{beta}{probability of infection per contact}

  \item{cont}{average number of contacts per person per
  unit time}

  \item{nu}{average rate of recovery (1/duration of
  disease)}

  \item{b}{birth rate}

  \item{ms}{mortality rate for susceptibles}

  \item{mi}{mortality rate for infecteds}

  \item{mr}{mortality rate for recovered}

  \item{dt}{either scalar or vector of timesteps for model}

  \item{new.mod}{if not \code{NULL}, insert a new model
  definition}

  \item{print.mod}{extract and print the current model from
  the function}

  \item{verbose}{print progress}

  \item{...}{additional arguments to pass to model}
}
\description{
  This function solves a deterministic model of an
  infectious disease with random mixing in the population.
}
\examples{
# SI model
out <- epiDCM(type='SI', s.num=500, i.num=1, beta=0.2, cont=0.25, dt=500)
par(mar=c(3.5,3,1,1), mgp=c(2,1,0))
plotDCM(out, c('s.num', 'i.num'), leg='full')

# SIR model
out <- epiDCM(type='SIR', s.num=500, i.num=1, r.num=0,
              beta=0.2, cont=0.25, nu=1/50, dt=500)
plotDCM(out, c('s.num', 'i.num', 'r.num'), leg='full')

# SIS model
out <- epiDCM(type='SIS', s.num=500, i.num=1, r.num=0,
              beta=0.2, cont=0.25, nu=1/50, dt=500)
plotDCM(out, c('s.num', 'i.num'), leg='full')

# SIR model with vital dynamics
out <- epiDCM(type='SIR', s.num=1000, i.num=1, r.num=0,
              beta=0.2, cont=5, nu=1/3,
              b=1/90, ms=1/100, mi=1/35, mr=1/100, dt=500)
plotDCM(out, c('s.num', 'i.num', 'r.num'), leg='full')

# SIS model with sensitivity parameter for contact rate
out <- epiDCM(type='SIS', s.num=500, i.num=1, r.num=0,
              beta=0.2, cont=seq(0.05,0.5,0.05),
              nu=1/50, dt=500, verbose=TRUE)
plotDCM(out, 's.num', leg='lim', alpha=0.6)
}
\author{
  Samuel M. Jenness <sjenness@uw.edu>
}
\seealso{
  \code{\link{plotDCM}}
}
\keyword{model}

