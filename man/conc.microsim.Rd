\name{conc.microsim}
\alias{conc.microsim}
\title{Concurrency Microsimulation Model for HIV-1 Transmission Dynamics}
\usage{
conc.microsim(s.num.f, i.num.f, s.num.m, i.num.m, monog.f = TRUE,
  monog.m = TRUE, meandeg, part.duration, nsteps, nsims = 1,
  verbose = TRUE)
}
\arguments{
  \item{s.num.f}{number of initial susceptible females in
  the population.}

  \item{i.num.f}{number of initial infected females in the
  population.}

  \item{s.num.m}{number of initial susceptible males in the
  population.}

  \item{i.num.m}{number of initial infected females in the
  population.}

  \item{monog.f}{if \code{TRUE}, enforce a momentary degree
  constraint of monogamy for females (females not allowed
  concurrent partnerships).}

  \item{monog.m}{if \code{TRUE}, enforce a momentary degree
  constraint of monogamy for males (males not allowed
  concurrent partnerships).}

  \item{meandeg}{average momentary mean degree (number of
  current partnerships) in the population.}

  \item{part.duration}{average length of partnerships in
  months.}

  \item{nsteps}{number of time steps to simulate the model
  over. This must be a positive integer.}

  \item{nsims}{number of simulations to run.}

  \item{verbose}{print model simulation progress to the
  console.}
}
\description{
This function simulates an HIV epidemic in a population of
men and women with purely heterosexual mixing under varying
scenarios of sexual partnership concurrency.
}
\details{


This is a microsimulation model of a dynamic network to see
how the presence or absence of relational concurrency
affects the prevalence of HIV prevalence at the population
level.

HIV infection is simulated based on a four-stage disease
progression model in which persons transition from acute to
latent to pre-AIDS to AIDS stages. These transitions occur
at deterministic intervals based on estimates of the
average time per stage. Also, the transmission probability
to uninfected partners varies by stage of the infected
partner: it is highest in the acute stage and lowest in the
AIDS stage when no sexual acts are simulated to occur. See
the Hollingsworth reference below for further details.

The main parameters of the model include the initial number
of susceptible females and males, the initial number of
infected females and males, whether males and females are
allowed concurrency, the mean degree of all persons in the
population, and the average duration of partnerships (in
momths). As the four examples below show, we can hold
constant all parameters but toggle whether men, women, or
both sexes exhibit concurrency.

This model makes several simplifying assumptions about
partnership formation and dissolution, such as the
phenomenon of all dissolving partnerships being immediately
replaced by a new partnership in the network. Additionally,
the user may specify whether concurrency is allowed, but
not the level of concurrency (it is automatically
calculated here based on a binomial distribution with the
probability set by the mean degree parameter; see the full
tutorial link below). Therefore, this model as an
introduction to network modeling featured in the
\code{epiNet} class of functions in \code{EpiModel}: there
the user has much more control over the network
parameterization and evolution.
}
\examples{
\dontrun{
# No concurrency model
no.conc <- conc.microsim(
   s.num.f = 1000,
   i.num.f = 50,
   s.num.m = 1000,
   i.num.m = 50,
   monog.f = TRUE,
   monog.m = TRUE,
   meandeg = 0.8,
   part.duration = 10,
   nsteps = 2000,
   nsims = 10,
   verbose = TRUE)

# Male concurrency only model
male.conc <- conc.microsim(
   s.num.f = 1000,
   i.num.f = 50,
   s.num.m = 1000,
   i.num.m = 50,
   monog.f = TRUE,
   monog.m = FALSE,
   meandeg = 0.8,
   part.duration = 10,
   nsteps = 2000,
   nsims = 10,
   verbose = TRUE)

# Female concurrency only model
feml.conc <- conc.microsim(
   s.num.f = 1000,
   i.num.f = 50,
   s.num.m = 1000,
   i.num.m = 50,
   monog.f = FALSE,
   monog.m = TRUE,
   meandeg = 0.8,
   part.duration = 10,
   nsteps = 2000,
   nsims = 10,
   verbose = TRUE)

# Both sexes concurrency model
both.conc <- conc.microsim(
   s.num.f = 1000,
   i.num.f = 50,
   s.num.m = 1000,
   i.num.m = 50,
   monog.f = FALSE,
   monog.m = FALSE,
   meandeg = 0.8,
   part.duration = 10,
   nsteps = 2000,
   nsims = 10,
   verbose = TRUE)

# Plot the results
par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(2,1,0))
plot(no.conc, alpha=0.5, ylim=c(0, 0.5), main="No Concurrency")
plot(male.conc, alpha=0.5, ylim=c(0, 0.5), main="Male Concurrency")
plot(feml.conc, alpha=0.5, ylim=c(0, 0.5), main="Female Concurrency")
plot(both.conc, alpha=0.5, ylim=c(0, 0.5), main="Both Concurrency")
}
}
\references{
A web-based implementation of this model is available at
\url{http://statnet.org/apps/Conc}. The background and
details of this model are explained in a full tutorial on
concurrency at
\url{https://statnet.csde.washington.edu/trac/wiki/ConcurrencyIndex}.

Hollingsworth TD, Anderson RM, Fraser C. HIV-1
transmission, by stage of infection. Journal of Infectious
Diseases. 2008; 198(5): 687-693.
}
\keyword{model}

