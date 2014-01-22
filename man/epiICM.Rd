\name{epiICM}
\alias{epiICM}
\title{Individual Contact Epidemic Models}
\usage{
epiICM(type, groups = 1, s.num, i.num, r.num = 0, trans.rate, act.rate,
  rec.rate = 0, b.rate = 0, ds.rate = 0, di.rate = 0, dr.rate = 0,
  s.num.g2, i.num.g2, r.num.g2 = 0, trans.rate.g2, act.rate.g2,
  rec.rate.g2 = 0, b.rate.g2 = 0, ds.rate.g2 = 0, di.rate.g2 = 0,
  dr.rate.g2 = 0, balance, nsteps, rec.rand = TRUE, b.rand = TRUE,
  d.rand = TRUE, nsims = 1, verbose = TRUE)
}
\arguments{
  \item{type}{disease type to be modeled, with the choice
  of \code{"SI"} for Susceptible-Infected diseases,
  \code{"SIR"} for Susceptible-Infected-Recovered diseases,
  and \code{"SIS"} for Susceptible-Infected-Susceptible
  diseases.}

  \item{groups}{number of mixing groups to model, with the
  default of \code{1} for purely random mixing in the
  population, or \code{2} for purely heterogenous mixing
  (e.g., two sexes with purely heterosexual contact).}

  \item{s.num}{number of initial susceptibles in the
  population. For two-group models, this is the number of
  initial group 1 susceptible in the population.}

  \item{i.num}{number of initial infecteds in the
  population. For two-group models, this is the number of
  initial group 1 infected in the population.}

  \item{r.num}{number of initial recovereds in the
  population. For two-group models, this is the number of
  initial group 1 recovered in the populations. This
  parameter is only used for the \code{SIR} model type.}

  \item{trans.rate}{probability of transmission given an
  act or contact between a susceptible and an infected
  person in the population. In two-group models this is the
  rate governing the probability of transmission to the
  group 1 members.}

  \item{act.rate}{average number of acts governing
  transmission per person per unit time, regardless of
  disease status. For two-group models, this is the number
  of acts per group 1 persons per unit time; note that a
  balance between the acts in groups 1 and 2 is necessary,
  and set using the \code{balance} parameter.}

  \item{rec.rate}{average rate of recovery with immunity
  (in \code{SIR} models) or re-susceptibility (in
  \code{SIS} models). The recovery rate is a function of
  the disease duration as 1/duration of disease. For
  two-group models, this is the recovery rate for group 1
  persons only. This parameter is only used for \code{SIR}
  and \code{SIS} models.}

  \item{b.rate}{birth rate into the population. For
  one-group models, the birth rate is the probability of a
  new birth per person per unit time. For two-group models,
  the birth rate may be parameterized as a rate per group 1
  person per unit time (with group 1 persons presumably
  representing females), and with the \code{b.rate.g2} rate
  set as described below.}

  \item{ds.rate}{mortality rate for susceptibles, i.e. the
  probability of death per susceptible person per unit
  time. For two-group models, it is the rate for the group
  1 susceptibles only.}

  \item{di.rate}{mortality rate for infecteds, i.e., the
  probability of death per infected person per unit time.
  For two-group models, it is the rate for the group 1
  infecteds only.}

  \item{dr.rate}{mortality rate for recovereds, i.e., the
  probability of death per recovered person per unit time.
  For two-group models, it is the rate for the group 1
  recovered only. This parameter is only used for
  \code{SIR} models.}

  \item{s.num.g2}{number of initial susceptibles in group 2
  in the population. This parameter is only used for
  two-group models.}

  \item{i.num.g2}{number of initial infecteds in group 2 in
  the population. This parameter is only used for two-group
  models.}

  \item{r.num.g2}{number of initial recovereds in group 2
  in the population. This parameter is only used for
  two-group \code{SIR} models.}

  \item{trans.rate.g2}{probability of transmission given an
  act or contact between a susceptible group 2 person and
  an infected group 1 person in the population. It is the
  rate governing the probability of transmission to group 2
  members.}

  \item{act.rate.g2}{average number of acts capable of
  transmission per group 2 person per unit time. Note that
  a balance between the acts in groups 1 and 2 is
  necessary, and set using the \code{balance} parameter.}

  \item{rec.rate.g2}{average rate of recovery with immunity
  (in \code{SIR} models) or re-susceptibility (in
  \code{SIS} models) for group 2 persons. This parameter is
  only used for two-group \code{SIR} and \code{SIS}
  models.}

  \item{b.rate.g2}{birth rate into the population for group
  2. This may either be specified numerically as the number
  of new births per group 2 persons per unit time, or as
  \code{NA} in which case the group 1 rate, \code{b.rate},
  governs the birth process. The latter is used in cases
  where the first group is conceptualized as female, and
  the group 1 population size determines the birth rate.
  Currently births are evenly allocated to the two groups.}

  \item{ds.rate.g2}{mortality rate for group 2
  susceptibles, i.e., the probability of death per group 2
  susceptible person per unit time.}

  \item{di.rate.g2}{mortality rate for group 2 infecteds,
  i.e., the probability of death per group 2 infected
  person per unit time.}

  \item{dr.rate.g2}{mortality rate for group 2 recovereds,
  i.e., the probability of death per group 2 recovered
  person per unit time. This parameter is only used for
  \code{SIR} model types.}

  \item{balance}{for two-group models, balance the
  \code{act.rate} to the rate set for group 1 (with
  \code{balance="g1"}) or group 2 (with
  \code{balance="g2"}). Further details provided below.}

  \item{nsteps}{number of time steps to solve the model
  over. This must be a positive integer.}

  \item{rec.rand}{if \code{TRUE}, use a stochastic recovery
  model, with the number of recovered a random draw from a
  Poisson distribution for the product of the recovery rate
  and the number infected. If \code{FALSE}, then a
  deterministic rounded count of that product.}

  \item{b.rand}{if \code{TRUE}, use a stochastic birth
  model, with the number of new births a random draw from a
  Poisson distribution for the product of the birth rate
  and population size. If \code{FALSE}, then a
  deterministic rounded count of that product.}

  \item{d.rand}{if \code{TRUE}, use a stochastic mortality
  model, with the number of new deaths a random draw from a
  Poisson distribution for the product of the
  state-specific death rates and the state size. If
  \code{FALSE}, then a deterministic rounded count of those
  products.}

  \item{nsims}{number of simulations to run.}

  \item{verbose}{print model simulation progress to the
  console.}
}
\value{
This function returns a list object of class \code{epiICM}
with the following elements: \itemize{ \item \strong{type:}
disease type as specified in model parameter.  \item
\strong{groups:} groups as specified in model parameter.
\item \strong{params:} list of model parameters as
specified in model parameter.  \item \strong{time:} a
vector of time steps over which the model was solved.
\item \strong{vital:} logical, whether vital dynamics were
specified in the parameterization.  \item \strong{nsims:}
number of model simulations.  \item \strong{call:} exact
model call.  \item \strong{*.num:} a vector (if
\code{nruns=1}) or data frame (if \code{nruns>1}) of
compartment or state sizes over time, for each model run,
where * may equal \code{s}, \code{i}, or \code{r} for
susceptible, infected, or recovered, respectively. Two
group models have equivalent elements with a \code{.g2}
suffix.  \item \strong{*.flow:} a vector (if
\code{nruns=1}) or data frame (if \code{nruns>1}) of
transition sizes between compartments, for each model run,
where * may equal \code{si} for susceptible to infected
flows, \code{is} for infected to susceptible flows,
\code{ir} for infected to recovered flows, \code{b} for
birth in-flows, \code{ds} for susceptible death out-flows,
\code{di} for infected death out-flows, \code{dr} for
recovered death out-flows.  Two group models have
equivalent elements with a \code{.g2} suffix. }
}
\description{
This function simulates individual contact transmission
models for infectious disease.
}
\details{


The \code{epiICM} function is a stochastic, discrete-time
representation of the deterministic compartmental disease
models that may be solved with \code{epiDCM}. Random
variation is simulated in all transition processes:
infection, recovery, and vital dynamics. The stochasticity
in each transition except infection may be toggled on or
off using the control parameters \code{rec.rand},
\code{b.rand}, and \code{d.rand}.

One important caveat is the current parameter structure of
specifing initial state sizes and individual parameters as
primary \code{epiICM} parameters will be deprecated in the
\code{EpiModel} v1.0 release. More flexibile passing of
model behavior and parameters will be implemented.
}
\section{The act.rate Parameter}{
  The \code{act.rate} parameter represents the average
  number of acts capable of transmission per person per
  unit time. In the sexually transmitted disease modeling
  literature, the term ``contact'' if often used, as in the
  standard formulation \eqn{R_0=\beta c D}. But contact may
  either mean the number of independent acts or number of
  partnerships per person per unit time. In the latter
  case, one often obtains the risk of transmission within a
  partnership per unit time by solving the formula
  \eqn{\tau = 1-(1-\alpha)^n}, where \eqn{\tau} is the
  probability of transmission per partnership, \eqn{\alpha}
  is the probability of transmission per act, and \eqn{n}
  is the number of acts per partnership. \code{act.rate}
  parameter is to be distinguished from ``effective
  contact,'' which here is the product of \code{trans.rate}
  and \code{act.rate}.

  In two-group models, it is necessary to maintain a
  balance between the number of acts for group 1 members
  and those for group 2 members. Mathematically in a purely
  heterogenous mixing framework, the product of one group
  size and act rate must equal the product of the other
  group size and act rate: \eqn{N_1 \alpha_1 = N_2
  \alpha_2}, where \eqn{N_i} is the group size and
  \eqn{\alpha_i} the group-specific act rates at time
  \eqn{t}.
}
\examples{
\dontrun{
# SI Model
mod <- epiICM(type="SI", s.num=500, i.num=1,
              trans.rate=0.2, act.rate=0.25, nsteps=500, nsims=10)
plot(mod)

# SIR Model
mod <- epiICM(type="SIR", s.num=500, i.num=1,
              trans.rate=0.2, act.rate=0.25, rec.rate=1/50,
              nsteps=500, nsims=10)
plot(mod, y="r.num", sim.lines=TRUE, sim.alpha=0.25)

# SIS Model
mod <- epiICM(type="SIS", s.num=500, i.num=1,
              trans.rate=0.2, act.rate=0.25, rec.rate=1/50,
              nsteps=500, nsims=10)
plot(mod, sim.lines=T, qnts=FALSE, sim.col=c("steelblue", "red"),
     leg.cex=0.75, popfrac=FALSE)
}
}
\seealso{
\code{\link{plot.epiICM}}, \code{\link{summary.epiICM}},
\code{\link{as.data.frame.epiICM}}
}
\keyword{model}

