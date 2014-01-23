\name{epiNet.simTrans}
\alias{epiNet.simTrans}
\title{Stochastic Network Epidemic Models}
\usage{
epiNet.simTrans(x, type, vital = FALSE, i.num, i.ids, i.rand = TRUE,
  trans.rate, trans.rate.m2, act.rate = 1, rec.rate, b.rate, ds.rate, di.rate,
  dr.rate, nsteps, sims.per.nw = 1, tea = FALSE, save.statmat = TRUE,
  save.stats = FALSE, stats.formula, save.trans = TRUE,
  save.network = TRUE, verbose = TRUE, plot.prog = FALSE, ...)
}
\arguments{
  \item{x}{an \code{EpiModel} object of class
  \code{epiNet.est} for dependent disease-network
  simulations or class \code{epiNet.simNet} for independent
  disease-network simulations (see details).}

  \item{type}{disease type to be modeled, with the choice
  of \code{"SI"} for Susceptible-Infected diseases,
  \code{"SIR"} for Susceptible-Infected-Recovered diseases,
  and \code{"SIS"} for Susceptible-Infected-Susceptible
  diseases.}

  \item{vital}{if \code{TRUE}, incorporate vital dynamics
  (births and deaths) into the model. This requires passing
  a birth rate and death rate for each compartment in the
  model (see details).}

  \item{i.num}{number in the population initially infected.
  In bipartite models, infections are randomly allocated
  across modes with equal probability.}

  \item{i.ids}{a vector of node ID numbers to set those
  specific nodes as initially infected at baseline. One
  must use either \code{i.num} or \code{i.ids} to set
  initial infected. Setting infected IDs here overrides
  \code{i.rand=TRUE}.}

  \item{i.rand}{if \code{TRUE} and using \code{i.num}, sets
  infection based on random binomial draws from
  distribution implied by i.num.}

  \item{trans.rate}{probability of transmission given an
  act or contact between a susceptible and an infected
  person in the population. In bipartite models this is the
  rate governing the infection to the mode 1 members.}

  \item{trans.rate.m2}{probability of transmission given an
  act or contact between a susceptible mode 2 person and an
  infected mode 1 person in the population. It is the rate
  governing the infection to the mode 2 members.}

  \item{act.rate}{average number of acts governing
  transmission \emph{per partnership} per unit time,
  regardless of disease status (see details).}

  \item{rec.rate}{average rate of recovery with immunity
  (in \code{SIR} models) or re-susceptibility (in
  \code{SIS} models). The recovery rate is a function of
  the disease duration as 1/duration of disease. The
  recovery rate is a function of the disease duration as
  1/duration of disease.}

  \item{b.rate}{birth rate into the network. For one-mode
  models, the birth rate is the probability of a new birth
  per person per unit time. For bipartite models, the birth
  rate is the rate per "mode 1 persons" per unit time.}

  \item{ds.rate}{mortality rate for susceptibles, i.e. the
  probability of death per susceptible person per unit
  time.}

  \item{di.rate}{mortality rate for infecteds, i.e., the
  probability of death per infected person per unit time.}

  \item{dr.rate}{mortality rate for recovered, i.e., the
  probability of death per recovered person per unit time.}

  \item{nsteps}{number of time steps to simulate the model
  over. This must be a positive integer. In independent
  models (\code{vital=FALSE}), \code{nsteps} is set by
  default through the \code{\link{epiNet.simNet}} network
  simulation, although \code{nsteps} may be set here to a
  smaller number of steps. For dependent models
  (\code{vital=FALSE}), \code{nsteps} must be set since
  there is no default.}

  \item{sims.per.nw}{In independent models
  (\code{vital=FALSE}), the number of disease simulations
  per network simulated in \code{\link{epiNet.simNet}}.  In
  dependent models, (\code{vital=TRUE}), the total number
  of disease simulations, since the network is resimulated
  at each time step.}

  \item{tea}{if \code{TRUE}, use temporally extended
  attributes to store disease status information. This
  tends to slow the simulations but allows for more
  flexible network plotting and animations.}

  \item{save.statmat}{if \code{TRUE} and \code{tea=FALSE},
  save a disease status matrix to the main output object.
  This is a time by n sized matrix containing the disease
  status of all nodes at each time; it may be used instead
  of temporally extended attributes to store all disease
  history for later plotting and analysis.}

  \item{save.stats}{if \code{TRUE} and \code{vital=TRUE},
  save network statistics in a \code{data.frame} out to the
  main output object. The specific statistics to be
  requested are set in the \code{stats.formula} parameter.
  Note that in independent models, \code{save.stats=TRUE}
  is overridden because the networks are their statistics
  of interest have already been generated during the
  \code{\link{epiNet.simNet}} simulation.}

  \item{stats.formula}{formula to specify diagnostic
  statistics, where the default is formation formula
  (specified during the estimation phase).}

  \item{save.trans}{if \code{TRUE}, save a transmission
  \code{data.frame} for each simulation out to the main
  object. This object contains one row for each
  transmission event that has occurred (see
  \code{\link{discord.edgelist}}).}

  \item{save.network}{if \code{TRUE}, save a
  \code{networkDynamic} object containing full edge and
  nodal history for each simulation out to the main
  object.}

  \item{verbose}{if \code{TRUE}, print simulation progress
  to the console.}

  \item{plot.prog}{if \code{TRUE}, dynamically plot
  simulation progress to the plot window.}

  \item{...}{additional arguments to pass to
  \code{simulate.stergm}.}
}
\value{
This function returns a list object of class
\code{epiNet.simTrans} with the following elements:
\itemize{ \item \strong{type:} disease type as specified in
model parameter.  \item \strong{modes:} network modes as
set in the starting network for estimation in the
\code{\link{epiNet.est}} call.  \item \strong{time:} a
vector of time steps over which the model was solved.
\item \strong{nruns:} number of independent model runs.
\item \strong{*.num:} a vector (if \code{nruns=1}) or data
frame (if \code{nruns>1}) of compartment or state sizes
over time, for each model run, where * may equal \code{s},
\code{i}, or \code{r} for susceptible, infected, or
recovered, respectively. Bipartite models have equivalent
elements with a \code{.m2} suffix.  \item \strong{*.flow:}
a vector (if \code{nruns=1}) or data frame (if
\code{nruns>1}) of transition sizes between compartments,
for each model run, where * may equal \code{si} for
susceptible to infected flows, \code{is} for infected to
susceptible flows, \code{ir} for infected to recovered
flows, \code{b} for birth in-flows, \code{ds} for
susceptible death out-flows, \code{di} for infected death
out-flows, \code{dr} for recovered death out-flows.
Bipartite models have equivalent elements with a \code{.m2}
suffix.  \item \strong{stat.mat:} a list of matrices, one
for each simulation, containing the disease status history,
specified to be saved with the \code{save.statmat}
parameter.  \item \strong{trans:} list of data frames, one
for each simulation, containing the transmission history,
specified to be saved with the \code{save.trans} parameter.
\item \strong{network:} a list of \code{networkDynamic}
objects, one for each simulation, containing the complete
edge and nodal histories, specified to be saved with the
\code{save.network} parameter.  \item \strong{vital:}
logical, whether vital dynamics were specified in the
parameterization.  \item \strong{nsims:} number of disease
simulations conducted.  \item \strong{call:} exact model
call. }
}
\description{
This function simulates stochastic network-based epidemic
models for infectious disease.
}
\details{


This is the main disease simulation function for stochastic
network-based epidemic models in \code{EpiModel}. This
class of models differs from both deterministic
compartmental models in \code{epiDCM} and their stochastic
counterparts in \code{epiICM} by explicitly tracking
partnership dyads and changes within those dyads over time.
All three disease types (SI, SIR, and SIS) may be modeled
in \code{epiNet.simTrans}.

If one is interested in extended these network models for
novel research, it will be necessary to understand what
goes on "under-the-hood" of this simulation function. To
start, it is recommended that you review the
\code{\link{epiNetModules}} help page, which gives an
overview of all the functions that update the disease
simulation at each time step. Since these functions are not
called directly by the end-user, they are not shown on the
main help index.
}
\section{Dependent versus Independent Models}{
  One important modelling choice is between a dependent and
  independent disease model. In an independent model, the
  disease simulation does not change the network structure,
  but in a dependent model demographics and disease
  influence the network structure. See Section 4 of
  \href{../doc/Tutorial.pdf}{EpiModel Tutorial} for further
  clarification.

  In this current version of \code{EpiModel}, dependent
  models is synonymous with incorporating vital dynamics
  (births and deaths), however that is only one type of
  dependence, and other versions will be incorporated in
  future \code{EpiModel} releases. One incorporates vital
  dynamics into an \code{epiNet.simTrans} simulation by
  setting \code{vital=TRUE}. In these dependent
  simulations, the \code{x} parameter is the network
  estimation object from \code{\link{epiNet.est}}, whereas
  in independent simulations the \code{x} parameter is the
  network simulation object from
  \code{\link{epiNet.simNet}}. With the latter, the dynamic
  network is simulated in advance and outside of this
  disease simulation; with the former, the network is
  resimulated at each timestep, and only the model fit is
  needed to get started.
}

\section{The act.rate Parameter}{
  A key difference between these network models and the
  DCM/ICM models is their treatment of transmission events.
  In DCM and ICM, contacts or partnerships are
  mathematically instantaneous events: they have no
  duration in time, and thus no changes may occur within
  them over time. In contrast, network models allow for
  partnership durations defined by the dynamic network
  model, summarized in the model dissolution coefficients
  calculated in \code{\link{dissolution.coefs}}. Therefore,
  the \code{act.rate} parameter has a different
  interpretation here, where it is the number of acts
  defining transmission per partnership per unit time.
}
\examples{
\dontrun{
## See the EpiModel Tutorial vignette for more details and examples ##

## Independent SI Model ##
# Initialize network and set network model parameters
nw <- network.initialize(n=500, bipartite=250, directed=FALSE)
formation <- ~ edges + b1degree(0:1) + b2degree(0:1)
target.stats <- c(165, 100, 137.5, 120, 102.5)
dissolution <- ~ offset(edges)
coef.diss <- dissolution.coefs(dissolution, duration=25)

# Estimate the ERGM models (see help for epiNet.est)
# Skipping model diagnostics for this, but should always run
est1 <- epiNet.est(
  nw,
  formation,
  dissolution,
  target.stats,
  coef.diss)

# Simulate 10 networks from model fit (see help for epiNet.simNet)
# Skipping model diagnostics for this, but should always run
nwsims <- epiNet.simNet(est1, nsteps = 250, nsims = 10)

# Independent epidemic model simulation for SI disease
# Transmissibility is higher to mode 1 than mode 2
sim.SI <- epiNet.simTrans(
  nwsims,
  type = "SI",
  vital = FALSE,
  trans.rate = 0.3,
  trans.rate.m2 = 0.15,
  i.num = 50,
  sims.per.nw = 1)

# Print, plot, and summarize the results
sim.SI
plot(sim.SI, type="sim")
plot(sim.SI, y=c("i.num", "i.num.m2"), sim.lines=TRUE,
     sim.col=c("steelblue", "firebrick"), qnts=0.5)
plot(sim.SI, type="network", at=50, col.inf=TRUE, shp.bip="square")
summary(sim.SI, time=50)

## Dependent SIR Model ##
coef.diss <- dissolution.coefs(dissolution, duration=25, d.rate=0.0021)
est2 <- epiNet.est(
  nw,
  formation,
  dissolution,
  target.stats,
  coef.diss)

sim.SIR <- epiNet.simTrans(
  est2,
  type = "SIR",
  vital = TRUE,
  i.num = 50,
  trans.rate = 0.3,
  trans.rate.m2 = 0.15,
  rec.rate = 0.02,
  b.rate = 0.0021*2,
  ds.rate = 0.0021,
  di.rate = 0.0021*1.1,
  dr.rate = 0.0021,
  nsteps = 250,
  sims.per.nw = 10)

# Print, plot, and summarize the results
sim.SIR
plot(sim.SIR, type="sim")
plot(sim.SIR, type="network")
summary(sim.SIR, time=100)
}
}
\references{
Goodreau SM, Carnegie NB, et al. What drives the US and
Peruvian HIV epidemics in men who have sex with men (MSM)?
PloS One. 2012; 7(11): e50522.
}
\seealso{
\code{\link{epiNet.est}} for network estimation,
\code{\link{epiNet.simNet}} for network simulation,
\code{\link{plot.epiNet.simTrans}} for plots, and
\code{\link{summary.epiNet.simTrans}} for data summaries.
}
\keyword{model}

