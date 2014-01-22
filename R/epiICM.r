#' 
#' @title Individual Contact Epidemic Models 
#'
#' @description This function simulates individual contact transmission models
#'   for infectious disease.
#'
#' @param type disease type to be modeled, with the choice of \code{"SI"} for
#'   Susceptible-Infected diseases, \code{"SIR"} for Susceptible-Infected-Recovered
#'   diseases, and \code{"SIS"} for Susceptible-Infected-Susceptible diseases.
#' @param groups number of mixing groups to model, with the default of \code{1} for
#'   purely random mixing in the population, or \code{2} for purely heterogenous
#'   mixing (e.g., two sexes with purely heterosexual contact).
#' @param s.num number of initial susceptibles in the population. For two-group
#'   models, this is the number of initial group 1 susceptible in the population.
#' @param i.num number of initial infecteds in the population. For two-group models,
#'   this is the number of initial group 1 infected in the population.
#' @param r.num number of initial recovereds in the population. For two-group
#'   models, this is the number of initial group 1 recovered in the populations. This
#'   parameter is only used for the \code{SIR} model type.
#' @param trans.rate probability of transmission given an act or contact between
#'   a susceptible and an infected person in the population. In two-group models
#'   this is the rate governing the probability of transmission to the group 1 
#'   members.
#' @param act.rate average number of acts governing transmission per person  
#'   per unit time, regardless of disease status. For two-group models, this
#'   is the number of acts per group 1 persons per unit time; note that a 
#'   balance between the acts in groups 1 and 2 is necessary, and set using
#'   the \code{balance} parameter.
#' @param rec.rate average rate of recovery with immunity (in \code{SIR} models) 
#'   or re-susceptibility (in \code{SIS} models). The recovery rate is a function 
#'   of the disease duration as 1/duration of disease. For two-group models, 
#'   this is the recovery rate for group 1 persons only. This parameter is only 
#'   used for \code{SIR} and \code{SIS} models.
#' @param b.rate birth rate into the population. For one-group models, the birth 
#'   rate is the probability of a new birth per person per unit time. For 
#'   two-group models, the birth rate may be parameterized as a rate per group 1 
#'   person per unit time (with group 1 persons presumably representing females), 
#'   and with the \code{b.rate.g2} rate set as described below.
#' @param ds.rate mortality rate for susceptibles, i.e. the probability of death 
#'   per susceptible person per unit time. For two-group models, it is the rate 
#'   for the group 1 susceptibles only.
#' @param di.rate mortality rate for infecteds, i.e., the probability of death per 
#'   infected person per unit time. For two-group models, it is the rate for the 
#'   group 1 infecteds only.
#' @param dr.rate mortality rate for recovereds, i.e., the probability of death 
#'   per recovered person per unit time. For two-group models, it is the rate for 
#'   the group 1 recovered only. This parameter is only used for \code{SIR} models.
#' @param s.num.g2 number of initial susceptibles in group 2 in the population. This
#'   parameter is only used for two-group models. 
#' @param i.num.g2 number of initial infecteds in group 2 in the population. This
#'   parameter is only used for two-group models. 
#' @param r.num.g2 number of initial recovereds in group 2 in the population. This
#'   parameter is only used for two-group \code{SIR} models. 
#' @param trans.rate.g2 probability of transmission given an act or contact between
#'   a susceptible group 2 person and an infected group 1 person in the 
#'   population. It is the rate governing the probability of transmission to 
#'   group 2 members.
#' @param act.rate.g2 average number of acts capable of transmission 
#'   per group 2 person per unit time. Note that a balance between the acts in 
#'   groups 1 and 2 is necessary, and set using the \code{balance} parameter.
#' @param rec.rate.g2 average rate of recovery with immunity (in \code{SIR} models) 
#'   or re-susceptibility (in \code{SIS} models) for group 2 persons. This 
#'   parameter is only used for two-group \code{SIR} and \code{SIS} models.
#' @param b.rate.g2 birth rate into the population for group 2. This may either 
#'   be specified numerically as the number of new births per group 2 persons per 
#'   unit time, or as \code{NA} in which case the group 1 rate, \code{b.rate}, 
#'   governs the birth process. The latter is used in cases where the first group 
#'   is conceptualized as female, and the group 1 population size determines the 
#'   birth rate. Currently births are evenly allocated to the two groups.
#' @param ds.rate.g2 mortality rate for group 2 susceptibles, i.e., the probability 
#'   of death per group 2 susceptible person per unit time.
#' @param di.rate.g2 mortality rate for group 2 infecteds, i.e., the probability 
#'   of death per group 2 infected person per unit time.
#' @param dr.rate.g2 mortality rate for group 2 recovereds, i.e., the probability 
#'   of death per group 2 recovered person per unit time. This parameter is only 
#'   used for \code{SIR} model types. 
#' @param balance for two-group models, balance the \code{act.rate} to the rate 
#'   set for group 1 (with \code{balance="g1"}) or group 2 (with 
#'   \code{balance="g2"}). Further details provided below. 
#' @param nsteps number of time steps to solve the model over. This must be a 
#'   positive integer.
#' @param rec.rand if \code{TRUE}, use a stochastic recovery model, with the number 
#'   of recovered a random draw from a Poisson distribution for the product of the 
#'   recovery rate and the number infected. If \code{FALSE}, then a deterministic 
#'   rounded count of that product.
#' @param b.rand if \code{TRUE}, use a stochastic birth model, with the number of 
#'   new births a random draw from a Poisson distribution for the product of the 
#'   birth rate and population size. If \code{FALSE}, then a deterministic 
#'   rounded count of that product.
#' @param d.rand if \code{TRUE}, use a stochastic mortality model, with the number 
#'   of new deaths a random draw from a Poisson distribution for the product of 
#'   the state-specific death rates and the state size. If \code{FALSE}, then a 
#'   deterministic rounded count of those products.
#' @param nsims number of simulations to run.
#' @param verbose print model simulation progress to the console.
#' 
#' @details 
#' The \code{epiICM} function is a stochastic, discrete-time representation of the 
#' deterministic compartmental disease models that may be solved with \code{epiDCM}. 
#' Random variation is simulated in all transition processes: infection, recovery, 
#' and vital dynamics. The stochasticity in each transition except infection may 
#' be toggled on or off using the control parameters \code{rec.rand}, 
#' \code{b.rand}, and \code{d.rand}.
#' 
#' One important caveat is the current parameter structure of specifing initial state
#' sizes and individual parameters as primary \code{epiICM} parameters will be 
#' deprecated in the \code{EpiModel} v1.0 release. More flexibile passing of model
#' behavior and parameters will be implemented. 
#' 
#' @section The act.rate Parameter:
#' The \code{act.rate} parameter represents the average number of acts capable of 
#' transmission per person per unit time. In the sexually transmitted disease
#' modeling literature, the term ``contact'' if often used, as in the standard
#' formulation \eqn{R_0=\beta c D}. But contact may either mean the number of 
#' independent acts or number of partnerships per person per unit time. In the 
#' latter case, one often obtains the risk of transmission within a partnership
#' per unit time by solving the formula \eqn{\tau = 1-(1-\alpha)^n}, where 
#' \eqn{\tau} is the probability of transmission per partnership, \eqn{\alpha} 
#' is the probability of transmission per act, and \eqn{n} is the number of acts
#' per partnership. \code{act.rate} parameter is to be distinguished from 
#' ``effective contact,'' which here is the product of \code{trans.rate} and 
#' \code{act.rate}. 
#' 
#' In two-group models, it is necessary to maintain a balance between the number 
#' of acts for group 1 members and those for group 2 members. Mathematically in 
#' a purely heterogenous mixing framework, the product of one group size and act 
#' rate must equal the product of the other group size and act rate: \eqn{N_1 
#' \alpha_1 = N_2 \alpha_2}, where \eqn{N_i} is the group size and \eqn{\alpha_i} 
#' the group-specific act rates at time \eqn{t}.
#' 
#' @return
#' This function returns a list object of class \code{epiICM} with the following 
#' elements:
#' \itemize{
#'  \item \strong{type:} disease type as specified in model parameter.
#'  \item \strong{groups:} groups as specified in model parameter.
#'  \item \strong{params:} list of model parameters as specified in model parameter.
#'  \item \strong{time:} a vector of time steps over which the model was solved.
#'  \item \strong{vital:} logical, whether vital dynamics were specified in the 
#'  parameterization.
#'  \item \strong{nsims:} number of model simulations.
#'  \item \strong{call:} exact model call.
#'  \item \strong{*.num:} a vector (if \code{nruns=1}) or data frame (if 
#'    \code{nruns>1}) of compartment or state sizes over time, for each model run, 
#'    where * may equal \code{s}, \code{i}, or \code{r} for susceptible, infected, 
#'    or recovered, respectively. Two group models have equivalent elements with a
#'    \code{.g2} suffix.
#'  \item \strong{*.flow:} a vector (if \code{nruns=1}) or data frame (if 
#'    \code{nruns>1}) of transition sizes between compartments, for each model run, 
#'    where * may equal \code{si} for susceptible to infected flows, \code{is} for 
#'    infected to susceptible flows, \code{ir} for infected to recovered flows, 
#'    \code{b} for birth in-flows, \code{ds} for susceptible death out-flows, 
#'    \code{di} for infected death out-flows, \code{dr} for recovered death out-flows.
#'    Two group models have equivalent elements with a \code{.g2} suffix.
#' }
#' 
#' @keywords model
#' @seealso \code{\link{plot.epiICM}}, \code{\link{summary.epiICM}}, 
#'  \code{\link{as.data.frame.epiICM}}
#' @export
#' 
#' @examples
#' \dontrun{
#' # SI Model
#' mod <- epiICM(type="SI", s.num=500, i.num=1, 
#'               trans.rate=0.2, act.rate=0.25, nsteps=500, nsims=10)
#' plot(mod)
#' 
#' # SIR Model
#' mod <- epiICM(type="SIR", s.num=500, i.num=1, 
#'               trans.rate=0.2, act.rate=0.25, rec.rate=1/50, 
#'               nsteps=500, nsims=10)
#' plot(mod, y="r.num", sim.lines=TRUE, sim.alpha=0.25)
#' 
#' # SIS Model
#' mod <- epiICM(type="SIS", s.num=500, i.num=1, 
#'               trans.rate=0.2, act.rate=0.25, rec.rate=1/50, 
#'               nsteps=500, nsims=10)
#' plot(mod, sim.lines=T, qnts=FALSE, sim.col=c("steelblue", "red"), 
#'      leg.cex=0.75, popfrac=FALSE)
#' }
#' 
epiICM <- function(type, 
                   groups = 1,
                   s.num, 
                   i.num, 
                   r.num = 0, 
                   trans.rate, 
                   act.rate, 
                   rec.rate = 0, 
                   b.rate = 0, 
                   ds.rate = 0, 
                   di.rate = 0, 
                   dr.rate = 0,
                   s.num.g2, 
                   i.num.g2, 
                   r.num.g2 = 0,
                   trans.rate.g2, 
                   act.rate.g2, 
                   rec.rate.g2 = 0, 
                   b.rate.g2 = 0, 
                   ds.rate.g2 = 0, 
                   di.rate.g2 = 0, 
                   dr.rate.g2 = 0,
                   balance, 
                   nsteps, 
                   rec.rand = TRUE, 
                   b.rand = TRUE, 
                   d.rand = TRUE,
                   nsims = 1, 
                   verbose = TRUE) { 
  
  ## WARNINGS ##
  {
    # Model type check
    if ( missing(type) | !(type %in% c('SI','SIS','SIR')) ) {
      stop('\nAvailable model types are: SI SIS SIR')
    }
    
    # Group number check
    if ( !(groups %in% 1:2) ) {
      stop('\nOnly groups=1 and groups=2 are allowed')
    }
    
    # Balance check
    if ( groups == 2 && missing(balance) ) {
      stop("\nMust specify balance='m1' or balance='m2' with 2-group models")
    }
    
    if ( type == 'SI' && (rec.rate > 0 | rec.rate.g2 > 0) ) {
      stop('\nSet rec.rate=0 since this is an SI model')
    }
    
    # dr.rate warning for SI/SIS models
    if (type %in% c('SI', 'SIS') && (dr.rate > 0 | dr.rate.g2 > 0)) {
      warning('dr.rate parameter specified, but not used in SI/SIS models')
    }
  }
  
  if (verbose == TRUE) {
    if (nsims == 1) {
      cat('===============================')
      cat('\nStarting 1 Disease Simulation')
      cat('\n===============================\n')
    } else {
      cat('===============================')
      cat('\nStarting', nsims, 'Disease Simulations')
      cat('\n===============================\n')
    }
  }
  
  all.p <- as.list(match.call(expand.dots = TRUE)[-1]) 
  all.p <- all.p[-(grep('num', names(all.p)))]
  all.p <- all.p[-which(names(all.p) %in% c('type', 'groups', 
                                            'balance', 'rec.rand',
                                            'b.rand', 'd.rand',
                                            'nsteps', 'nsims', 
                                            'verbose'))]
  all.p <- sapply(all.p, function(x) eval(x))
  all.pout <- as.list(all.p)
  vd.params <- c('b.rate', 'ds.rate', 'di.rate', 'dr.rate',
                 'b.rate.g2', 'ds.rate.g2', 'di.rate.g2', 'dr.rate.g2')
  if (any(names(all.pout) %in% vd.params)) vital <- TRUE  else  vital <- FALSE
  
  # Compartment t0 values
  s.num.t0 <- s.num
  i.num.t0 <- i.num
  r.num.t0 <- r.num
  if (groups == 2) {
    s.num.t0.g2 <- s.num.g2
    i.num.t0.g2 <- i.num.g2
    r.num.t0.g2 <- r.num.g2
  }
  
  # Simulation loop start
  for (s in 1:nsims) {
    
    # Starting values
    {
      s.num <- s.num.t0
      i.num <- i.num.t0
      r.num <- r.num.t0
      num <- s.num + i.num + r.num
      ids <- 1:num
      status <- rep(0, num)
      status[ssample(ids, i.num.t0)] <- 1
      duration <- rep(NA, num)
      duration[status == 1] <- 0      
      si.flow <- ir.flow <- 0
      b.flow <- ds.flow <- di.flow <- dr.flow <- 0
      if (groups == 2) {
        s.num.g2 <- s.num.t0.g2
        i.num.g2 <- i.num.t0.g2
        r.num.g2 <- r.num.t0.g2
        num.g2 <- s.num.g2 + i.num.g2 + r.num.g2
        ids.g2 <- 1:num.g2
        status.g2 <- rep(0, num.g2)
        status.g2[ssample(ids.g2, i.num.t0.g2)] <- 1
        duration.g2 <- rep(NA, num.g2)
        duration.g2[status.g2 == 1] <- 0      
        si.flow.g2 <- ir.flow.g2 <- 0
        b.flow.g2 <- ds.flow.g2 <- di.flow.g2 <- dr.flow.g2 <- 0
      }
    }
    
    
    # Timestep loop
    for (ts in 2:nsteps) {
      
      ## Setup ##
      {
      # Reset flows
      si.flow[ts] <- ir.flow[ts] <- 0
      if (groups == 2) si.flow.g2[ts] <- ir.flow.g2[ts] <- 0
      
      ## Act rate per person per ts to edges per ts
      if (groups == 1) {
        acts <- round(act.rate * num[ts-1] / 2)
      }
      if (groups == 2) {
        if (balance == 'g1') {
          acts <- round(act.rate * (num[ts-1]+num.g2[ts-1]) / 2)
        }
        if (balance == 'g2') {
          acts <- round(act.rate.g2 * (num[ts-1]+num.g2[ts-1]) / 2)
        }
      }
      
      # NULL parameters for 1 group models
      if (groups == 1) {
        ids.g2 <- status.g2 <- duration.g2 <- trans.rate.g2 <-
          rec.rate.g2 <- ds.rate.g2 <- di.rate.g2 <- dr.rate.g2 <- 
          b.rate.g2 <- num.g2 <- NULL
      }
      }
      
      ### MODULES ###
      
      ## Infection Module
      {
      if (groups == 1) {
        p1 <- ssample(ids[status != 9], acts, replace=TRUE) 
        p2 <- ssample(ids[status != 9], acts, replace=TRUE)
        df <- data.frame(p1, p2)
        while (any(df$p1 == df$p2)) {
          df$p2 <- ifelse(df$p1 == df$p2, ssample(ids[status != 9], 1), df$p2)
        }
        df$p1.stat <- status[df$p1]
        df$p2.stat <- status[df$p2]
        df$serodis <- (df$p1.stat == 0 & df$p2.stat == 1) | 
          (df$p1.stat == 1 & df$p2.stat == 0)
        df <- df[df$serodis,]
        df$transmit <- runif(nrow(df)) <= trans.rate
        df <- df[df$transmit,]
        new.inf <- unique(ifelse(df$p1.stat == 0, df$p1, df$p2))
        status[new.inf] <- 1
        duration[new.inf] <- 0
        temp.si.flow <- length(new.inf)
      }
      
      if (groups == 2) { 
        p1 <- ssample(ids[status != 9], acts, replace=TRUE)
        p2 <- ssample(ids.g2[status.g2 != 9], acts, replace=TRUE)
        df <- data.frame(p1, p2)
        df$p1.stat <- status[df$p1]
        df$p2.stat <- status.g2[df$p2]
        df$serodis <- (df$p1.stat == 0 & df$p2.stat == 1) | 
          (df$p1.stat == 1 & df$p2.stat == 0)
        df <- df[df$serodis,]
        df$trans.rate <- ifelse(df$p1.stat == 0, trans.rate, trans.rate.g2)
        df$transmit <- runif(nrow(df)) <= df$trans.rate
        df <- df[df$transmit,]
        new.inf.g1 <- unique(df$p1[df$p1.stat == 0])
        status[new.inf.g1] <- 1
        duration[new.inf.g1] <- 0
        temp.si.flow <- length(new.inf.g1)
        new.inf.g2 <- unique(df$p2[df$p2.stat == 0])
        status.g2[new.inf.g2] <- 1
        duration.g2[new.inf.g2] <- 0
        temp.si.flow.g2 <- length(new.inf.g2)
      }
      }
      
      ## Recovery Module
      # Stochastic
      if (type != 'SI') {
        if (rec.rand == TRUE) {
          if (type == 'SIS') rec.flow <- 0
          if (type == 'SIR') rec.flow <- 2
          rec.elig <- sum(status == 1)
          rec.count <- rec.rate * rec.elig
          recs <- rpois(1, rec.count)
          if (recs > rec.elig) recs <- rec.elig
          temp.ir.flow <- recs
          if (rec.elig > 0) {
            status[ssample(ids[status == 1], recs)] <- rec.flow
          } 
          if (groups == 2) {
            rec.elig.g2 <- sum(status.g2 == 1)
            rec.count.g2 <- rec.rate.g2 * rec.elig.g2
            recs.g2 <- rpois(1, rec.count.g2)
            if (recs.g2 > rec.elig.g2) recs.g2 <- rec.elig.g2
            temp.ir.flow.g2 <- recs.g2
            if (rec.elig.g2 > 0) {
              status.g2[ssample(ids.g2[status.g2 == 1], recs.g2)] <- rec.flow
            } 
          }
        }
        # Deterministic
        if (rec.rand == FALSE) {
          dur.inf <- round(1/rec.rate)
          if (type == 'SIS') status[duration == dur.inf] <- 0
          if (type == 'SIR') status[duration == dur.inf] <- 2
          temp.ir.flow <- sum(duration == dur.inf, na.rm=T)
          duration[duration == dur.inf] <- NA
          duration[!is.na(duration)] <- duration[!is.na(duration)] + 1
          if (groups == 2) {
            dur.inf.g2 <- round(1/rec.rate.g2)
            temp.ir.flow.g2 <- sum(duration.g2 == dur.inf.g2, na.rm=T)
            duration.g2[duration.g2 == dur.inf.g2] <- NA
            duration.g2[!is.na(duration.g2)] <- duration.g2[!is.na(duration.g2)] + 1
          }
        }
      } else {
        temp.ir.flow <- temp.ir.flow.g2 <- 0
      }
      
      
      ## Mortality Module
      {
        ds.deaths <- di.deaths <- dr.deaths <- 0
        ds.deaths.g2 <- di.deaths.g2 <- dr.deaths.g2 <- 0
        if (ds.rate > 0) {
          ds.elig <- sum(status == 0)
          if (ds.elig > 0) {
            ds.dcnt <- ds.rate * ds.elig
            if (d.rand == TRUE) ds.deaths <- rpois(1, ds.dcnt)
            if (d.rand == FALSE) ds.deaths <- round(ds.dcnt)
            if (ds.deaths > ds.elig) ds.deaths <- ds.elig
            status[ssample(ids[status == 0], ds.deaths)] <- 9
          }
        }
        if (di.rate > 0) {
          di.elig <- sum(status == 1)
          if (di.elig > 0) {
            di.dcnt <- di.rate * di.elig
            if (d.rand == TRUE) di.deaths <- rpois(1, di.dcnt)
            if (d.rand == FALSE) di.deaths <- round(di.dcnt)
            if (di.deaths > di.elig) di.deaths <-  di.elig
            status[ssample(ids[status == 1], di.deaths)] <- 9
          }
        }
        if (dr.rate > 0) {
          dr.elig <- sum(status == 2)
          if (dr.elig > 0) {
            dr.dcnt <- dr.rate * dr.elig
            if (d.rand == TRUE) dr.deaths <- rpois(1, dr.dcnt)
            if (d.rand == FALSE) dr.deaths <- round(dr.dcnt)
            if (dr.deaths > dr.elig) dr.deaths <- dr.elig
            status[ssample(ids[status == 2], dr.deaths)] <- 9
          }
        }
        if (groups == 2) {
          if (ds.rate.g2 > 0) {
            ds.elig.g2 <- sum(status.g2 == 0)
            if (ds.elig.g2 > 0) {
              ds.dcnt.g2 <- ds.rate.g2 * ds.elig.g2
              if (d.rand == TRUE) ds.deaths.g2 <- rpois(1, ds.dcnt.g2)
              if (d.rand == FALSE) ds.deaths.g2 <- round(ds.dcnt.g2)
              if (ds.deaths.g2 > ds.elig.g2) ds.deaths.g2 <- ds.elig.g2
              status.g2[ssample(ids.g2[status.g2 == 0], ds.deaths.g2)] <- 9
            }
          }
          if (di.rate.g2 > 0) {
            di.elig.g2 <- sum(status.g2 == 1)
            if (di.elig.g2 > 0) {
              di.dcnt.g2 <- di.rate.g2 * di.elig.g2
              if (d.rand == TRUE) di.deaths.g2 <- rpois(1, di.dcnt.g2)
              if (d.rand == FALSE) di.deaths.g2 <- round(di.dcnt.g2)
              if (di.deaths.g2 > di.elig.g2) di.deaths.g2 <- di.elig.g2
              status.g2[ssample(ids.g2[status.g2 == 1], di.deaths.g2)] <- 9
            }
          }
          if (dr.rate.g2 > 0) {
            dr.elig.g2 <- sum(status.g2 == 2)
            if (dr.elig.g2 > 0) {
              dr.dcnt.g2 <- dr.rate.g2 * dr.elig.g2
              if (d.rand == TRUE) dr.deaths.g2 <- rpois(1, dr.dcnt.g2)
              if (d.rand == FALSE) dr.deaths.g2 <- round(dr.dcnt.g2)
              if (dr.deaths.g2 > dr.elig.g2) dr.deaths.g2 <- dr.elig.g2
              status.g2[ssample(ids.g2[status.g2 == 2], dr.deaths.g2)] <- 9
            }
          }
        }
      }
      
      
      ## Birth Module
      {
        if (groups == 1) {
          birth.count <- b.rate * num[ts-1]
        }
        if (groups == 2) {
          if (is.na(b.rate.g2)) {
            birth.count <- b.rate * num[ts-1] 
            birth.count.g2 <- b.rate * num[ts-1] 
          } else {
            birth.count <- b.rate * num[ts-1]
            birth.count.g2 <- b.rate.g2 * num.g2[ts-1]
          }
        }
        new.births <- new.births.g2 <- 0
        if (b.rate > 0) {
          if (b.rand == TRUE) new.births <- rpois(1, birth.count)
          if (b.rand == FALSE) new.births <- round(birth.count)
          if (new.births > 0) {
            ids <- c(ids, (max(ids)+1):(max(ids)+new.births))
          }
          status <- c(status, rep(0, new.births))
          duration <- c(duration, rep(NA, new.births))
        }
        if (groups == 2 && (b.rate.g2 > 0 | is.na(b.rate.g2))) {
          if (b.rand == TRUE) new.births.g2 <- rpois(1, birth.count.g2)
          if (b.rand == FALSE) new.births.g2 <- round(birth.count.g2)
          if (new.births.g2 > 0) {
            ids.g2 <- c(ids.g2, (max(ids.g2)+1):(max(ids.g2)+new.births.g2))
          }
          status.g2 <- c(status.g2, rep(0, new.births.g2))
          duration.g2 <- c(duration.g2, rep(NA, new.births.g2))
        }    
      }
      
      
      ### BOOKKEEPING ###
      {
        s.num[ts] <- sum(status == 0)
        i.num[ts] <- sum(status == 1)
        r.num[ts] <- sum(status == 2)
        num[ts] <- s.num[ts] + i.num[ts] + r.num[ts]
        si.flow[ts] <- temp.si.flow
        ir.flow[ts] <- temp.ir.flow
        b.flow[ts] <- new.births
        ds.flow[ts] <- ds.deaths
        di.flow[ts] <- di.deaths
        dr.flow[ts] <- dr.deaths
        if (groups == 2) {
          s.num.g2[ts] <- sum(status.g2 == 0)
          i.num.g2[ts] <- sum(status.g2 == 1)
          r.num.g2[ts] <- sum(status.g2 == 2)
          num.g2[ts] <- s.num.g2[ts] + i.num.g2[ts] + r.num.g2[ts]
          si.flow.g2[ts] <- temp.si.flow.g2
          ir.flow.g2[ts] <- temp.ir.flow.g2
          b.flow.g2[ts] <- new.births.g2
          ds.flow.g2[ts] <- ds.deaths.g2
          di.flow.g2[ts] <- di.deaths.g2
          dr.flow.g2[ts] <- dr.deaths.g2
        }
      }
      
    } # Timestep loop end #
    
    # Save simulation output
    if (s == 1) {
      out <- list(type = type,
                  groups = groups,
                  params = all.pout,
                  time = 1:nsteps,
                  s.num = data.frame(s.num),
                  i.num = data.frame(i.num),
                  si.flow = data.frame(si.flow))
      if (type == "SIR") {
        out$r.num <- data.frame(r.num)
        out$ir.flow <- data.frame(ir.flow)
      }
      if (type == "SIS") {
        out$is.flow <- data.frame(ir.flow)
      }
      if (vital == TRUE) {
        out$b.flow <- data.frame(b.flow)
        out$ds.flow <- data.frame(ds.flow)
        out$di.flow <- data.frame(di.flow)
        if (type == "SIR") {
          out$dr.flow <- data.frame(dr.flow)
        }
      }
      if (groups == 2) {
        out$s.num.g2 = data.frame(s.num.g2)
        out$i.num.g2 = data.frame(i.num.g2)
        out$si.flow.g2 = data.frame(si.flow.g2)
        if (type == "SIR") {
          out$r.num.g2 <- data.frame(r.num.g2)
          out$ir.flow.g2 <- data.frame(ir.flow.g2)
        }
        if (type == "SIS") {
          out$is.flow.g2 <- data.frame(ir.flow.g2)
        }
        if (vital == TRUE) {
          out$b.flow.g2 <- data.frame(b.flow.g2)
          out$ds.flow.g2 <- data.frame(ds.flow.g2)
          out$di.flow.g2 <- data.frame(di.flow.g2)
          if (type == "SIR") {
            out$dr.flow.g2 <- data.frame(dr.flow.g2)
          }
        }
      }
    }
    if (s > 1) {
      out$s.num[,s] <- s.num
      out$i.num[,s] <- i.num
      out$si.flow[,s] <- si.flow
      if (type == "SIR") {
        out$r.num[,s] <- r.num
        out$ir.flow[,s] <- ir.flow
      }
      if (type == "SIS") {
        out$is.flow[,s] <- ir.flow
      }
      if (vital == TRUE) {
        out$b.flow[,s] <- b.flow
        out$ds.flow[,s] <- ds.flow
        out$di.flow[,s] <- di.flow
        if (type == "SIR") {
          out$dr.flow[,s] <- dr.flow
        }
      }
      if (groups == 2) {
        out$s.num.g2[,s] <- s.num.g2
        out$i.num.g2[,s] <- i.num.g2
        out$si.flow.g2[,s] <- si.flow.g2
        if (type == "SIR") {
          out$r.num.g2[,s] <- r.num.g2
          out$ir.flow.g2[,s] <- ir.flow.g2
        }
        if (type == "SIS") {
          out$is.flow.g2[,s] <- ir.flow.g2
        }
        if (vital == TRUE) {
          out$b.flow.g2[,s] <- b.flow.g2
          out$ds.flow.g2[,s] <- ds.flow.g2
          out$di.flow.g2[,s] <- di.flow.g2
          if (type == "SIR") {
            out$dr.flow.g2[,s] <- dr.flow.g2
          }
        }
      }
    }
    
    # Track progress
    if (verbose == TRUE) {
      cat(paste('SIM=', s, '/', nsims, '\n', sep=''))
    }
    
  } # Simulation loop end
  
  
  # Set column names for varying list elements
  for (i in as.vector(which(lapply(out, class) == 'data.frame'))) {
    colnames(out[[i]]) <- paste('sim', 1:nsims, sep='')
  }
  # If only 1 sim, then change output from df to vectors
  if (nsims == 1) {
    for (i in as.vector(which(lapply(out, class) == 'data.frame'))) {
      out[[i]] <- out[[i]][,1]
    }
  }
  
  out$vital <- vital
  out$nsims <- nsims
  out$call <- match.call()
  
  class(out) <- 'epiICM'
  return(out)
}

