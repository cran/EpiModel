#' 
#' @title Deterministic Compartmental Epidemic Models 
#'
#' @description This function solves deterministic compartmental epidemic models 
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
#' @param dt fractional time unit for ODE solutions, with the default of 1. Model 
#'   solutions for fractional timesteps may be obtained by setting this to a 
#'   number between 0 and 1. 
#' @param odemethod ODE integration method,  with the default for the Runge-Kutta 
#'   4 method (see \code{\link{ode}} for other options)
#' @param new.mod if not \code{NULL}, this parameter requires an object containing 
#'   a function specifying a new mathematical model form to be solved over. This 
#'   functionality allows the user to specify a different system of differential 
#'   equations than those currently embedded in \code{epiDCM}. 
#' @param print.mod if \code{TRUE}, running \code{epiDCM} will print the model form
#'   given the specified \code{type} and \code{groups}. This is used to see the 
#'   mathematical structure of the current model form, as well as obtain a base 
#'   functional form to modify and input using the \code{new.mod} parameter. 
#' @param verbose if \code{TRUE}, print modelling progress to the console
#' @param ... additional parameters to pass to model, currently not used 
#' 
#' @details 
#' The \code{epiDCM} function is uses the ordinary differential equation solver in
#' the \code{deSolve} package to model disease as a deterministic compartmental 
#' system. The current implementation of this function is limited to one- and two-
#' group models with disease types for Susceptible-Infected (SI), Susceptible-
#' Infected-Recovered (SIR), and Susceptible-Infected-Susceptible (SIS). The
#' one-group model, specified with the parameter \code{groups=1}, assumes purely
#' random mixing in the population. The two-group model, \code{groups=2}, models
#' purely heterogenous mixing in the population, wherein group 1 members may only
#' have contact with group 2 members (e.g., a heterosexual-only disease transmission).
#' This will be further expanded in future releases of \code{EpiModel}.
#'   
#' One important caveat is the current parameter structure of specifing initial state
#' sizes and individual parameters as primary \code{epiDCM} parameters will be 
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
#' @section Sensitivity Analyses:
#' The \code{epiDCM} function has been designed to facilitate DCM sensitivity 
#' analyses, wherein one may run a series of models while varying one or more
#' of the model parameters. This is possible, as shown in the example below,
#' by specifying any parameter as a vector of length greater than one. See the 
#' Section 2.3 of the \href{../doc/Tutorial.pdf}{EpiModel Tutorial vignette.}
#' for further explanation and example.
#' 
#' @return
#' This function returns a list object of class \code{epiDCM} with the following 
#' elements:
#' \itemize{
#'  \item \strong{type:} disease type as specified in model parameter.
#'  \item \strong{groups:} groups as specified in model parameter.
#'  \item \strong{params:} list of model parameters as specified in model parameter.
#'  \item \strong{time:} a vector of time steps over which the model was solved.
#'  \item \strong{vital:} logical, whether vital dynamics were specified in the 
#'  parameterization.
#'  \item \strong{nruns:} number of independent model runs.
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
#' @seealso 
#' Analysis, plotting, and extraction of model data are available through the 
#'   \code{\link{summary.epiDCM}}, \code{\link{plot.epiDCM}}, and 
#'   \code{\link{as.data.frame.epiDCM}} functions, respectively. For further 
#'   discussion of expanding DCMs past the current parameterizations, see the
#'   HTML vignette \href{../doc/ExtendepiDCM.html}{Extending DCMs Past EpiModel}.
#' @export
#' 
#' @examples
#' ## SI Model (One-Group)
#' mod <- epiDCM(type="SI", s.num = 500, i.num = 1, 
#'               trans.rate = 0.2, act.rate = 0.25, nsteps = 500)
#' plot(mod)
#' 
#' ## SIR Model with Vital Dynamics (One-Group)
#' mod <- epiDCM(type = "SIR", s.num=1000, i.num = 1, r.num = 0, 
#'               trans.rate = 0.2, act.rate = 5, rec.rate = 1/3,
#'               b.rate = 1/90, ds.rate = 1/100, di.rate = 1/35, 
#'               dr.rate = 1/100, nsteps = 500)
#' plot(mod)
#'
#' ## SIS Model with act.rate Sensitivity Parameter
#' mod <- epiDCM(type="SIS", s.num = 500, i.num = 1, r.num = 0, 
#'               trans.rate = 0.2, act.rate = seq(0.1, 0.5, 0.1), 
#'               rec.rate = 1/50, nsteps = 500, verbose = TRUE)
#' plot(mod)
#' 
#' ## SI Model with Vital Dynamics (Two Groups)
#' mod <- epiDCM(type = "SI", groups = 2, 
#'               s.num = 500, i.num = 1, s.num.g2 = 500, i.num.g2 = 0, 
#'               trans.rate = 0.4,  trans.rate.g2 = 0.1, 
#'               act.rate = 0.25, balance = "g1",
#'               b.rate = 1/100, b.rate.g2 = NA,
#'               ds.rate = 1/100, ds.rate.g2 = 1/100,
#'               di.rate = 1/50, di.rate.g2 = 1/50,
#'               nsteps = 500)
#' plot(mod)
#' 
epiDCM <- function(type, 
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
                   r.num.g2=0, 
                   trans.rate.g2, 
                   act.rate.g2, 
                   rec.rate.g2 = 0, 
                   b.rate.g2 = 0, 
                   ds.rate.g2 = 0, 
                   di.rate.g2 = 0, 
                   dr.rate.g2 = 0, 
                   balance, 
                   nsteps, 
                   dt, 
                   odemethod,
                   new.mod = NULL, 
                   print.mod = FALSE, 
                   verbose = FALSE, 
                   ...){
  
  # Model type check
  if (missing(type) | !(type %in% c('SI','SIS','SIR'))) {
    stop('Available model types are: SI SIS SIR')
  }

  # Group number check
  if (!(groups %in% 1:2)) {
    stop('Only groups=1 and groups=2 are allowed')
  }
    
  # Collect parameters and determine if any sensitivity parms
  all.p <- as.list(match.call(expand.dots = TRUE)[-1]) 
  if (missing(rec.rate)) all.p$rec.rate <- 0
  if (missing(rec.rate.g2) & groups == 2) all.p$rec.rate.g2 <- 0
  all.p <- all.p[-(grep('num', names(all.p)))]
  all.p <- all.p[-which(names(all.p) %in% c('type', 'groups', 
                                            'balance','dt', 'nsteps',
                                            'odemethod', 'new.mod', 
                                            'print.mod','verbose'))]
  all.p <- sapply(names(all.p), function(x) eval(parse(text=x)))
  #all.p <- sapply(all.p, function(x) eval(x))
  all.pout <- as.list(all.p)
  
  # Determine if any vital dynamics
  vd.params <- c('b.rate', 'ds.rate', 'di.rate', 'dr.rate',
                 'b.rate.g2', 'ds.rate.g2', 'di.rate.g2', 'dr.rate.g2')
  if (any(names(all.pout) %in% vd.params)) { 
    vital <- TRUE 
    if (groups == 1) {
      if (max(c(b.rate, ds.rate, di.rate, dr.rate)) == 0)
        vital <- FALSE
    }
    if (groups == 2) {
      if (max(c(b.rate, ds.rate, di.rate, dr.rate,
                b.rate.g2, ds.rate.g2, di.rate.g2, dr.rate.g2), na.rm=T)==0)
        vital <- FALSE
    }
  } else { 
    vital <- FALSE 
  }
  
  
  # Number of model runs
  nruns <- max(sapply(all.p, length))
  
  # Set up parameters and out matrix
  if (nruns > 1) {
    longv <- which(sapply(all.p, length) == max(sapply(all.p, length)))
    longvn <- names(longv)
    lim.p <- all.p[!(names(all.p) %in% longvn)] 
  }
  
  # Balance check
  if (groups == 2 && (missing(balance) || !(balance %in% c('g1', 'g2')))) {
    stop("Must specify balance='g1' or balance='g2' with 2-group models")
  }
   
  # SI parameter check
  if (type == 'SI' && (rec.rate > 0 | rec.rate.g2 > 0)) { 
    stop('Set rec.rate=0 since this is an SI model')
  }
  
  # dr.rate warning for SI/SIS models
  if (type %in% c('SI', 'SIS') && (dr.rate > 0 | dr.rate.g2 > 0)) {
    warning('dr.rate parameter specified, but not used in SI/SIS models')
  }
  
  # Time steps fraction for ODE
  if (missing(dt)) dt <- 1
  dt <- seq(1, nsteps, dt)
  
  # ODE method
  if (missing(odemethod)) odemethod <- 'rk4'
  
  # Use this model or sub in a new one
  if (is.null(new.mod)) {
    
    # 1. SIR/SI 1 group
    if (type %in% c('SI', 'SIR') && groups == 1) {
      model <- function(t, t0, parms) {
        with(as.list(c(t0, parms)), {
          
          # derived totals
          num <- s.num + i.num + r.num
          
          # varying parameters
          lambda <- trans.rate * act.rate * i.num / num
          
          # main ODEs
      		dS <- -lambda*s.num + b.rate*num - ds.rate*s.num
      		dI <- lambda*s.num - rec.rate*i.num - di.rate*i.num
      		dR <- rec.rate*i.num - dr.rate*r.num
      		 
      	list(c(dS, dI, dR), num=num, lambda=lambda)
      	})
      }
    }
    
    # 2. SIR/SI 2 group
    if (type %in% c('SI', 'SIR') && groups == 2) {
      model <- function(t, t0, parms) {
        with(as.list(c(t0, parms)), {
          
          # derived totals
          num.g1 <- s.num + i.num + r.num
          num.g2 <- s.num.g2 + i.num.g2 + r.num.g2
          num <- num.g1 + num.g2
          
          # act rate balancing
          if (balance == 'g1') {
            ct.g1 <- act.rate
            ct.g2 <- ct.g1 * num.g1 / num.g2
          }
          if (balance == 'g2') {
            ct.g2 <- act.rate.g2
            ct.g1 <- ct.g2 * num.g2 / num.g1
          }
          
          # group-specific foi
          lambda.g1 <- trans.rate * ct.g1 * i.num.g2 / num.g2
          lambda.g2 <- trans.rate.g2 * ct.g2 * i.num / num.g1
          
          # birth rates
          if (is.na(b.rate.g2)) {
            br.g1 <- b.rate*num.g1
            br.g2 <- b.rate*num.g1
          } else {
            br.g1 <- b.rate*num.g1
            br.g2 <- b.rate*num.g2
          }
          
          # main ODEs
          dSm1 <- -lambda.g1*s.num + br.g1 - ds.rate*s.num
          dIm1 <- lambda.g1*s.num - rec.rate*i.num - di.rate*i.num
          dRm1 <- rec.rate*i.num - dr.rate*r.num
          
          dSm2 <- -lambda.g2*s.num.g2 + br.g2 - ds.rate.g2*s.num.g2
          dIm2 <- lambda.g2*s.num.g2 - rec.rate.g2*i.num.g2 - di.rate.g2*i.num.g2
          dRm2 <- rec.rate.g2*i.num.g2 - dr.rate.g2*r.num.g2
          
          list(c(dSm1, dIm1, dRm1, dSm2, dIm2, dRm2), 
               num.g1=num.g1, num.g2=num.g2, num=num,
               lambda.g1=lambda.g1, lambda.g2=lambda.g2)
        })
      }
    }
    
    # 3. SIS 1 group
    if (type == 'SIS' & groups == 1) {
      model <- function(t, t0, parms) {
        with(as.list(c(t0, parms)), {
          
          # derived totals
          num <- s.num + i.num
          
          # varying parameters
          lambda <- trans.rate * act.rate * i.num / num
          
          # main ODEs
          dS <- -lambda*s.num + rec.rate*i.num + b.rate*num - ds.rate*s.num 
          dI <- lambda*s.num - rec.rate*i.num - di.rate*i.num
          dR <- 0
          
          list(c(dS, dI, dR), num=num, lambda=lambda)
        })
      }
    }
    
    # 4. SIS 2 group
    if (type == 'SIS' & groups == 2) {
      model <- function(t, t0, parms) {
        with(as.list(c(t0, parms)), {
          
          # derived totals
          num.g1 <- s.num + i.num
          num.g2 <- s.num.g2 + i.num.g2
          num <- num.g1 + num.g2
          
          # act rate balancing
          if (balance == 'g1') {
            ct.g1 <- act.rate
            ct.g2 <- ct.g1 * num.g1 / num.g2
          }
          if (balance == 'g2') {
            ct.g2 <- act.rate.g2
            ct.g1 <- ct.g2 * num.g2 / num.g1
          }
          
          # group-specific foi
          lambda.g1 <- trans.rate * ct.g1 * i.num.g2/num.g2
          lambda.g2 <- trans.rate.g2 * ct.g2 * i.num/num.g1
          
          # birth rates
          if (is.na(b.rate.g2)) {
            br.g1 <- b.rate*num.g1
            br.g2 <- b.rate*num.g1
          } else {
            br.g1 <- b.rate*num.g1
            br.g2 <- b.rate*num.g2
          }
          
          # main ODEs
          dSm1 <- -lambda.g1*s.num + rec.rate*i.num + br.g1 - ds.rate*s.num
          dIm1 <- lambda.g1*s.num - rec.rate*i.num - di.rate*i.num
          dRm1 <- 0
          
          dSm2 <- -lambda.g2*s.num.g2 + rec.rate.g2*i.num.g2 + br.g2 - ds.rate.g2*s.num.g2
          dIm2 <- lambda.g2*s.num.g2 - rec.rate.g2*i.num.g2 - di.rate.g2*i.num.g2
          dRm2 <- 0
          
          list(c(dSm1, dIm1, dRm1, dSm2, dIm2, dRm2), 
               num.g1=num.g1, num.g2=num.g2, num=num,
               lambda.g1=lambda.g1, lambda.g2=lambda.g2)
        })
      }
    }
  
  } else {
    # Sub in new model
    model <- new.mod
  }
  
  # Print out model or run it
  if (print.mod) {
    print(model)
  } else {
  
    # Initial conditions
    if (groups == 1) {
      t0 <- c(s.num = s.num,
              i.num = i.num,
              r.num = r.num)
    }
    if (groups == 2) {
      t0 <- c(s.num = s.num,
              i.num = i.num,
              r.num = r.num,
              s.num.g2 = s.num.g2,
              i.num.g2 = i.num.g2,
              r.num.g2 = r.num.g2)
    }
    
    # Model runs
    outdf <- list()
    for (i in 1:nruns){
        
        # Sub in sensitivity parameter
        if (nruns > 1) {
          sens.p <- sapply(longvn, function(x) eval(as.symbol(x))[i])
          all.p <- as.list(c(lim.p, sens.p))   
        }
        
        # Solve ODE
        df <- data.frame(ode(y=t0, times=dt, func=model, parms=all.p, method=odemethod))
        
          # Calculate flows
          if (groups == 1) {
            si.flow <- df$lambda * df$s.num
            if (type == 'SIR') ir.flow <- all.p[['rec.rate']] * df$i.num
            if (type == 'SIS') is.flow <- all.p[['rec.rate']] * df$i.num
          }
          if (groups == 2) {
              si.flow <- df$lambda.g1 * df$s.num
              si.flow.g2 <- df$lambda.g2 * df$s.num.g2
            if (type == 'SIR') {
              ir.flow <- all.p[['rec.rate']] * df$i.num
              ir.flow.g2 <- all.p[['rec.rate.g2']] * df$i.num.g2
            }
            if (type == 'SIS') {
              is.flow <- all.p[['rec.rate']] * df$i.num
              is.flow.g2 <- all.p[['rec.rate.g2']] * df$i.num.g2
            }
          }
        
          # Add vital dynamics flows as necessary
          if (groups == 1) {
            if (any(names(all.p) == 'b.rate')) b.flow <- all.p[['b.rate']] * df$num
          }
          if (any(names(all.p) == 'ds.rate')) ds.flow <- all.p[['ds.rate']] * df$s.num
          if (any(names(all.p) == 'di.rate')) di.flow <- all.p[['di.rate']] * df$i.num
          if (any(names(all.p) == 'dr.rate')) dr.flow <- all.p[['dr.rate']] * df$r.num
          if (groups == 2) {
            if (any(names(all.p) == 'b.rate')) b.flow <- all.p[['b.rate']] * df$num.g1
            if (any(names(all.p) == 'b.rate.g2') && !is.na(all.p[['b.rate.g2']]))
              b.flow.g2 <- all.p[['b.rate']] * df$num.g2
            if (any(names(all.p) == 'b.rate.g2') && is.na(all.p[['b.rate.g2']]))
              b.flow.g2 <- all.p[['b.rate']] * df$num.g1
            if (any(names(all.p) == 'ds.rate.g2')) ds.flow.g2 <- all.p[['ds.rate.g2']] * df$s.num.g2
            if (any(names(all.p) == 'di.rate.g2')) di.flow.g2 <- all.p[['di.rate.g2']] * df$i.num.g2
            if (any(names(all.p) == 'dr.rate.g2')) dr.flow.g2 <- all.p[['dr.rate.g2']] * df$r.num.g2
          }

        
        ## Collect in list
        if (groups == 1) { 
          if (type == 'SI') {
            outdf[[i]] <- data.frame(df, si.flow)
          }
          if (type == 'SIR') {
            outdf[[i]] <- data.frame(df, si.flow, ir.flow)
          }
          if (type == 'SIS') {
            outdf[[i]] <- data.frame(df, si.flow, is.flow)
          }
        }
        if (groups == 2) {
          if (type == 'SI') {
            outdf[[i]] <- data.frame(df, si.flow, si.flow.g2)
          }
          if (type == 'SIR') {
            outdf[[i]] <- data.frame(df, si.flow, si.flow.g2, ir.flow, ir.flow.g2)
          }
          if (type == 'SIS') {
            outdf[[i]] <- data.frame(df, si.flow, si.flow.g2, is.flow, is.flow.g2)
          }
        }
        # Save out vital dynamics if necessary
        if (any(names(all.p) == 'b.rate')) outdf[[i]]['b.flow'] <- b.flow
        if (any(names(all.p) == 'ds.rate')) outdf[[i]]['ds.flow'] <- ds.flow
        if (any(names(all.p) == 'di.rate')) outdf[[i]]['di.flow'] <- di.flow
        if (any(names(all.p) == 'dr.rate')) outdf[[i]]['dr.flow'] <- dr.flow
        if (any(names(all.p) == 'b.rate.g2')) outdf[[i]]['b.flow.g2'] <- b.flow.g2
        if (any(names(all.p) == 'ds.rate.g2')) outdf[[i]]['ds.flow.g2'] <- ds.flow.g2
        if (any(names(all.p) == 'di.rate.g2')) outdf[[i]]['di.flow.g2'] <- di.flow.g2
        if (any(names(all.p) == 'dr.rate.g2')) outdf[[i]]['dr.flow.g2'] <- dr.flow.g2
        # Save out sensitivity parameter
        if (nruns > 1) {
          outdf[[i]][longvn] <- eval(as.symbol(longvn))[i]
        }
    
        # Track progress
        if (verbose) { cat('Run=', i, '/', nruns, '\n', sep='') }
    } # end of run loop
    
    # Write out all data to list
    all <- do.call(cbind, outdf)
    if (groups == 1) {
      if (type == 'SI') {
        out <- list(type=type,
                    groups=groups,
                    params=all.pout,
                    time=all[, 1],
                    s.num=all[, names(all) == 's.num'],
                    i.num=all[, names(all) == 'i.num'],
                    si.flow=all[, names(all) == 'si.flow'])
      }
      if (type == 'SIR') {
        out <- list(type=type,
                    groups=groups,
                    params=all.pout,
                    time=all[, 1],
                    s.num=all[, names(all) == 's.num'],
                    i.num=all[, names(all) == 'i.num'],
                    r.num=all[, names(all) == 'r.num'],
                    si.flow=all[, names(all) == 'si.flow'],
                    ir.flow=all[, names(all) == 'ir.flow'])
      }
      if (type == 'SIS') {
        out <- list(type=type,
                    groups=groups,
                    params=all.pout,
                    time=all[, 1],
                    s.num=all[, names(all) == 's.num'],
                    i.num=all[, names(all) == 'i.num'],
                    si.flow=all[, names(all) == 'si.flow'],
                    is.flow=all[, names(all) == 'is.flow'])
      }
    }
    if (groups == 2) {
      if (type == 'SI') {
        out <- list(type=type,
                    groups=groups,
                    params=all.pout,
                    time=all[, 1],
                    s.num=all[, names(all) == 's.num'],
                    s.num.g2=all[, names(all) == 's.num.g2'],
                    i.num=all[, names(all) == 'i.num'],
                    i.num.g2=all[, names(all) == 'i.num.g2'],
                    si.flow=all[, names(all) == 'si.flow'],
                    si.flow.g2=all[, names(all) == 'si.flow.g2'])
      }
      if (type == 'SIR') {
        out <- list(type=type,
                    groups=groups,
                    params=all.pout,
                    time=all[, 1],
                    s.num=all[, names(all) == 's.num'],
                    s.num.g2=all[, names(all) == 's.num.g2'],
                    i.num=all[, names(all) == 'i.num'],
                    i.num.g2=all[, names(all) == 'i.num.g2'],
                    r.num=all[, names(all) == 'r.num'],
                    r.num.g2=all[, names(all) == 'r.num.g2'],
                    si.flow=all[, names(all) == 'si.flow'],
                    si.flow.g2=all[, names(all) == 'si.flow.g2'],
                    ir.flow=all[, names(all) == 'ir.flow'],
                    ir.flow.g2=all[, names(all) == 'ir.flow.g2'])
      }
      if (type == 'SIS') {
        out <- list(type=type,
                    groups=groups,
                    parms=all.pout,
                    time=all[, 1],
                    s.num=all[, names(all) == 's.num'],
                    s.num.g2=all[, names(all) == 's.num.g2'],
                    i.num=all[, names(all) == 'i.num'],
                    i.num.g2=all[, names(all) == 'i.num.g2'],
                    si.flow=all[, names(all) == 'si.flow'],
                    si.flow.g2=all[, names(all) == 'si.flow.g2'],
                    is.flow=all[, names(all) == 'is.flow'],
                    is.flow.g2=all[, names(all) == 'is.flow.g2'])
      }
    }
    if (any(b.rate > 0)) 
      out$b.flow = all[, names(all) == 'b.flow']
    if (is.na(b.rate.g2) || any(b.rate.g2 > 0)) 
      out$b.flow.g2 = all[, names(all) == 'b.flow.g2']
    if (any(ds.rate > 0)) 
      out$ds.flow = all[, names(all) == 'ds.flow']
    if (any(ds.rate.g2 > 0)) 
      out$ds.flow.g2 = all[, names(all) == 'ds.flow.g2']
    if (any(di.rate > 0)) 
      out$di.flow = all[, names(all) == 'di.flow']
    if (any(di.rate.g2 > 0)) 
      out$di.flow.g2 = all[, names(all) == 'di.flow.g2']
    if (any(dr.rate > 0)) 
      out$dr.flow = all[, names(all) == 'dr.flow']
    if (any(dr.rate.g2 > 0)) 
      out$dr.flow.g2 = all[, names(all) == 'dr.flow.g2']
    
    if (nruns > 1) {
      out$longvn=all[1, names(all) == names(all)[ncol(all)]]
      names(out)[length(out)] <- names(all)[ncol(all)]
      
      for (i in as.vector(which(lapply(out, class) == 'data.frame'))) {
        colnames(out[[i]]) <- paste(names(out)[length(out)], 
                                    round(out[[length(out)]], 3), sep='')
      }
    }
    
    out$vital <- vital
    out$nruns <- nruns
    out$call <- match.call()
    
    class(out) <- 'epiDCM'
    invisible(out)
    }
}