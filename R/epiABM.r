#' 
#' @title Stochastic Agent-Based SI-SIS-SIR Models 
#'
#' @description This function simulates an agent-based model of an infectious
#' disease with random mixing in the population.
#'
#' @param type model type (choice of 'SI', 'SIR', or 'SIS')
#' @param s.num number of initial susceptible in population
#' @param i.num number of initial infected in population
#' @param r.num number of initial recovered in poulation
#' @param beta probability of infection per contact
#' @param cont average number of contacts per person per unit time
#' @param nu average rate of recovery (1/duration of disease)
#' @param nsteps number of time steps to simulate over
#' @param nsims number of simulations to run
#' @param verbose print progress
#' 
#' @details This function provides flexibility in model choice of a susceptible-infected (SI),
#' susceptible-infected-recovered (SIR), or a susceptible-infected-susceptible (SIS) epidemic
#' model. For an SI model, use the defaults of 0 for \code{r.num}, \code{nu}, and \code{mu}. For
#' an SIR model, specify the recovery rate with the \code{nu} parameter, but leave \code{mu} to
#' its 0 default. For an SIS model, specify the re-susceptiblity with the \code{mu} parameter,
#' but leave the \code{nu} to its 0 default.
#' 
#' The model has two sources of stochasticity. First, the contact structure is random, with 
#' new partnerships selected at random from among the actors. The contact rate argument 
#' \code{cont} is a rate per person per timestep, and is transformed to the number of new 
#' partnerships per timestep by cont*(s.num+i.num)/2. Second, the probability of transmission 
#' within each partnership is based on a standard uniform random deviate falling below \code{beta}.
#' 
#' @author Samuel M. Jenness <sjenness@@uw.edu> and Steven M. Goodreau <goodreau@@uw.edu>
#' @keywords model
#' @seealso \code{\link{plotABM}}
#' @export
#' 
#' @examples
#' # SI model
#' out <- epiABM(type='SI', s.num=500, i.num=1, beta=0.2, cont=0.25, nsteps=250, nsims=10)
#' par(mar=c(3,3,1,1), mgp=c(2,1,0))
#' plotABM(out, 'i.num', medline=FALSE)
#' 
#' # SIR model
#' out <- epiABM(type='SIR', s.num=500, i.num=1, beta=0.2, cont=0.25, nu=1/50, nsteps=250, nsims=10)
#' plotABM(out, 'i.num', medline=FALSE)
#' 
#' # SIS model
#' out <- epiABM(type='SIS', s.num=500, i.num=1, beta=0.2, cont=0.25, nu=1/50, nsteps=250, nsims=10)
#' plotABM(out, 'i.num', medline=FALSE)
#' 
epiABM <- function(type, s.num, i.num, r.num=0, beta, cont, nu=0,
                   nsteps, nsims, verbose=TRUE) {
    
    s.num.t0 <- s.num
    i.num.t0 <- i.num
    r.num.t0 <- r.num
    num.t0 <- s.num + i.num + r.num
    ids <- 1:num.t0
    
    contacts <- round(cont*num.t0/2)
    
    # Model type check
    if (!(type %in% c('SI', 'SIS', 'SIR'))) stop('Available model types are SI, SIS, or SIR')
        
    dur.inf <- 1/nu
    
    out.s.num <- out.i.num <- out.r.num <- data.frame(rep(NA, nsteps))
    out.i.incid <- out.r.incid <- data.frame(rep(NA, nsteps))
    out.eq.time <- rep(NA, nsims)
    
    for (s in 1:nsims) {
      
      status <- rep(0, num.t0)
      status[1:i.num.t0] <- 1
      duration <- rep(NA, num.t0)
      duration[1:i.num.t0] <- 0
      
      s.num <- s.num.t0
      i.num <- i.num.t0
      r.num <- r.num.t0
      i.incid <- 0
      r.incid <- 0
      eq.time <- NA
      
      for (ts in 2:nsteps) {
        
        i.incid[ts] <- 0
        
        if (!is.na(eq.time))  { 
          s.num[ts] <- s.num[ts-1]
          i.num[ts] <- i.num[ts-1]
          r.num[ts] <- r.num[ts-1]
          if (type == 'SIR') {
            r.incid[ts] <- sum(duration == dur.inf, na.rm=T)
          }
          next
        } else {
          for (cnt in 1:contacts) {
            partners <- sample(ids, 2)
            if (status[partners[1]] + status[partners[2]] == 1) {
              if (runif(1) <= beta) {
                i.incid[ts] <- i.incid[ts]+1    
                if (status[partners[1]] == 1) {
                  status[partners[2]] <- 1
                  duration[partners[2]] <- 0
                } else {
                  status[partners[1]] <- 1
                  duration[partners[1]] <- 0
                }
              }
            }  
          }
        
          if (type == 'SIS') {
            status[duration == dur.inf] <- 0
          }
          if (type == 'SIR') {
            status[duration == dur.inf] <- 2
            r.incid[ts] <- sum(duration == dur.inf, na.rm=T)
          }
          
          duration[duration == dur.inf] <- NA
          duration[!is.na(duration)] <- duration[!is.na(duration)]+1
          
          s.num[ts] <- sum(status == 0)
          i.num[ts] <- sum(status == 1)
          r.num[ts] <- sum(status == 2)
          
          if (i.num[ts] == 0) eq.time <- ts
          
        }
      }
      
      out.s.num[,s] <- s.num
      out.i.num[,s] <- i.num
      out.r.num[,s] <- r.num
      out.i.incid[,s] <- i.incid
      out.r.incid[,s] <- r.incid    
      out.eq.time[s] <- eq.time
      
      if (verbose){
        cat(paste('SIM=', s, '/', nsims, '\n', sep=''))
      }
    }
    
    out <- list(type=type,
                time=1:nsteps, 
                s.num=out.s.num, 
                i.num=out.i.num, 
                r.num=out.r.num, 
                i.incid=out.i.incid, 
                r.incid=out.r.incid, 
                eq.time=out.eq.time)
    
    vartype <- vector()
      for (i in 1:length(out)) vartype[i] <- class(out[[i]])
      dfvars <- which(vartype=='data.frame')
      for (i in dfvars) colnames(out[[i]]) <- paste('sim', 1:nsims, sep='')
    
    
  return(out)
}
