#' 
#' @title Deterministic Compartmental SI-SIS-SIR Models 
#'
#' @description This function solves a deterministic model of an infectious
#' disease with random mixing in the population.
#'
#' @param type model type (choice of 'SI', 'SIR', or 'SIS')
#' @param s.num number of initial susceptible in population
#' @param i.num number of initial infected in population
#' @param r.num number of initial recovered in population
#' @param beta probability of infection per contact
#' @param cont average number of contacts per person per unit time
#' @param nu average rate of recovery (1/duration of disease)
#' @param b birth rate
#' @param ms mortality rate for susceptibles
#' @param mi mortality rate for infecteds
#' @param mr mortality rate for recovered
#' @param dt either scalar or vector of timesteps for model
#' @param new.mod if not \code{NULL}, insert a new model definition
#' @param print.mod extract and print the current model from the function
#' @param verbose print progress
#' @param ... additional arguments to pass to model
#' 
#' @author Samuel M. Jenness <sjenness@@uw.edu>
#' @keywords model
#' @seealso \code{\link{plotDCM}}
#' @export
#' 
#' @examples
#' # SI model
#' out <- epiDCM(type='SI', s.num=500, i.num=1, beta=0.2, cont=0.25, dt=500)
#' par(mar=c(3.5,3,1,1), mgp=c(2,1,0))
#' plotDCM(out, c('s.num', 'i.num'), leg='full')
#' 
#' # SIR model
#' out <- epiDCM(type='SIR', s.num=500, i.num=1, r.num=0, 
#'               beta=0.2, cont=0.25, nu=1/50, dt=500)
#' plotDCM(out, c('s.num', 'i.num', 'r.num'), leg='full')
#' 
#' # SIS model
#' out <- epiDCM(type='SIS', s.num=500, i.num=1, r.num=0, 
#'               beta=0.2, cont=0.25, nu=1/50, dt=500)
#' plotDCM(out, c('s.num', 'i.num'), leg='full')
#' 
#' # SIR model with vital dynamics
#' out <- epiDCM(type='SIR', s.num=1000, i.num=1, r.num=0, 
#'               beta=0.2, cont=5, nu=1/3,
#'               b=1/90, ms=1/100, mi=1/35, mr=1/100, dt=500)
#' plotDCM(out, c('s.num', 'i.num', 'r.num'), leg='full')
#'
#' # SIS model with sensitivity parameter for contact rate
#' out <- epiDCM(type='SIS', s.num=500, i.num=1, r.num=0, 
#'               beta=0.2, cont=seq(0.05,0.5,0.05), 
#'               nu=1/50, dt=500, verbose=TRUE)
#' plotDCM(out, 's.num', leg='lim', alpha=0.6)
#' 
epiDCM <- function(type, s.num, i.num, r.num=0, beta, cont, nu=0,
                    b=0, ms=0, mi=0, mr=0, dt, new.mod=NULL, print.mod=FALSE, 
                    verbose=FALSE, ...){
  
  # External libraries
  require(deSolve)
  
  # Model type check
  if (!(type %in% c('SI', 'SIS', 'SIR'))) stop('Available model types are SI, SIS, or SIR')
  
  # Time steps
  if (length(dt)==1) dt <- 1:dt
  
  # Use this model or sub in a new one
  if (is.null(new.mod)) {
    if (type %in% c('SI', 'SIR')) {
      model <- function(t, t0, parms) {
        with(as.list(c(t0, parms)), {
          
          # derived totals
          num <- s.num + i.num + r.num
          
          # varying parameters
          lambda <- beta*cont*i.num/num
          
          # main ODEs
      		dS <- -lambda*s.num + b*num - ms*s.num
      		dI <- lambda*s.num- nu*i.num - mi*i.num
      		dR <- nu*i.num - mr*r.num
      		 
      	list(c(dS, dI, dR), num=num)
      	})
      }
    }
    if (type == 'SIS') {
      model <- function(t, t0, parms) {
        with(as.list(c(t0, parms)), {
          
          # derived totals
          num <- s.num + i.num + r.num
          
          # varying parameters
          lambda <- beta*cont*i.num/num
          
          # main ODEs
          dS <- -lambda*s.num + nu*i.num + b*num - ms*s.num 
          dI <- lambda*s.num - nu*i.num - mi*i.num
          dR <- 0
          
          list(c(dS, dI, dR), num=num)
        })
      }
    }
  } else {
    model <- new.mod
  }
  
  # Print out model or run it!
  if (print.mod) {
    print(model)
  } else {
  
    # Initial conditions
    t0 <- c(s.num = s.num,
            i.num = i.num,
            r.num = r.num)

    # Collect parameters and determine if any sensitivity parms
    all.p <- as.list(match.call(expand.dots = TRUE)[-1]) 
      if (missing(nu)) all.p$nu <- 0
      all.p[c('type', 'dt', 'new.mod', 'print.mod', 'verbose')] <- NULL
      all.p <- sapply(all.p, function(x) eval(x))
    nruns <- max(sapply(all.p, length))
    
    # Set up parameters and out matrix
    if (nruns > 1) {
      longv <- which(sapply(all.p, length) == max(sapply(all.p, length)))
      longvn <- names(longv)
      lim.p <- all.p[names(all.p) != longvn] 
    }
    
    # Main model runs
    outdf <- list()
    for (i in 1:nruns){
        
        # Sub in sensitivity parameter
        if (nruns > 1) {
          all.p <- as.list(c(lim.p, var=eval(as.symbol(longvn))[i]))   
          names(all.p)[which(names(all.p)=='var')] <- longvn
        }
        
        # Solve ODE
        df <- data.frame(ode(y=t0, times=dt, func=model, parms=all.p, method='rk4'))
        
        # Calculate model statistics
        s.prev <- df$s.num / df$num
        i.prev <- df$i.num / df$num
        r.prev <- df$r.num / df$num
        incid <- (all.p[['beta']] * (all.p[['cont']]/df$num)) * df$s.num * df$i.num
        R0 <- all.p[['beta']] * all.p[['cont']] * (1/all.p[['nu']])
        Rn <- all.p[['beta']] * all.p[['cont']] * (1/all.p[['nu']]) * (df$s.num/df$num)
        
        # Collect in list
        if (nruns > 1){
          outdf[[i]] <- data.frame(df, s.prev, i.prev, r.prev, incid, R0, Rn)
          outdf[[i]][longvn] <- eval(as.symbol(longvn))[i]
        } else {
          outdf[[i]] <- data.frame(df, s.prev, i.prev, r.prev, incid, R0, Rn)
        }
        
        # Track progress
        if (verbose) { cat('Run=', i, '/', nruns, '\n', sep='') }
    }
    
    # Write out all data to list
    all <- do.call(cbind, outdf)
    out <- list(type=type,
                time=all[,1],
                s.num=all[,names(all)=='s.num'],
                i.num=all[,names(all)=='i.num'],
                r.num=all[,names(all)=='r.num'],
                s.prev=all[,names(all)=='s.prev'],
                i.prev=all[,names(all)=='i.prev'],
                r.prev=all[,names(all)=='r.prev'],
                i.incid=all[,names(all)=='incid'],
                R0=all[1,names(all)=='R0'],
                Rn=all[,names(all)=='Rn'])
    if (nruns > 1) {
      out$longvn=all[1,names(all)==names(all)[ncol(all)]]
      names(out)[length(out)] <- names(all)[ncol(all)]
      
      vartype <- vector()
      for (i in 1:length(out)) vartype[i] <- class(out[[i]])
      dfvars <- which(vartype=='data.frame')
      for (i in dfvars) {
        colnames(out[[i]]) <- paste(names(out)[length(out)], 
                                    round(out[[length(out)]], 3), sep='')
      }
    }
    
    invisible(out)
    }
}

