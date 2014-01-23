
#' @title Extract Model Run for epiDCM Object
#'
#' @description This function extracts a model run from an object of class
#'   \code{epiDCM} into a data frame using the generic \code{as.data.frame} 
#'   function.
#'
#' @param x an \code{EpiModel} object of class \code{\link{epiDCM}}.
#' @param run run number for model, for multiple-run sensitivity models 
#'   run in \code{\link{epiDCM}}.
#' @param row.names see \code{\link{as.data.frame.default}}.
#' @param optional see \code{\link{as.data.frame.default}}.
#' @param ...  see \code{\link{as.data.frame.default}}.
#' 
#' @details
#' Model output from an \code{epiDCM} run are available through a data frame
#' with this helper function. The output data frame will include columns for
#' time, the size of each compartment, the overall population size 
#' (the sum of compartment sizes), and the size of each flow.
#' 
#' @method as.data.frame epiDCM
#' @keywords extract
#' @export
#' 
#' @examples
#' # One group SIS model with varying act.rate
#' mod <- epiDCM(type="SIS", s.num=500, i.num=1, r.num=0, 
#'               trans.rate=0.2, act.rate=seq(0.05,0.5,0.05), 
#'               rec.rate=1/50, nsteps=500, verbose=TRUE)
#' head(as.data.frame(mod, run=3))
#' 
#' # Two group SIR model with vital dynamics
#' mod <- epiDCM(type = "SIR", groups = 2,
#'               s.num = 500, i.num = 1, s.num.g2 = 500, i.num.g2 = 1,
#'               trans.rate = 0.2, trans.rate.g2 = 0.1, act.rate = 3,
#'               balance = "g1", rec.rate = 1/50, rec.rate.g2 = 1/50,
#'               b.rate = 1/100, b.rate.g2 = NA, ds.rate = 1/100,
#'               ds.rate.g2 = 1/100, di.rate = 1/90, di.rate.g2 = 1/90,
#'               dr.rate = 1/100, dr.rate.g2 = 1/100, nsteps = 500)
#' head(as.data.frame(mod))
#'  
as.data.frame.epiDCM <- function(x, 
                                 row.names=NULL, 
                                 optional=FALSE, 
                                 run=1,
                                 ...) {
  
  new.x <- data.frame(time=x$time)
  nruns <- x$nruns
  groups <- x$groups
  
  out.vars <- c(grep("num", names(x)), grep("flow", names(x)))
  out.vars.n <- names(x)[out.vars]
  
  # Output for models with 1 run
  if (nruns == 1) {
    if (run > 1) stop("Specify run = 1")
    for (i in seq_along(out.vars)) new.x[,i+1] <- x[[ out.vars[i] ]]
  }
  
  # Output for models with multiple runs
  if (nruns > 1) {
    if (run > nruns) {
      stop(paste("Specify run between 1 and", nruns))
    }
    for (i in seq_along(out.vars)) new.x[,i+1] <- x[[ out.vars[i] ]][run]
  } 
  
  names(new.x)[2:ncol(new.x)] <- out.vars.n
  
  # Calculate group totals
  new.x$num <- rowSums(new.x[,setdiff(grep("num", names(new.x)), 
                                      grep(".g2", names(new.x)))])
  if (groups == 2) {
    new.x$num.g2 <- rowSums(new.x[,intersect(grep("num", names(new.x)), 
                                             grep(".g2", names(new.x)))])
  }
  
  return(new.x)
}


#' @title Extract Model Run for Stochastic Models
#'
#' @description This function extracts model simulations for objects of classes
#'   \code{epiICM} and \code{epiNet.simTrans} into a data frame using the 
#'   generic \code{as.data.frame} function. 
#'
#' @param x an \code{EpiModel} object of class \code{epiICM} or \code{epiNet.simTrans}.
#' @param sim simulation number from model, used only if \code{out="vals"}.
#' @param out type of output to data frame: use \code{"mean"} for mean across 
#'   simulations, \code{"sd"} for standard deviations across all simulations, 
#'   and \code{"vals"} for values from one simulation (with simulation set 
#'   with \code{sim}).
#' @param row.names see \code{\link{as.data.frame.default}}.
#' @param optional see \code{\link{as.data.frame.default}}.
#' @param ...  see \code{\link{as.data.frame.default}}.
#' 
#' @details
#' These generics work from both \code{epiICM} and \code{epiNet.simTrans} class
#' objects, and output data frames containing the time-specific means, 
#' standard deviations, or model values (compartment and flow sizes from individual
#' simulations) from these stochastic model classes. For summary output types, 
#' row mean or standard deviations are calculated across all simulations. 
#' Total population size for one and two-group models are also calculated and appended.
#' 
#' @method as.data.frame epiICM
#' @keywords extract
#' @export
#' 
#' @examples
#' # Individual contact SIS model with 2 simulations
#' mod <- epiICM(type="SIS", s.num=500, i.num=1, 
#'               trans.rate=0.8, act.rate=2, rec.rate=0.1, 
#'               nsteps=20, nsims=2)
#'               
#' # Default output is mean across simulations
#' as.data.frame(mod)
#' 
#' # Standard deviations of simulations
#' as.data.frame(mod, out="sd")
#' 
#' # Individual simulation runs, with default sim=1
#' as.data.frame(mod, out="vals")
#' as.data.frame(mod, out="vals", sim=2)
#' 
as.data.frame.epiICM <- function(x, 
                                 row.names = NULL, 
                                 optional = FALSE, 
                                 sim,
                                 out = "mean",
                                 ...) {
  
  new.x <- data.frame(time=x$time)
  nsims <- x$nsims
  nts <- max(x$time)
  if (class(x) == "epiICM") groups <- x$groups
  if (class(x) == "epiNet.simTrans") groups <- x$modes
  
  out.vars <- c(grep("num", names(x)), grep("flow", names(x)))
  out.vars.n <- names(x)[out.vars]
  
  if (out == "vals") {
    if (missing(sim)) sim <- 1

    # Output for models with 1 sim
    if (nsims == 1) {
      if (sim > 1) stop("Specify sim = 1")
      for (i in seq_along(out.vars)) new.x[,i+1] <- x[[ out.vars[i] ]]
    }
    
    # Output for models with multiple sims
    if (nsims > 1) {
      if (sim > nsims) {
        stop(paste("Specify sim between 1 and", nsims))
      }
      for (i in seq_along(out.vars)) new.x[,i+1] <- x[[ out.vars[i] ]][sim]
    }
    names(new.x)[2:ncol(new.x)] <- out.vars.n
    # Calculate group totals
    if (class(x) == "epiICM") {
      new.x$num <- rowSums(new.x[,setdiff(grep("num", names(new.x)), 
                                          grep(".g2", names(new.x)))])
    }
    if (class(x) == "epiNet.simTrans") {
      new.x$num <- rowSums(new.x[,setdiff(grep("num", names(new.x)), 
                                          grep(".m2", names(new.x)))])
    }
    if (groups == 2) {
      if (class(x) == "epiICM") {
        new.x$num.g2 <- rowSums(new.x[,intersect(grep("num", names(new.x)), 
                                                 grep(".g2", names(new.x)))])
      }
      if (class(x) == "epiNet.simTrans") {
        new.x$num.m2 <- rowSums(new.x[,intersect(grep("num", names(new.x)), 
                                                 grep(".m2", names(new.x)))])
      }
    }
  }
  
  ## Output means
  if (out == "mean") {
    if (nsims == 1) {
      for (i in seq_along(out.vars)) new.x[,i+1] <- x[[ out.vars[i] ]]
    }
    if (nsims > 1) {
      for (i in seq_along(out.vars)) new.x[,i+1] <- rowMeans(x[[ out.vars[i] ]])
    }
    names(new.x)[2:ncol(new.x)] <- out.vars.n
    # Calculate mean of group totals
    if (class(x) == "epiICM") {
      numvars <- Reduce("+", x[setdiff(grep("num", names(x)), grep(".g2", names(x)))])
      if (nsims == 1) new.x$num <- numvars
      if (nsims > 1) new.x$num <- rowMeans(numvars)
    }
    if (class(x) == "epiNet.simTrans") {
      numvars <- Reduce("+", x[setdiff(grep("num", names(x)), grep(".m2", names(x)))])
      if (nsims == 1) new.x$num <- numvars
      if (nsims > 1) new.x$num <- rowMeans(numvars)
    }
    if (groups == 2) {
      if (class(x) == "epiICM") {
        numvars <- Reduce("+", x[intersect(grep("num", names(x)), grep(".g2", names(x)))])
        if (nsims == 1) new.x$num.g2 <- numvars
        if (nsims > 1) new.x$num.g2 <- rowMeans(numvars)
      }
      if (class(x) == "epiNet.simTrans") {
        numvars <- Reduce("+", x[intersect(grep("num", names(x)), grep(".m2", names(x)))])
        if (nsims == 1) new.x$num.m2 <- numvars
        if (nsims > 1) new.x$num.m2 <- rowMeans(numvars)
      }
    } 
  }
  
  ## Output standard deviations
  if (out == "sd") {
    if (nsims == 1) {
      for (i in seq_along(out.vars)) new.x[,i+1] <- 0
    }
    if (nsims > 1) {
      for (i in seq_along(out.vars)) new.x[,i+1] <- apply((x[[ out.vars[i] ]]), 1, sd)
    }
    names(new.x)[2:ncol(new.x)] <- out.vars.n
    # Calculate standard deviations of group totals
    if (class(x) == "epiICM") {
      numvars <- Reduce("+", x[setdiff(grep("num", names(x)), grep(".g2", names(x)))])
      if (nsims == 1) new.x$num <- 0
      if (nsims > 1) new.x$num <- apply(numvars, 1, sd)
    }
    if (class(x) == "epiNet.simTrans") {
      numvars <- Reduce("+", x[setdiff(grep("num", names(x)), grep(".m2", names(x)))])
      if (nsims == 1) new.x$num <- 0
      if (nsims > 1) new.x$num <- apply(numvars, 1, sd)
    }
    if (groups == 2) {
      if (class(x) == "epiICM") {
        numvars <- Reduce("+", x[intersect(grep("num", names(x)), grep(".g2", names(x)))])
        if (nsims == 1) new.x$num.g2 <- 0
        if (nsims > 1) new.x$num.g2 <- apply(numvars, 1, sd)
      }
      if (class(x) == "epiNet.simTrans") {
        numvars <- Reduce("+", x[intersect(grep("num", names(x)), grep(".m2", names(x)))])
        if (nsims == 1) new.x$num.m2 <- 0
        if (nsims > 1) new.x$num.m2 <- apply(numvars, 1, sd)
      }
    } 
  }
  
  return(new.x)
}


#' @method as.data.frame epiNet.simTrans
#' @export
#' @rdname as.data.frame.epiICM
as.data.frame.epiNet.simTrans <- function(x, 
                                          row.names = NULL, 
                                          optional = FALSE, 
                                          sim,
                                          out = "mean",
                                          ...) {
  
  new.x <- as.data.frame.epiICM(x, 
                                row.names = row.names, 
                                optional = optional, 
                                sim = sim,
                                out = out,
                                ...) 
  return(new.x)
  
}