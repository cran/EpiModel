#' 
#' @title Simulate Networks from STERGM Model Fit
#'
#' @description This function simulates a set of dynamic networks of class
#'  \code{networkDynamic} given a STERGM model fit with \code{\link{epiNet.est}}. 
#'
#' @param x an \code{EpiModel} object of class \code{\link{epiNet.est}}.
#' @param nsteps number of time steps over which to simulate relational dynamics.
#' @param nsims number of network simulations.
#' @param save.stats if \code{TRUE}, run simulation-based model diagnostics for the
#'   formation and dissolution processes.
#' @param stats.formula a right-hand sided ERGM formula that includes 
#'   network statistics of interest, with the default to the formation formula 
#'   terms.
#' @param verbose if \code{TRUE}, print simulation progress to the console.
#' @param ... additional arguments to pass to \code{simulate.stergm}.
#'
#' @details
#' This function is a wrapper function around the \code{simulate} functions in the 
#' \code{tergm} package. The purpose is to simulate a series of dynamic relational
#' networks for a set number of time steps, and to provide access to access to 
#' diagnostics on those network simulations. 
#' 
#' The intended use of these simulated networks is in independent network disease
#' models simulated in \code{\link{epiNet.simTrans}}. With dependent network 
#' models, one bypasses this network simulation step as that processes is 
#' coincident with the disease simulation. See Section 4 of 
#' \href{../doc/Tutorial.pdf}{EpiModel Tutorial} for further clarification.
#'
#' @seealso \code{\link{plot.epiNet.simNet}}, \code{\link{epiNet.est}}, 
#'  \code{\link{epiNet.simTrans}}
#' @keywords model
#' @export
#'
#' @examples
#' \dontrun{
#' # Initialize bipartite network
#' num.m1 <- 100 
#' num.m2 <- 100
#' nw <- network.initialize(num.m1 + num.m2, 
#'                          bipartite=num.m1, directed=FALSE)
#' 
#' # Set mode-specific fractional degree distributions 
#' deg.dist.m1 <- c(0.40, 0.55, 0.04, 0.01)
#' deg.dist.m2 <- c(0.48, 0.41, 0.08, 0.03)
#' 
#' # Set model parameters and control arguments
#' formation <- ~ edges + b1degree(0:1) + b2degree(0:1)
#' target.stats <- c(66, 40, 55, 48, 41)
#' dissolution <- ~ offset(edges)
#' duration <- 25
#' coef.diss <- dissolution.coefs(dissolution, duration)
#' 
#' # Estimate the model with all the parameters set above
#' est <- epiNet.est(
#'   nw, 
#'   formation, 
#'   dissolution, 
#'   target.stats, 
#'   coef.diss)
#' est
#' 
#' # Simulate 5 networks over 100 timesteps
#' nwsims <- epiNet.simNet(est, nsteps = 100, nsims = 5)
#' 
#' # Print the object, plot the formation summary, plot the duration summary
#' nwsims
#' plot(nwsims, sim=2)
#' plot(nwsims, type = "duration", sim = 2)
#' }
#' 
epiNet.simNet <- function(x, 
                          nsteps, 
                          nsims = 1,
                          save.stats = TRUE,
                          stats.formula,
                          verbose = TRUE,
                          ...) {
  
  if (verbose == TRUE) {
    if (nsims == 1) {
      cat("\n======================")
      cat("\nSimulating 1 Network")
      cat("\n======================\n")
    } else {
      cat("\n======================")
      cat("\nSimulating ", nsims, " Networks", sep="")
      cat("\n======================\n")
    }
  }
  
  nwsims <- list()
  nwstats <- list()
  nwpages <- list()
  for (s in 1:nsims) {
  
    ## Extract base network information
    base.nw <- x$fit$network
    formation <- x$formation
    dissolution <- x$dissolution
    coef.form <- x$coef.form
    coef.diss <- x$coef.diss$coef.crude
    constraints <- x$constraints
    modes <- ifelse(base.nw %n% "bipartite", 2, 1)
    
    ## Simulate base network to get non-empty network
    nw <- simulate(x$fit)
    if (class(x$fit) == "stergm") {
      nw <- network.collapse(nw, at=1)
    } 
    base.nw <- activate.vertices(base.nw, onset=1, terminus=Inf)
    base.nw <- activate.edges(base.nw, onset=1, terminus=Inf)
    
    ## Default network statistics is formation formula
    if (missing(stats.formula)) stats.formula <- formation
    
    ## Simulate dynamic network
    nwsim <- simulate(nw,
                      formation = formation, 
                      dissolution = dissolution,
                      coef.form = coef.form, 
                      coef.diss = coef.diss,
                      constraints = constraints,
                      time.start = 1,
                      time.slices = nsteps,
                      time.offset = 0,
                      monitor = stats.formula, 
                      ...)
  
    ## Save simulate out to output list
    nwsims[[s]] <- nwsim
    names(nwsims)[[s]] <- paste("nwsim", s, sep="")
    
    ## Save simulated network statistics separately
    if (save.stats == TRUE) {
      
      # Store out summary stats separately
      nwstats[[s]] <- as.data.frame(attributes(nwsim)$stats)
      
      el <- as.data.frame(nwsim)
      pages <- edgelist.meanage(el=el, dissolution=dissolution, nw=base.nw)
      nwpages[[s]] <- pages
    }
  
    if (verbose == TRUE) cat("SIM = ", s, "/", nsims, "\n", sep="")
  }
  if (verbose == TRUE) cat("\n")
  
  ## Other data to save out to highest object level
  nwsims$base.nw <- x$fit$network
  nwsims$nsteps <- nsteps
  nwsims$nsims <- nsims
  nwsims$formation <- formation
  nwsims$dissolution <- dissolution
  nwsims$coef.form <- coef.form
  nwsims$coef.diss <- coef.diss
  nwsims$constraints <- constraints
  nwsims$target.stats <- x$target.stats
  if (save.stats == TRUE) {
    nwsims$stats <- nwstats
    nwsims$pages <- nwpages
  }
  
  ## List element names for simulated networks
  simnames <- paste("sim", 1:nsims, sep="")
  if (save.stats == TRUE) {
    names(nwsims$stats) <- simnames
    names(nwsims$pages) <- simnames
  }
  
  class(nwsims) <- "epiNet.simNet"
  return(nwsims)
}

