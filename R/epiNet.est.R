#' 
#' @title Stochastic Network Estimation with STERGM
#'
#' @description Estimates and diagnoses a statistical network
#'   model using the exponential random graph modeling (ERGM) framework with
#'   extensions for dynamic/temporal models (STERGM)
#'
#' @param nw an object of class \code{network}.
#' @param formation a right-hand sided STERGM formation formula in the form 
#'   \code{~ edges + ...}, where \code{...} are additional network statistics.
#' @param dissolution a right-hand sided STERGM dissolution formula in the form
#'   \code{~ offset(edges) + ...}, where \code{...} are additional network 
#'   statistics.
#' @param target.stats a vector of target statistics for formation model, with
#'   one number for each network statistic in the model (see \code{\link{stergm}}).
#' @param coef.diss an object of class \code{dissolution.coefs} that is calculated
#'   with the function of that name containing the target coefficients to target.
#' @param constraints a right-hand sided formula specifying contraints for the
#'   modeled network, in the form \code{~...}, where \code{...} are constraint
#'   terms described in \code{\link{stergm}}. By default, no constraints are set.
#' @param edapprox if \code{TRUE}, use the indirect edges dissolution approximation 
#'   method for dynamic modeling fitting rather than the more time-intensive full
#'   STERGM estimation (available with \code{FALSE}). For this to function,
#'   the \code{coef.diss} terms should be in the same order as the formation 
#'   formula terms.
#' @param save.stats if \code{TRUE}, run simulation-based model diagnostics for 
#'   the formation and dissolution processes.
#' @param stats.formula a right-hand sided ERGM formula that includes network 
#'   statistics of interest, with the default to the formation formula terms.
#' @param stats.start starting timestep for calculating the summary measure 
#'   diagnostics.
#' @param stats.end ending timestep for calculating the summary measure diagnostics.
#' @param verbose if \code{TRUE}, print simulation progress to the console.
#' @param ... additional arguments to pass to \code{ergm} or \code{stergm}.
#' 
#' @details
#' This function is a wrapper function around \code{ergm} and \code{stergm} to
#' facilitate estimating a dynamic/temporal ERGM from target statistics. This
#' framework is suggested for parameterizing a stochastic network epidemic
#' model from empirical partnership formation and dissolution data collected
#' in an egocentric format (i.e., study subjects asked about their recent
#' partners, including questions on partnership duration).
#'   
#' Additional functionality here includes time-varying diagnostics for a temporal 
#' simulation from the fitted model, to ensure that target statistics are
#' approximated by the network simulation. This is available through the 
#' \code{save.stats} argument, with the \code{stats.formula} argument allowing 
#' flexibility in network statistics of interest outside those in the formation
#' formula. The diagnostics are available by printing and plotting the fitted
#' \code{epiNet.est} object.
#'   
#' With a fitted and diagnosed network model, one may proceed to 
#' \code{\link{epiNet.simNet}} for simulating a series of networks for use in
#' an independent stochastic network epidemic model, or straight to the epidemic 
#' model function \code{\link{epiNet.simTrans}} if there are dependencies between
#' the network model and the disease model. See Section 4 of 
#' \href{../doc/Tutorial.pdf}{EpiModel Tutorial} for further clarification.
#' 
#' @references
#' Krivitsky PN, Handcock MS (2010). A Separable Model for Dynamic Networks. 
#'   \url{http://arxiv.org/abs/1011.1937}.
#' 
#' @keywords model
#' @seealso \code{\link{stergm}}, \code{\link{epiNet.simNet}}, 
#'   \code{\link{epiNet.simTrans}}
#' @export
#' 
#' @examples
#' \dontrun{
#' # Initialize an empty network of 500 nodes divided equally into two races
#' nw <- network.initialize(n=500, directed=F)  
#' nw <- set.vertex.attribute(nw, "race", value = rep(0:1, each=250))
#' 
#' # Set formation and dissolution formulas
#' formation <- ~ edges + nodematch("race") + degree(0) + concurrent
#' dissolution <- ~ offset(edges)
#' 
#' # Set target statistics for formation 
#' target.stats <- c(225, 187, 180, 90)
#' 
#' # Set a vector of durations (here in months, but timestep units are arbitrary)
#' #  Then use the dissolution.coefs functions to obtain the offset coefficients
#' durations <- 20
#' coef.diss <- dissolution.coefs(dissolution, durations)
#' coef.diss
#'
#' # Set the stats.formula to include more degree tergms than the formation formula
#' dx.stats <- ~edges + nodematch("race") + degree(0:5)
#' 
#' # Estimate the STERGM with all the parameters set above with diagnostic
#' #  simulations through 2000 timesteps, using the edges dissolution approximation
#' est <- epiNet.est(
#'   nw, 
#'   formation, 
#'   dissolution, 
#'   target.stats, 
#'   coef.diss,
#'   edapprox = TRUE,
#'   stats.formula = dx.stats,
#'   dx.end = 2000)
#' 
#' # Estimate the STERGM directly and use the default for model statistics
#' est2 <- epiNet.est(
#'   nw,
#'   formation,
#'   dissolution,
#'   target.stats,
#'   coef.diss,
#'   edapprox = FALSE)
#' est2
#' 
#' # Print and plot summaries
#' est
#' plot(est)
#' plot(est, type="duration")
#' }
#' 
epiNet.est <- function(nw,
                       formation,
                       dissolution,
                       target.stats, 
                       coef.diss,
                       constraints,
                       edapprox = TRUE,
                       save.stats = TRUE, 
                       stats.formula, 
                       stats.start = 1,
                       stats.end = 1000,
                       verbose = TRUE,
                       ...) {
  
  ###
  
	formation.nw <- update(formation, nw ~.)
  environment(formation.nw) <- environment()
  
  if (missing(constraints)) constraints	<- ~. 
  environment(constraints) <- environment()
  
  ###

  if (edapprox == FALSE) {
    if (verbose == TRUE) {
      cat("\n======================")
      cat("\nFitting STERGM")
      cat("\n======================\n")
    }
    
    fit <- stergm(nw,
                  formation = formation,
                  dissolution = dissolution,
                  targets = "formation",
                  target.stats = target.stats,
                  offset.coef.diss = coef.diss$coef.crude,
                  constraints = constraints, 
                  estimate = "EGMME", 
                  ...)
    
    coef.form <- fit$formation.fit
    
    result <- list()
    result$fit <- fit
    result$formation <- formation
    result$coef.form <- coef.form$coef
    result$dissolution <- dissolution
    result$coef.diss <- coef.diss
    result$constraints <- constraints
    result$edapprox <- edapprox
    
  } else {
    if (verbose == TRUE) {
      cat("\n======================")
      cat("\nFitting ERGM")
      cat("\n======================\n")
    }
        
    fit <- ergm(formation.nw, 
                target.stats=target.stats, 
                constraints=constraints,
                ...)

    ## TODO: enforce a check here to make sure dissolution formula embedded in formation
    coef.form <- fit$coef 
    for (i in 1:length(coef.diss$coef.crude)) {
      coef.form[i] <- coef.form[i] - coef.diss$coef.crude[i]
    }
    
    result <- list()
    result$fit <- fit
    result$formation <- formation
    result$coef.form <- coef.form
    result$dissolution <- dissolution
    result$coef.diss <- coef.diss
    result$constraints <- constraints
    result$edapprox <- edapprox
    
  }
  
  if (save.stats == TRUE) {
    
    if (verbose == TRUE) {
      cat("\n======================")
      cat("\nRunning Diagnostics")
      cat("\n======================\n")
    }
    
    if (missing(stats.formula)) stats.formula <- formation
    
    if (verbose == TRUE) cat("* Simulating network\n")
    if (edapprox == FALSE) {
      diag.sim <- simulate(fit, 
                           time.slices = stats.end,
                           monitor = stats.formula)
      diag.sim.ts <- simulate(fit, 
                              time.slices = 1,
                              monitor = formation,
                              output = "stats")
    } else {
      fit.sim <- simulate(fit)
      diag.sim <- simulate(fit.sim,
                           formation = formation, 
                           dissolution = dissolution,
                           coef.form = coef.form, 
                           coef.diss = coef.diss$coef.crude,
                           time.slices = stats.end, 
                           constraints = constraints,
                           monitor = stats.formula)
      diag.sim.ts <- simulate(nw, 
                              formation = formation, 
                              dissolution = dissolution,
                              coef.form = coef.form, 
                              coef.diss = coef.diss$coef.crude, 
                              constraints = constraints,
                              monitor = formation,
                              output = "stats")
    }
    
    if (verbose == TRUE) cat("* Calculating formation statistics\n")
    # Extract model stats and calculate mean/sd
    result$sim.stats <- attributes(diag.sim)$stats
    sim.stats <- attributes(diag.sim)$stats[stats.start:stats.end, ]
    if (class(sim.stats) == "numeric") {
      stats.means <- mean(sim.stats)
      stats.sd <- sd(sim.stats)
      stats <- data.frame(sorder=1, 
                          names=attributes(attributes(diag.sim)$stats)$dimnames[[2]],
                          stats.means=stats.means,
                          stats.sd=stats.sd)
    } else {
      stats.means <- colMeans(sim.stats)
      stats.sd <- apply(sim.stats, 2, sd)
      stats <- data.frame(sorder=1:length(names(stats.means)), 
                          names=names(stats.means), 
                          stats.means, stats.sd)
    }
    
    # Get stats form for target statistics
    ts.out <- data.frame(names=attributes(diag.sim.ts)$dimnames[[2]], 
                         targets=target.stats)
    result$target.stats <- ts.out
    
    # Create stats table for output
    stats.table <- merge(ts.out, stats, all=T)
    stats.table <- stats.table[do.call("order",
                      stats.table[, "sorder", drop = FALSE]), , drop = FALSE]
    rownames(stats.table) <- stats.table$names
    stats.formation <- stats.table[,-c(1,3)]
    result$stats.formation <- stats.formation
    
    if (verbose == TRUE) cat("* Calculating duration statistics\n\n")
    # Duration calculations
    sim.df <- as.data.frame(diag.sim)
    result$edgelist <- sim.df
    
    # Create duration table for "dissolution = ~ offset(edges)"
    if (dissolution == ~offset(edges)) {
      duration.mean <- mean(sim.df$duration)
      duration.sd <- sd(sim.df$duration)
      duration.expected <- exp(coef.diss$coef.crude[1]) + 1
      result$stats.duration <- c(target.dur = duration.expected,
                                 sim.mean.dur = duration.mean,
                                 sim.sd.dur = duration.sd)
    } else {
      stop('Only ~offset(edges) dissolution models currently supported')
    }
    
#     # Create duration table for "dissolution = ~ offset(edges) + offset(nodematch("foo"))
#     if (length(all.names(dissolution)) > 3 && any(all.names(dissolution) == "nodematch")) {
#       dform <- strsplit(as.character(dissolution), "[+]")[[2]]
#       dform2 <- dform[grep("nodematch", dform)]
#       nmatchvar <- strsplit(dform2, "[\"]")[[1]][2]
#       nmatchvals <- nw %v% nmatchvar
#       sim.df$nmatch <- nmatchvals[sim.df$head] == nmatchvals[sim.df$tail]
#       duration.mean.nomatch <- mean(sim.df$duration[sim.df$nmatch == FALSE])
#       duration.mean.match <- mean(sim.df$duration[sim.df$nmatch == TRUE])
#       duration.sd.nomatch <- sd(sim.df$duration[sim.df$nmatch == FALSE])
#       duration.sd.match <- sd(sim.df$duration[sim.df$nmatch == TRUE])
#       
#       expected.dur <- c(nomatch=exp(coef.diss[1]) + 1, 
#                         match=exp(coef.diss[1] + coef.diss[2]) + 1)
#       result$stats.duration <- cbind(expected.dur,
#                                 sim.mean.dur=c(duration.mean.nomatch, duration.mean.match),
#                                 sim.mean.sd=c(duration.sd.nomatch, duration.sd.match))                 
#     }
    
    # Calculate mean partnership age from edgelist
    result$pages <- edgelist.meanage(el=sim.df, dissolution=dissolution, nw=fit$network)
    
    
  } # End diagnostic stats condition
  
  class(result) <- "epiNet.est"
  return(result) 
}

