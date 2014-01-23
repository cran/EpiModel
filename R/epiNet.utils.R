

# Exported Functions ------------------------------------------------------

#' @title Vertex Attributes for Bipartite Network
#'
#' @description This function outputs static vertex attributes for a 
#'  bipartite network for one specified mode.
#' 
#' @param nw an object of class \code{network} or \code{networkDynamic}.
#' @param mode mode number to extract values for.
#' @param val static attribute values to return.
#' 
#' @export
#' @keywords epiNetUtils internal
#' 
#' @examples
#' nw <- network.initialize(10, bipartite=5)
#' nw %v% "sex" <- rep(c("M", "F"), each=5)
#' bipvals(nw, mode=1, "sex")
#' 
bipvals <- function(nw, 
                    mode, 
                    val
                    ) {
  
  if (!is.numeric(nw$gal$bipartite))
    stop("nw must be a bipartite network")
  if (missing(mode))
    stop("Specify mode=1 or mode=2")
  
  nw %s% modeids(nw, mode) %v% val
}


#' @title Check Degree Distribution for Bipartite Target Statistics
#'
#' @description This function checks for consistency in the degree distributions
#' between the two modes in a bipartite network, and if balance is present returns
#' the implied network statistics for use in a ERGM.
#' 
#' 
#' @param num.m1 number of vertices in mode 1.
#' @param num.m2 number of vertices in mode 2.
#' @param deg.dist.m1 vector with fractional degree distribution for mode 1.
#' @param deg.dist.m2 vector with fractional degree distribution for mode 2.
#' 
#' @details
#' This function outputs the number of nodes of degree 0 to m, where m is the
#' length of a fractional degree distribution vector, given that vector and the
#' size of the mode. This utility is used to check for balance in implied degree
#' given that fractional distribution within bipartite network simulations, in which
#' the degree-contrained counts must match across modes. 
#' 
#' @seealso
#' For a detailed explanation of this function, see the included HTML vignette: 
#' \href{../doc/epiNetUtils.html}{EpiModel Network Utility Functions.}
#' 
#' @export
#' @keywords epiNetUtils
#' 
#' @examples
#' # An imbalanced distribution          
#' bip.degdist.check(num.m1 = 500, num.m2 = 500, 
#'                   deg.dist.m1 = c(0.40, 0.55, 0.03, 0.02),
#'                   deg.dist.m2 = c(0.48, 0.41, 0.08, 0.03))
#'                   
#' # A balanced distribution
#' targets <- bip.degdist.check(num.m1 = 500, num.m2 = 500, 
#'                   deg.dist.m1 = c(0.40, 0.55, 0.04, 0.01),
#'                   deg.dist.m2 = c(0.48, 0.41, 0.08, 0.03))
#' targets
#' 
bip.degdist.check <- function(num.m1, 
                              num.m2, 
                              deg.dist.m1, 
                              deg.dist.m2
                              ) {
  
  num <- num.m1 + num.m2
  
  deg.counts.m1 <- deg.dist.m1*num.m1
  deg.counts.m2 <- deg.dist.m2*num.m2
  
  tot.deg.m1 <- sum(deg.counts.m1 * (1:length(deg.dist.m1)-1))
  tot.deg.m2 <- sum(deg.counts.m2 * (1:length(deg.dist.m2)-1))
  
  mat <- matrix(c(deg.dist.m1, deg.counts.m1, 
                  deg.dist.m2, deg.counts.m2), ncol=4)
  mat <- rbind(mat, c(sum(deg.dist.m1), tot.deg.m1, sum(deg.dist.m2), tot.deg.m2))
  
  colnames(mat) <- c("m1.dist", "m1.cnt", "m2.dist", "m2.cnt")
  rownames(mat) <- c(paste("Deg", 0:(length(deg.dist.m1)-1), sep=""), "TOTAL")
  
  cat("Bipartite Degree Distribution Check\n")
  cat("=============================================\n")
  print(mat, print.gap=3)
  cat("=============================================\n")
  
  if (sum(deg.dist.m1) != 1 | sum(deg.dist.m2) != 1 | round(tot.deg.m1) != round(tot.deg.m2)) {
    if (sum(deg.dist.m1) != 1) cat("** deg.dist.m1 TOTAL != 1 \n")
    if (sum(deg.dist.m2) != 1) cat("** deg.dist.m2 TOTAL != 1 \n")
    
    if (round(tot.deg.m1) != round(tot.deg.m2)) cat("** m1.cnt TOTAL != m2.cnt TOTAL \n")
  } else {
    cat("** distributions balanced \n")
  }
  invisible(c(tot.deg.m1, deg.counts.m1, deg.counts.m2))
}


#' @title Calcuate Dissolution Coefficients for Stochastic Network Models
#'
#' @description This function calculates dissolution coefficients from a 
#'   duration to pass as terms to an ERGM/STERGM fit and applies a coefficient
#'   correction to account for deaths.
#' 
#' @param dissolution a right-hand sided STERGM dissolution formula 
#'  (see \code{\link{epiNet.est}}); currently limited only to simple
#'  \code{~offset(edges)} dissolution models.
#' @param duration an edge duration in arbitrary time units.
#' @param d.rate background death rate in the absence of disease.
#'   
#' @details
#' This function performs two duties prior to the estimation phase of \code{epiNet}
#' class stochastic network models.
#' \enumerate{
#'  \item \strong{Transforming duration to coefficient form:} an average duration
#'    of edges in a population must be mathematically transformed to 
#'    logit coefficients. The theory and mathematics of the transformation are
#'    further developed in the \code{networkDynamic} package vignette.
#'  \item \strong{Adjusting coefficients to account for death:} In a dynamic 
#'    network simulation in an open population (in which there are births and 
#'    deaths), it is necessary to adjust the dissolution coefficient for the 
#'    STERGM simulations to account for the death as a competing risk to 
#'    edge dissolution.
#' }
#'
#' Future releases of this software will allow for more flexibility in the possible
#' dissolution models that may be calculated here, including models that have 
#' heterogenous dissolution probabilities conditional on nodal or edge attributes.
#' 
#' @return
#' A list of with the following elements: 
#' 
#' @seealso 
#' For a detailed explanation of this function, please see the HTML vignette: 
#' \href{../doc/epiNetUtils.html}{EpiModel Network Utility Functions}. 
#'   
#' @export
#' @keywords epiNetUtils
#'
#' @examples
#' dissolution <- ~offset(edges)
#' duration <- 25
#' dissolution.coefs(dissolution, duration)
#' dissolution.coefs(dissolution, duration, d.rate=0.001)
#'
dissolution.coefs <- function(dissolution,
                              duration, 
                              d.rate = 0
                              ) {
  
  
  # Check form of dissolution formula
  form.length <- length(strsplit(as.character(dissolution)[2], "[+]")[[1]])
  t1.edges <- grepl("offset[(]edges", 
                    strsplit(as.character(dissolution)[2], "[+]")[[1]][1])
  
  # Log transformation of duration to coefficent  
  if (t1.edges == TRUE && form.length == 1) {
    coef.diss <- log(duration[1] - 1)
  } else {
    stop('Only ~offset(edges) dissolution models currently supported')
  }
  
  if (d.rate > 0) {
    # Exogenous death correction to coefficient
    exp.dur <- 1 + exp(coef.diss) 
    prob.diss <- 1 / exp.dur
    
    prob.neither.dying <- (1 - d.rate)^2
    prob.either.dying <- 2*d.rate - d.rate^2
    
    prob <- 1 - ((prob.diss - prob.either.dying) / prob.neither.dying)
    if (prob >= 1) stop('The competing risk of mortality is too high for the given duration. Specify a lower d.rate')
    coef.diss.adj <- logit(prob)
  } else {
    coef.diss.adj <- coef.diss
  }
  
  out <- list()
  out$dissolution <- dissolution
  out$duration <- duration
  out$coef.adj <- coef.diss.adj
  out$coef.crude <- coef.diss
  
  class(out) <- 'dissolution.coefs'
  return(out)
}




#' @title Table of Partnership Censoring
#'
#' @description This function outputs a table of the number and percent
#'  of partnerships that are left, right, both, or uncensored for a
#'  \code{networkDynamic} object.
#' 
#' @param el a timed edgelist with start and end times from extracted from
#' a \code{networkDynamic} object using the \code{as.data.frame.networkDynamic}
#' function.
#' 
#' @export
#' @keywords epiNetUtils
#' 
#' @details
#' Given a STERGM simulation over a specified number of time steps, the edges
#' within that simulation may be left-censored (started before the first step),
#' right-censored (continued after the last step), right and left-censored, or
#' uncensored. The amount of censoring will tend to increase when the average
#' partnership duration approaches the length of the simulation.
#' 
#' @examples
#' nw <- network.initialize(n = 100, directed = FALSE)
#' nw <- set.vertex.attribute(nw, "race", value = rep(0:1, each = 50))
#' formation <- ~ edges + nodematch("race")
#' target.stats <- c(45, 25)
#' dissolution <- ~ offset(edges)
#' coef.diss <- dissolution.coefs(dissolution, duration = 20)
#' 
#' # Model estimation, then simulation
#' est <- epiNet.est(nw,
#'                   formation,
#'                   dissolution,
#'                   target.stats,
#'                   coef.diss, 
#'                   save.stats = FALSE,
#'                   verbose = FALSE)
#' nwsims <- epiNet.simNet(est, nsteps = 100, 
#'                         nsims = 1, verbose = FALSE)
#' 
#' # Extract edgelist from first simulation, then calculate censoring
#' el.sim1 <- as.data.frame(nwsims$nwsim1)
#' edgelist.censor(el.sim1)
#'
edgelist.censor <- function(el) {
  
  time.steps <- max(el$terminus)
  min.step <- min(el$onset)
  
  # left censored
  leftcens <- el$onset == min.step
  leftcens.num <- sum(leftcens)
  leftcens.pct <- leftcens.num/nrow(el)
  
  # right censored
  rightcens <- el$terminus == time.steps
  rightcens.num <- sum(rightcens)
  rightcens.pct <- rightcens.num/nrow(el)
  
  # partnership lasts for entire window (left and right censored)
  lrcens <- el$onset == min.step & el$terminus == time.steps
  lrcens.num <- sum(lrcens)
  lrcens.pct <- lrcens.num/nrow(el)
  
  # fully observed
  nocens <- el$onset > min.step & el$terminus < time.steps
  nocens.num <- sum(nocens)
  nocens.pct <- nocens.num/nrow(el)
  
  ## Table
  nums <- rbind(leftcens.num, rightcens.num, lrcens.num, nocens.num)
  pcts <- rbind(leftcens.pct, rightcens.pct, lrcens.pct, nocens.pct)  
  out <- cbind(nums, pcts)
  rownames(out) <- c("Left Cens.", "Right Cens.", "Both Cens.", "No Cens.")
  colnames(out) <- c("num", "pct")
  
  return(out)
}


#' @title Mean Age of Partnerships over Time
#'
#' @description This function outputs a vector of mean ages of partnerships
#'   at a series of timesteps
#' 
#' @param x an \code{EpiModel} object of class \code{\link{epiNet.est}}.
#' @param nw if not passing \code{x}, a \code{networkDynamic} object.
#' @param el if not passing \code{x}, a timed edgelist from a \code{networkDynamic}
#' object extracted with the \code{as.data.frame.networkDynamic} function.
#' @param dissolution if not passing \code{x}, a right-hand sided STERGM 
#'  dissolution formula (see \code{\link{epiNet.est}}).
#' 
#' @details
#' This function calculates the mean partnership age at each time step over 
#' a dynamic network simulation from \code{\link{epiNet.est}}. These objects 
#' contain the network, edgelist, and dissolution objects needed for the 
#' calculation. Alternatively, one may pass in these objects separately if
#' \code{epiNet.est} was not used, or statistics were not run requested after
#' the estimation. 
#' 
#' Currently, the calculations are limited to those dissolution formulas with a single 
#' homogenous dissolution (\code{~offset(edges)}). This functionality will be 
#' expanded in future releases.
#' 
#' @export
#' @keywords epiNetUtils internal
#' 
#' @examples
#' nw <- network.initialize(n = 100, directed = FALSE)
#' nw <- set.vertex.attribute(nw, "race", value = rep(0:1, each = 50))
#' formation <- ~ edges + nodematch("race")
#' target.stats <- c(45, 25)
#' dissolution <- ~ offset(edges)
#' coef.diss <- dissolution.coefs(dissolution, duration = 20)
#' 
#' # Model estimation
#' est <- epiNet.est(
#'   nw,
#'   formation,
#'   dissolution,
#'   target.stats,
#'   coef.diss, 
#'   save.stats = FALSE,
#'   verbose = FALSE)
#' 
#' # Get the formation coefficients, then simulate a dynamic network
#' # from a simulated static network
#' coef.form <- est$coef.form
#' sim <- simulate(
#'   simulate(est$fit),
#'   formation = formation, 
#'   dissolution = dissolution,
#'   coef.form = coef.form, 
#'   coef.diss = coef.diss[[3]],
#'   time.slices = 200)
#' 
#' # Gather the three objects needed for the age calculation
#' nw <- est$fit$network
#' el <- as.data.frame(sim)
#' dissolution <- est$dissolution
#' edgelist.meanage(nw=nw, el=el, dissolution=dissolution)
#' 
#' # Alternatively, epiNet.est automatically calculates these with stats = TRUE
#' est <- epiNet.est(
#'   nw, 
#'   formation, 
#'   dissolution, 
#'   target.stats, 
#'   coef.diss, 
#'   edapprox = TRUE,
#'   save.stats = TRUE,
#'   stats.end = 200)
#' est$pages
#'
edgelist.meanage <- function(x, 
                             nw, 
                             el, 
                             dissolution
                             ) {
  
  if (!(missing(x))) {
    nw <- x$fit$network
    el <- x$edgelist
    dissolution <- x$dissolution
  }
  
  time.steps <- max(el$terminus)
  
  # Edges only dissolution model
  if (dissolution == ~ offset(edges)) {
    mean.page <- rep(NA, time.steps)
    for (at in 1:(time.steps-1)) {
      active <- el[(el$onset <= at & el$terminus > at) | 
                     (el$onset == at & el$terminus == at), ]
      p.ages <- at - active$onset
      mean.page[at] <- mean(p.ages, na.rm=T) 
    }
    mean.page <- mean.page[2:(length(mean.page)-1)]
    return(mean.page)
  }
  
  # Edges + nodematch dissolution model
  if (length(all.names(dissolution)) > 3 && any(all.names(dissolution) == "nodematch")) {
    
    # Decompose dissolution formula for nodematch term
    dform <- strsplit(as.character(dissolution), "[+]")[[2]]
    dform2 <- dform[grep("nodematch", dform)]
    nmatchvar <- strsplit(dform2, "[\"]")[[1]][2]
    nmatchvals <- nw %v% nmatchvar
    el$nmatch <- nmatchvals[el$head] == nmatchvals[el$tail]
    
    # Empty output vectors
    mean.page.match <- rep(NA, time.steps)
    mean.page.nomatch <- rep(NA, time.steps)
    
    for (at in 1:(time.steps-1)) {
      active.match <- el[el$nmatch == TRUE &
                           ((el$onset <= at & el$terminus > at) | 
                              (el$onset == at & el$terminus == at)), ]
      active.nomatch <- el[el$nmatch == FALSE &
                             ((el$onset <= at & el$terminus > at) | 
                                (el$onset == at & el$terminus == at)), ]
      
      p.ages.match <- at - active.match$onset
      p.ages.nomatch <- at - active.nomatch$onset
      
      mean.page.match[at] <- mean(p.ages.match, na.rm=T) 
      mean.page.nomatch[at] <- mean(p.ages.nomatch, na.rm=T)
      
    }         
    mean.page.match <- mean.page.match[2:(length(mean.page.match)-1)]
    mean.page.nomatch <- mean.page.nomatch[2:(length(mean.page.nomatch)-1)]
    return(list(mean.page.match=mean.page.match,
                mean.page.nomatch=mean.page.nomatch))
  }
}


#' @title Get Nodematch Term Proportions
#'
#' @description This function outputs a list of values and proportions of
#'   static vertex attributes present in the formation formula as nodematch
#'   terms.
#' 
#' @param nw an object of class \code{network} or \code{networkDynamic}.
#' @param formation a right-hand sided STERGM formation formula
#'  (see \code{\link{epiNet.est}}).
#' 
#' @details
#' Incoming nodes in a dynamic network simulation require values for
#' any attributes that are used in the formation formula through nodematch
#' network terms. This function gets the proportional distribution of those attribute
#' values of existing nodes in the network. The distributions may then be fed into
#' \code{set.nodematch.attributes} to update the attributes of new nodes.
#' 
#' @seealso \code{\link{set.nodematch.attributes}}
#' 
#' @export
#' @keywords epiNetUtils internal
#'
#' @examples
#' # Set network with two vertex attributes
#' nw <- network.initialize(100)
#' nw <- set.vertex.attribute(nw, "race", value = rep(0:1, each=50))
#' nw <- set.vertex.attribute(nw, "sex", 
#'                            value = c(rep("M", 10), rep("F", 90)))
#' 
#' # Formation formula uses both attributes
#' formation <- ~ edges + nodematch("race") + nodematch("sex")
#' 
#' # Get existing proportions of those two terms
#' get.nodematch.attributes(nw, formation)
#'
get.nodematch.attributes <- function(nw, 
                                     formation
                                     ) {
  
  # First find the nodematch term variables
  fterms <- (attr(terms.formula(formation), "term.labels"))
  fterms2 <- fterms[grepl("nodematch", fterms)]
  if (length(fterms2) == 0) {
    return(NULL)
  } else {
    fterms3 <- substr(fterms2, 12, nchar(fterms2))
    fterms4 <- substr(fterms3, 1, nchar(fterms3)-2)
    
    # Second match against nw vertex attributes
    vnames <- names(nw$val[[1]])[names(nw$val[[1]]) %in% fterms4]
    
    for (i in seq_along(vnames)) {
      v.vals <- get.vertex.attribute(nw, vnames[i])
      v.vals.uni <- unique(v.vals)[!is.na(unique(v.vals))]
      v.props <- as.vector(sort(prop.table(table(v.vals))))
      if (i == 1) {
        vals <- list(v.vals.uni)
        prop <- list(v.props)
      }
      if (i > 1) {
        vals[[i]] <- v.vals.uni
        prop[[i]] <- v.props
      }
    }
    v.out <- list(vnames, vals, prop)
    return(v.out)
  }
}


#' @title Get State Size from a NetworkDynamic Object
#'
#' @description This function provides all active model state sizes from 
#'   the network at the specified time step, output to a list of vectors.
#' 
#' @param all a list object containing a \code{networkDynamic} object and other
#'   initialization information passed from \code{\link{epiNet.simTrans}}.
#' @param at time step to query.
#' @param set.all if \code{TRUE}, sets output on all object, otherwise will
#'  output directly.
#' 
#' @details
#' This network utility is used during the \code{\link{epiNet.simTrans}}
#' simulation process to efficiently query the current size of each state
#' or compartment in the model at any given timestep. For a bipartite network,
#' the current state size for each mode, and overall is provided. A more flexible, 
#' but inefficient status query function is \code{get.stat}. 
#' 
#' @seealso \code{\link{get.stat}}
#' 
#' @export
#' @keywords epiNetUtils internal
#'
get.prev <- function(all, 
                     at, 
                     set.all = TRUE
                     ) {
  
  nw <- all$nw
  tea <- ifelse(any(names(nw$val[[1]]) %in% "status.active"), TRUE, FALSE)
  
  if (tea == FALSE) {
    status <- all$status
  }
  if (tea == TRUE) {
    if (class(nw)[1] != "networkDynamic")
      stop("Must supply a networkDynamic object")
    status <- get.vertex.attribute.active(nw, "status", 
                                          at=at, require.active=TRUE)
  }
  
  if (set.all == FALSE) {
    if (nw %n% "bipartite" == FALSE) {
      out <- list()
      out$s.num.m1 <- sum(status == 0, na.rm=TRUE)
      out$i.num.m1 <- sum(status == 1, na.rm=TRUE)
      if (nw$type == "SIR") {
        out$r.num.m1 <- sum(status == 2, na.rm=TRUE)
      }
    } else {
      out <- list()
      ids.m1 <- modeids(nw, 1)
      ids.m2 <- modeids(nw, 2)
      out$s.num.m1 <- sum(status[ids.m1] == 0, na.rm=TRUE)
      out$i.num.m1 <- sum(status[ids.m1] == 1, na.rm=TRUE)
      out$s.num.m2 <- sum(status[ids.m2] == 0, na.rm=TRUE)
      out$i.num.m2 <- sum(status[ids.m2] == 1, na.rm=TRUE)
      out$s.num <- out$s.num.m1 + out$s.num.m2
      out$i.num <- out$i.num.m1 + out$i.num.m2
      if (nw$type == "SIR") {
        out$r.num.m1 <- sum(status[ids.m1] == 2, na.rm=TRUE)
        out$r.num.m2 <- sum(status[ids.m2] == 2, na.rm=TRUE)
        out$r.num <- out$r.num.m1 + out$r.num.m2
      }
    }
    return(out)
  }
  if (set.all == TRUE) {
    if (nw %n% "bipartite" == FALSE) {
      if (at == 1) {
        all$s.num.m1 <- sum(status == 0, na.rm=TRUE)
        all$i.num.m1 <- sum(status == 1, na.rm=TRUE)
        if (all$type == "SIR") {
          all$r.num.m1 <- sum(status == 2, na.rm=TRUE)
        }
      } else {
        all$s.num.m1[at] <- sum(status == 0, na.rm=TRUE)
        all$i.num.m1[at] <- sum(status == 1, na.rm=TRUE)
        if (all$type == "SIR") {
          all$r.num.m1[at] <- sum(status == 2, na.rm=TRUE)
        }
      }
    } else {
      ids.m1 <- modeids(nw, 1)
      ids.m2 <- modeids(nw, 2)
      if (at == 1) {
        all$s.num.m1 <- sum(status[ids.m1] == 0, na.rm=TRUE)
        all$i.num.m1 <- sum(status[ids.m1] == 1, na.rm=TRUE)
        all$s.num.m2 <- sum(status[ids.m2] == 0, na.rm=TRUE)
        all$i.num.m2 <- sum(status[ids.m2] == 1, na.rm=TRUE)
        all$s.num <- all$s.num.m1 + all$s.num.m2
        all$i.num <- all$i.num.m1 + all$i.num.m2
        if (all$type == "SIR") {
          all$r.num.m1 <- sum(status[ids.m1] == 2, na.rm=TRUE)
          all$r.num.m2 <- sum(status[ids.m2] == 2, na.rm=TRUE)
          all$r.num <- all$r.num.m1 + all$r.num.m2
        }
      } else {
        all$s.num.m1[at] <- sum(status[ids.m1] == 0, na.rm=TRUE)
        all$i.num.m1[at] <- sum(status[ids.m1] == 1, na.rm=TRUE)
        all$s.num.m2[at] <- sum(status[ids.m2] == 0, na.rm=TRUE)
        all$i.num.m2[at] <- sum(status[ids.m2] == 1, na.rm=TRUE)
        all$s.num[at] <- all$s.num.m1[at] + all$s.num.m2[at]
        all$i.num[at] <- all$i.num.m1[at] + all$i.num.m2[at]
        if (all$type == "SIR") {
          all$r.num.m1[at] <- sum(status[ids.m1] == 2, na.rm=TRUE)
          all$r.num.m2[at] <- sum(status[ids.m2] == 2, na.rm=TRUE)
          all$r.num[at] <- all$r.num.m1[at] + all$r.num.m2[at]
        }
      }
    }
    return(all)
  }
  
}


#' @title Get Active Status from NetworkDynamic Object
#'
#' @description This function outputs information on the disease or
#'  outcome status stored in the networkDyamic status TEA used 
#'  by \code{\link{epiNet.simTrans}}.
#' 
#' @param all a list object containing a \code{networkDynamic} object and other
#'   initialization information passed from \code{\link{epiNet.simTrans}}.
#' @param stat a vector of status values in numeric to query (the standard 
#'  status representation is 0 = negative, 1 = infected, 2 = recovered).
#' @param at time step to query.
#' @param out function output, with options of \code{out="vec"} for a vector of status, 
#'  \code{out="ids"} for a vector of IDs with status equal to \code{stat}, and
#'  \code{out="prev"} for the number current with status equal to \code{stat},
#'   and \code{out="all"} to return a list of the prior three elements.
#' @param mode if \code{nw} is bipartite, the mode number for status (may be ignored
#'  if requesting output for both modes).
#' @param require.active if \code{TRUE}, NA will be returned instead of an attribute value 
#'  if the associated vertex or edge is inactive for the query period
#'  (see \code{\link{get.vertex.attribute.active}}).
#' 
#' @details 
#' This network utility is used during the \code{\link{epiNet.simTrans}}
#' simulation process to query the ID numbers and size of currently active
#' nodes with a specific disease status. A more efficient query just for
#' state sizes of status = 1 is found in \code{get.prev}.
#' 
#' @seealso \code{\link{get.prev}}
#' 
#' @export
#' @keywords epiNetUtils internal
#'
get.stat <- function(all, 
                     stat = 1, 
                     at, 
                     out, 
                     mode, 
                     require.active = TRUE
                     ) {
  
  nw <- all$nw
  if (!(missing(mode)) && !is.numeric(nw$gal$bipartite))
    stop("nw must be bipartite if mode argument is used")
  if (missing(out))
    stop("Supply an output type")
  
  tea <- ifelse(any(names(nw$val[[1]]) %in% "status.active"), TRUE, FALSE)
  
  lstat <- length(stat)
  
  if (tea == FALSE) {
    status <- all$status
    
      if (missing(mode)) {
        if (out == "vec") {
          return(status)
        } else {
          if (lstat == 1 && is.na(stat)) {
            act.stat <- is.na(status)
            if (out == "ids") {
              return(which(act.stat))
            }
            if (out == "prev") {
              return(sum(act.stat))
            }
            if (out == "all") {
              return(list(status=status,
                          ids=which(act.stat),
                          prev=sum(act.stat)))
            }
          } else {
            act.stat <- status %in% stat
            if (out == "ids") {
              return(which(act.stat))
            }
            if (out == "prev") {
              return(sum(act.stat, na.rm=T))
            }
            if (out == "all") {
              return(list(status=status,
                          ids=which(act.stat),
                          prev=sum(act.stat, na.rm=TRUE)))
            }
          } 
        }
      } else {
        ids <- modeids(nw, mode)
        if (out == "vec") {
          return(status[ids])
        } else {
          if (lstat == 1 && is.na(stat)) {
            act.stat <- is.na(status)
            if (out == "ids") {
              return(intersect(which(act.stat), ids))
            }
            if (out == "prev") {
              return(length(intersect(which(act.stat), ids)))
            }
            if (out == "all") {
              return(list(status=status[ids],
                          ids=intersect(which(act.stat), ids),
                          prev=length(intersect(which(act.stat), ids))))
            }
          } else {
            act.stat <- status %in% stat
            if (out == "ids") {
              return(intersect(which(act.stat), ids))
            }
            if (out == "prev") {
              return(length(intersect(which(act.stat), ids)))
            }
            if (out == "all") {
              return(list(status=status[ids],
                          ids=intersect(which(act.stat), ids),
                          prev=length(intersect(which(act.stat), ids))))
            }
          } 
        }
      }
    } # end tea == FALSE
    
  if (tea == TRUE) {
    if (class(nw)[1] != "networkDynamic")
      stop("Must supply a networkDynamic object")
    status <- get.vertex.attribute.active(nw, "status", 
                                          at=at, require.active=require.active)
      
    if (missing(mode)) {
      if (out == "vec") {
        return(status)
      } else {
        if (lstat == 1 && is.na(stat)) {
          act.stat <- is.na(status) & node.active(nw, at=at, out="vec")
          if (out == "ids") {
            return(which(act.stat))
          }
          if (out == "prev") {
            return(sum(act.stat))
          }
          if (out == "all") {
            return(list(status=status,
                        ids=which(act.stat),
                        prev=sum(act.stat)))
          }
        } else {
          act.stat <- status %in% stat & node.active(nw, at=at, out="vec")
          if (out == "ids") {
            return(which(act.stat))
          }
          if (out == "prev") {
            return(sum(act.stat, na.rm=T))
          }
          if (out == "all") {
            return(list(status=status,
                        ids=which(act.stat),
                        prev=sum(act.stat, na.rm=TRUE)))
          }
        } 
      }
    } else {
      ids <- modeids(nw, mode)
      if (out == "vec") {
        return(status[ids])
      } else {
        if (lstat == 1 && is.na(stat)) {
          act.stat <- is.na(status) & node.active(nw, at=at, out="vec")
          if (out == "ids") {
            return(intersect(which(act.stat), ids))
          }
          if (out == "prev") {
            return(length(intersect(which(act.stat), ids)))
          }
          if (out == "all") {
            return(list(status=status[ids],
                        ids=intersect(which(act.stat), ids),
                        prev=length(intersect(which(act.stat), ids))))
          }
        } else {
          act.stat <- status %in% stat & node.active(nw, at=at, out="vec")
          if (out == "ids") {
            return(intersect(which(act.stat), ids))
          }
          if (out == "prev") {
            return(length(intersect(which(act.stat), ids)))
          }
          if (out == "all") {
            return(list(status=status[ids],
                        ids=intersect(which(act.stat), ids),
                        prev=length(intersect(which(act.stat), ids))))
          }
        } 
      }
    }
  } # end tea == TRUE

}


#' @title Mode Numbers for Bipartite Network
#'
#' @description This function outputs mode numbers give ID numbers
#'   for a bipartite network.
#' 
#' @param nw an object of class \code{network} or \code{networkDynamic}.
#' @param ids a vector of ID numbers for which the mode number 
#'  should be returned.
#' 
#' @seealso \code{\link{modeids}} provides the reverse functionality.
#' 
#' @export
#' @keywords epiNetUtils internal
#'
#' @examples
#' nw <- network.initialize(10, bipartite = 5)
#' idmode(nw)
#' idmode(nw, ids = 3)
#'
idmode <- function(nw, 
                   ids
                   ) {
  
  if (!is.numeric(nw$gal$bipartite))
    stop("nw must be a bipartite network")
  
  if (missing(ids)) ids <- seq_len(network.size(nw))
  
  if (any(ids > network.size(nw))) 
    stop("Specify ids between 1 and ", network.size(nw))
  
  out <- 1 + as.integer(1 + ids > nw$gal$bipartite)
  
  return(out)
}


#' @title ID Numbers for Bipartite Network
#'
#' @description This function outputs ID numbers for a mode number
#'   for a bipartite network.
#' 
#' @param nw an object of class \code{network} or \code{networkDynamic}.
#' @param mode mode number to return ID numbers for.
#' 
#' @seealso \code{\link{idmode}} provides the reverse functionality.
#' 
#' @export
#' @keywords epiNetUtils internal
#' 
#' @examples
#' nw <- network.initialize(10, bipartite=5)
#' modeids(nw, mode=2)
#' 
modeids <- function(nw, 
                    mode
                    ) {
  
  if (!is.numeric(nw$gal$bipartite))
    stop("nw must be a bipartite network")
  if (missing(mode))
    stop("Specify mode=1 or mode=2")
  
  if (mode == 1) out <- 1:(nw %n% "bipartite")
  if (mode == 2) out <- ((nw %n% "bipartite")+1):(nw %n% "n")
  
  return(out)
}


#' @title Query Active Nodes in NetworkDynamic Object
#'
#' @description This function outputs information on the active nodes in
#'   a \code{networkDynamic} object.
#' 
#' @param nw an object of class \code{networkDynamic}.
#' @param at time step to query.
#' @param out function output, with options of \code{out="vec"} for 
#'  a T/F vector of whether the node is active, \code{out="ids"} for 
#'  a vector of IDs active, \code{out="prev"} for the number of 
#'  nodes that are active, and \code{out="all"} to return a list of 
#'  the prior three elements.
#' @param mode if \code{nw} is bipartite, the mode number for status (may 
#'  be ignored if requesting output for both modes).
#' @param active.default if \code{TRUE}, elements without an activity attribute 
#'  will be regarded as active.
#' 
#' @details
#' This is a specialized version of \code{\link{is.active}} from the 
#' \code{networkDynamic} package that allows for key output to be efficiently 
#' generated for use in \code{\link{epiNet.simTrans}} simulations. For functions 
#' that query active nodes with specific values of the disease status TEA, 
#' use \code{\link{get.stat}} or \code{\link{get.prev}}.
#' 
#' @seealso \code{\link{is.active}}, \code{\link{get.stat}}, \code{\link{get.prev}}
#' 
#' @export
#' @keywords epiNetUtils internal
#' 
#' @examples
#' # Initialize NW and activate vertices
#' nw <- network.initialize(20)
#' 
#' # Activate all vertices, then deactive half at time 5
#' activate.vertices(nw, onset = 1, terminus = 10)
#' deactivate.vertices(nw, onset = 5, terminus = 10, v = 1:10)
#' 
#' # Output all information for vertices at time 1 and time 5
#' node.active(nw, at = 1, out = "all")
#' node.active(nw, at = 5, out = "all")
#'
node.active <- function(nw, 
                        at, 
                        out, 
                        mode, 
                        active.default = FALSE
                        ) {
  
  if (!(missing(mode)) && !is.numeric(nw$gal$bipartite))
    stop("nw must be bipartite if mode argument is used")
  
  if (out %in% c("vec", "ids", "prev")) {
    if (missing(mode)) {
      node.active <- is.active(nw, v=seq_len(network.size(nw)), at=at,
                               active.default=active.default)
      out.vec <- node.active
      out.ids <- which(node.active)
      out.prev <- sum(node.active)
    } else {
      node.active <- is.active(nw, v=seq_len(network.size(nw)), at=at,
                               active.default=active.default)
      ids.m1 <- modeids(nw, 1)
      ids.m2 <- modeids(nw, 2)
      if (mode == 1) {
        out.vec <- node.active[ids.m1]
        out.ids <- intersect(which(node.active), ids.m1)
        out.prev <- sum(node.active[ids.m1])
      }
      if (mode == 2) {
        out.vec <- node.active[ids.m2]
        out.ids <- intersect(which(node.active), ids.m2)
        out.prev <- sum(node.active[ids.m2])
      }
    }
  }
  if (out == "all") {
    if (!is.numeric(nw$gal$bipartite)) {
      node.active <- is.active(nw, v=seq_len(network.size(nw)), at=at,
                               active.default=active.default)
      out.all <- list()
      out.all$vec$all <- node.active
      out.all$ids$all <- which(node.active)
      out.all$prev$all <- sum(node.active)
    } else {
      node.active <- is.active(nw, v=seq_len(network.size(nw)), at=at,
                               active.default=active.default)
      out.all <- list()
      ids.m1 <- modeids(nw, 1)
      ids.m2 <- modeids(nw, 2)
      out.all$vec$m1 <- node.active[ids.m1]
      out.all$vec$m2 <- node.active[ids.m2]
      out.all$vec$all <- node.active
      out.all$ids$m1 <- intersect(which(node.active), ids.m1)
      out.all$ids$m2 <- intersect(which(node.active), ids.m2)
      out.all$ids$all <- which(node.active)
      out.all$prev$m1 <- sum(node.active[ids.m1])
      out.all$prev$m2 <- sum(node.active[ids.m2])
      out.all$prev$all <- sum(node.active)
    }
  }
  
  if (out == "vec") return(out.vec)
  if (out == "ids") return(out.ids)
  if (out == "prev") return(out.prev)
  if (out == "all") return(out.all)
  
}


#' @title Set Nodematch Values from Proportions
#'
#' @description This function sets the values of static vertex attributes present 
#'   in the formation formula as nodematch terms.
#' 
#' @param nw an object of class \code{network} or \code{networkDynamic}.
#' @param nm.attr object output from \code{\link{get.nodematch.attributes}}.
#' @param new.nodes vector of new nodes with missing (that is, \code{NA}) 
#'   attribute values.
#'   
#' @details
#' Incoming nodes in a dynamic network simulation require values for
#' any attributes that are used in the formation formula through nodematch
#' network terms. This function sets the values of those attributes based on
#' the proportional distribution of the existing node attributes. The distribution
#' must first be obtained from \code{get.nodematch.attributes}.
#' 
#' @seealso \code{\link{get.nodematch.attributes}}
#' 
#' @export
#' @keywords epiNetUtils internal
#' 
#' @examples
#' # Set network with two vertex attributes and check distribution
#' nw <- network.initialize(100)
#' nw %v% "race" <- rep(0:1, each=50)
#' nw %v% "sex" <- c(rep("M", 10), rep("F", 90))
#' prop.table(table(nw %v% "race"))
#' prop.table(table(nw %v% "sex"))
#' 
#' # Formation formula uses both attributes
#' formation <- ~ edges + nodematch("race") + nodematch("sex")
#' 
#' # Need to get existing proportions of those two terms first
#' nm.attr <- get.nodematch.attributes(nw, formation)
#' 
#' # Now add 50 more nodes, and set the new attributes in proportion
#' add.vertices(nw, nv=50)
#' nw <- set.nodematch.attributes(nw, nm.attr, new.nodes=101:150)
#' 
#' # Some random variation due to sampling new node values
#' prop.table(table(nw %v% "race"))
#' prop.table(table(nw %v% "sex"))
#'
set.nodematch.attributes <- function(nw, 
                                     nm.attr, 
                                     new.nodes
                                     ) {
  
  nma <- nm.attr
  for (i in 1:length(nma[[1]])) {
    samp.vals <- nma[[2]][[i]] 
    samp.probs <- nma[[3]][[i]]
    samp.n <- length(new.nodes)
    samp.out <- ssample(samp.vals, samp.n, TRUE, samp.probs)
    set.vertex.attribute(nw, nma[[1]][i], value=samp.out, v=new.nodes)
  }
  invisible(nw)
}



# Non-exported functions --------------------------------------------------

## Runs the console progress tracker for epiNet.simTrans when verbose = TRUE
cons.prog <- function(all, 
                      s, 
                      ts
                      ) {

    cat("\014")
    cat("\nepiNet Disease Simulation")
    cat("\n----------------------------")
    cat("\nSimulation: ", s, "/", all$tot.sims, sep="")
    cat("\nTimestep: ", ts, "/", all$nsteps, sep="")
    if (all$modes == 1) cat("\nIncidence:", all$si.flow.m1[ts])
    if (all$modes == 2) cat("\nIncidence:", all$si.flow.m1[ts]+all$si.flow.m2[ts])
    if (all$type %in% c("SIR", "SIS")) {
      if (all$modes == 1) cat("\nRecoveries:", all$recovs.m1[ts])
      if (all$modes == 2) cat("\nRecoveries:", all$recovs.m1[ts]+all$recovs.m2[ts])
    }
    if (all$modes == 1) cat("\nPrevalence:", all$i.num.m1[ts])
    if (all$modes == 2) cat("\nPrevalence:", all$i.num.m1[ts]+all$i.num.m2[ts])
    if (all$type %in% c("SI", "SIS")) {
      if (all$modes == 1) cat("\nPopulation:", all$s.num.m1[ts]+all$i.num.m1[ts])
      if (all$modes == 2) cat("\nPopulation:", all$s.num.m1[ts]+all$s.num.m2[ts]+
                                              all$i.num.m1[ts]+all$i.num.m2[ts])
    }
    if (all$type == "SIR") {
      if (all$modes == 1) {
        cat("\nPopulation:", all$s.num.m1[ts]+all$i.num.m1[ts]+all$r.num.m1[ts])
      }
      if (all$modes == 2) {
        cat("\nPopulation:", all$s.num.m1[ts]+all$i.num.m1[ts]+all$r.num.m1[ts]+
                             all$s.num.m2[ts]+all$i.num.m2[ts]+all$r.num.m2[ts])
      }
    }
    if (all$vital == TRUE) {
      if (all$modes == 1) {
        cat("\nBirths:", all$b.flow.m1[ts])
        cat("\nDeaths, susceptibles:", all$ds.flow.m1[ts])
        cat("\nDeaths, infecteds:", all$di.flow.m1[ts]) 
        if (all$type == "SIR") {
          cat("\nDeaths, recovered:", all$dr.flow.m1[ts]) 
        }
      }
      if (all$modes == 2) {
        cat("\nBirths:", all$b.flow.m1[ts]+all$b.flow.m2[ts])
        cat("\nDeaths, susceptible:", all$ds.flow.m1[ts]+all$ds.flow.m2[ts])
        cat("\nDeaths, infected:", all$di.flow.m1[ts]+all$di.flow.m2[ts]) 
        if (all$type == "SIR") {
          cat("\nDeaths, recovered:", all$dr.flow.m1[ts]+all$dr.flow.m2[ts]) 
        }
      }
    }
    cat("\n----------------------------")

}


## Population size correction for epiNet.simTrans with vital dynamics
edges.correct <- function(all, 
                          ts, 
                          coef.form
) {
  
  if (all$modes == 1) {
    if (all$type %in% c("SI", "SIS")) {
      old.num <- all$s.num.m1[ts-1] + all$i.num.m1[ts-1]
      new.num <- all$s.num.m1[ts] + all$i.num.m1[ts]
    }
    if (all$type == "SIR") {
      old.num <- all$s.num.m1[ts-1] + all$i.num.m1[ts-1] + all$r.num.m1[ts-1]
      new.num <- all$s.num.m1[ts] + all$i.num.m1[ts] + all$r.num.m1[ts]
    }
    coef.form[1] <- coef.form[1] + log(old.num) - log(new.num)
  }
  if (all$modes == 2) {
    if (all$type %in% c("SI", "SIS")) {
      old.num.m1 <- all$s.num.m1[ts-1] + all$i.num.m1[ts-1]
      old.num.m2 <- all$s.num.m2[ts-1] + all$i.num.m2[ts-1]
      new.num.m1 <- all$s.num.m1[ts] + all$i.num.m1[ts]
      new.num.m2 <- all$s.num.m2[ts] + all$i.num.m2[ts]
    }
    if (all$type == "SIR") {
      old.num.m1 <- all$s.num.m1[ts-1] + all$i.num.m1[ts-1] + all$r.num.m1[ts-1]
      old.num.m2 <- all$s.num.m2[ts-1] + all$i.num.m2[ts-1] + all$r.num.m2[ts-1]
      new.num.m1 <- all$s.num.m1[ts] + all$i.num.m1[ts] + all$r.num.m1[ts]
      new.num.m2 <- all$s.num.m2[ts] + all$i.num.m2[ts] + all$r.num.m2[ts]
    }
    coef.form[1] <- coef.form[1] + log(2*old.num.m1*old.num.m2/(old.num.m1+old.num.m2)) - 
      log(2*new.num.m1*new.num.m2/(new.num.m1+new.num.m2))
  }
  
  return(coef.form)
}

# logit transformation of a probability
logit <- function(x) {
  log(x / (1 - x))
}

# Runs the progress plot for epiNet.simTrans when plot.prog = TRUE
prog.plot <- function(all, 
                      ts
                      ) {
  
  par(mfrow=c(1,2), mar=c(3,3,2,1), mgp=c(2,1,0), ask=FALSE)
  if (all$modes == 1) {
    y.prev <- all$i.num.m1/(all$s.num.m1+all$i.num.m1)
    if (all$type == "SIR") {
      y.prev <- all$i.num.m1/(all$s.num.m1+all$i.num.m1+all$r.num.m1)
    }
    plot(1:ts, y.prev, 
         type="l", lwd=2, col="firebrick",
         xlim=c(0, all$nsteps), ylim=0:1,
         xlab="Time", ylab="Prevalence", 
         main=paste("Disease Prevalence at Time", ts))
    plot(1:ts, all$si.flow.m1, type="h", xlim=c(0, all$nsteps),
         xlab="Time", ylab="Incidence", 
         main=paste("Disease Incidence at Time", ts))
  }
  if (all$modes == 2) {
    y.prev.m1 <- all$i.num.m1/(all$s.num.m1+all$i.num.m1)
    if (all$type == "SIR") {
      y.prev.m1 <- all$i.num.m1/(all$s.num.m1+all$i.num.m1+all$r.num.m1)
    }
    y.prev.m2 <- all$i.num.m2/(all$s.num.m2+all$i.num.m2)
    if (all$type == "SIR") {
      y.prev.m2 <- all$i.num.m2/(all$s.num.m2+all$i.num.m2+all$r.num.m2)
    }
    y.inci <- all$si.flow.m1+all$si.flow.m2
    plot(1:ts, y.prev.m1, type="l", 
         lwd=2, col="firebrick",
         xlim=c(0, all$nsteps), ylim=0:1,
         xlab="Time", ylab="Prevalence", 
         main=paste("Disease Prevalence at Time", ts))
    lines(1:ts, y.prev.m2, lwd=2, col="steelblue")
    legend("topleft", legend=c("m1", "m2"), lty=1, lwd=2, 
           col=c("firebrick", "steelblue"), cex=0.8)
    plot(1:ts, y.inci, type="h", xlim=c(0, all$nsteps),
         xlab="Time", ylab="Incidence", 
         main=paste("Total Disease Incidence at Time", ts))
  }
  
}

# When running batch mode, this function reads the end of the 
# potentially huge Rout to see verbose output
read.end.Rout <- function(Rout, 
                          last.lines = 10
                          ) {
  
  prog <- readLines(Rout, warn=F)
  out <- cbind(prog[(length(prog)-last.lines):length(prog)])
  rownames(out) <- rep('', nrow(out))
  colnames(out) <- ''
  print(out)
  
}

# Saves out object x to unique file with date/time stamp
sim.autosave <- function(x) {
  
  fn <- paste('sim', gsub(":", "-", 
                          gsub(" ", "_", as.character(Sys.time()))), 
              'Rdata', sep='.')
  
  save(x, file=fn)
}
