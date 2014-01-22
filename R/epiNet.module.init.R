##
## EpiModel Initialization Modules for epiNet.simTrans:
##   These modules govern the initial disease status and
##   conduct bookkeeping on the network structure at t0
##   before the main disease simulation loop is started.
##
##

#' @title Disease Status Initialization Module for epiNet.simTrans
#'
#' @description This function sets the initial disease status on the 
#'   network given the supplied prevalence.
#' 
#' @param all a list object containing a \code{networkDynamic} object and other
#'   initialization information passed from \code{\link{epiNet.simTrans}}.
#' @param tea if \code{TRUE}, set status in a temporally extended attribute, 
#'  otherwise use a static vertex attribute.
#' 
#' @details
#' This internal function sets, either randomly or deterministically, the nodes
#' that are infected at the starting time of \code{epiNet} simulations, \eqn{t_1}.
#' If the number to be initially infected is passed, this function may set the
#' initial number infected based on the number specified, either as a a set of 
#' random draws from a binomial distribution or as the exact number specified. In
#' either case, the specific nodes infected are a random sample from the network.
#' In contrast, a set of specific nodes may be infected by passing the vector to
#' \code{\link{epiNet.simTrans}}.
#' 
#' @seealso This is an initialization module for \code{\link{epiNet.simTrans}}.
#'  It precedes the infection time module in \code{\link{init.inf.time}}.
#' 
#' @export
#' @keywords epiNetModule internal
#'
init.status <- function(all, tea = FALSE) {

  i.num <- all$i.num.t0
  i.ids <- all$i.ids.t0
  i.rand <- all$i.rand
  
  ## Errors
  if (!is.null(i.num) & !is.null(i.ids)) {
    stop('Use i.num OR i.ids to set initial infected')
  }
  if (is.null(i.num) & is.null(i.ids)) {
    stop('Use i.num or i.ids to set initial infected')
  }
  
  ## Enforce deterministic if passing node IDs
  if (!is.null(i.ids)) i.rand <- FALSE
  
  num <- network.size(all$nw)
  
  ## Random prevalence
  if (i.rand == TRUE) {
    if (!is.null(i.num)) {
      if (i.num < 1) stop('i.num must be >= 1')
      status <- rbinom(num, 1, i.num/num)
    }
    # If none infected, then set at least 1 infected
    if (sum(status) == 0) {
      status[ssample(1:length(status), 1)] <- 1
    }
  }
  
  ## Deterministic prevalence
  if (i.rand == FALSE) {
    status <- rep(0, num)
    if (!is.null(i.num)) {
      if (i.num < 1) 
        stop('i.num must be >= 1')
      status[ssample(1:length(status), i.num)] <- 1
    }
    if (!is.null(i.ids)) {
      if (min(i.ids) < 1 || max(i.ids) > num)
        stop(paste('Set i.ids in range of 1 to', max(i.ids)))
      status[i.ids] <- 1
    }
  }

  all$status <- all$status.old <- status
  if (tea == TRUE) {
    all$nw <- activate.vertex.attribute(all$nw, prefix='status', value=0, 
                              onset=1, terminus=Inf)
    all$nw <- activate.vertex.attribute(all$nw, prefix='status', value=1, 
                              onset=1, terminus=Inf, v=which(status == 1))
  } else {
    if (all$save.statmat == TRUE) {
      all$stat.mat <- matrix(NA, ncol=num, nrow=all$nsteps)
      all$stat.mat[1, ] <- status
    }
  }
  
  ## Save summary output to all for timestep 1
  all <- get.prev(all, at = 1, set.all = TRUE)
  
  return(all)
}


#' @title Infection Time Module for epiNet.simTrans
#'
#' @description This function sets the initial time of infection
#'   for those with disease status infected.
#' 
#' @param all a list object containing a \code{networkDynamic} object and other
#'   initialization information passed from \code{\link{epiNet.simTrans}}.
#' 
#' @details
#' This module sets the time of infection for those nodes set infected
#' at the starting time of \code{epiNet} simulations, \eqn{t_1}. For vital 
#' dynamics models, the infection time for those nodes is a random draw from an 
#' exponential distribution with the rate parameter defined by the \code{di.rate} 
#' argument. For models without vital dynamics, the infection time is a random 
#' draw from a uniform distribution of integers with a minimum of 1 and a maximum 
#' of the number of time steps in the model. In both cases, to set the infection 
#' times to be in the past, these times are multiplied by -1, and 2 is added to 
#' allow for possible infection times up until step 2, when the disease simulation 
#' time loop starts.
#' 
#' @seealso This is an initialization module for \code{\link{epiNet.simTrans}}.
#'  It follows the infection status module in \code{\link{init.status}}.
#' 
#' @export
#' @keywords epiNetModule internal
#'
init.inf.time <- function(all) {
  
  nw <- all$nw
  
  ## Set up inf.time vector
  infecteds <- get.stat(all, stat=1, at=1, out='ids')
  inf.time <- rep(NA, network.size(nw))
  
  # If vital=TRUE, inf.time is a uniform draw over the duration of infection
  if (all$vital == TRUE) {
    if (all$di.rate <= 0)
      stop('Specify di.rate > 0 if vital=TRUE')
    inf.time[infecteds] <- round(-rexp(n=length(infecteds), rate=all$di.rate)+2)
  } 
  
  # If vital=FALSE, inf.time a uniform draw over the number of sim time steps
  if (all$vital == FALSE)  {
    inf.time[infecteds] <- ssample(1:(-all$nsteps+2), length(infecteds), replace=TRUE)
  }
  
  all$inf.time <- inf.time
  return(all)
}


#' @title Persistent ID Initialization
#'
#' @description This function initializes the persistent IDs for
#'   a \code{networkDynamic} object.
#' 
#' @param nw an object of class \code{networkDynamic}.
#' @param prefixes character string prefix for mode-specific ID, 
#'   with default to c('F', 'M').
#' 
#' @details
#' This function is used for \code{\link{epiNet.simTrans}} simulations
#' over bipartite networks for populations with vital dynamics. Persistent IDs are
#' required in this situation because when new nodes are added to the 
#' first mode in a bipartite network, the IDs for the second mode shift
#' upward. Persistent IDs allow for an analysis of disease transmission
#' chains for these simulations. These IDs are also invoked in the 
#' \code{\link{births}} module when the persistent IDs of incoming nodes
#' must be set.  
#' 
#' @export
#' @keywords epiNetModule internal
#' @seealso \code{\link{initialize.pids}}
#'
#' @examples
#' # Initialize network with 25 female and 75 male
#' nw <- network.initialize(100, bipartite=25)
#' 
#' # Set persistent IDs using the default F/M prefix
#' nw <- init.pids(nw)
#' nw %v% 'vertex.names'
#' 
#' # Use another prefix combination
#' nw <- init.pids(nw, c('A', 'B'))
#' nw %v% 'vertex.names'
#'
init.pids <- function(nw, prefixes=c('F', 'M')) {
  
  # Set persistent IDs
  t0.pids <- c(paste(prefixes[1], 1:length(modeids(nw, 1)), sep=''),
               paste(prefixes[2], 1:length(modeids(nw, 2)), sep='')) 
  
  # Initialize persistent IDs on network
  nw <- set.network.attribute(nw, 'vertex.pid', 'vertex.names')
  nw <- set.vertex.attribute(nw, 'vertex.names', t0.pids)
  
  return(nw)
}
