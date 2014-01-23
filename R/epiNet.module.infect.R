
##
## Infection modules for epiNet.simTrans
##

#' @title Primary Infection Module for epiNet.simTrans
#'
#' @description This function simulates the main infection process given the 
#'  current state of the partnerships and disease in the system.
#' 
#' @param all a list object containing a \code{networkDynamic} object and other
#'   initialization information passed from \code{\link{epiNet.simTrans}}.
#' @param at current time step.
#' 
#' @details
#' The main steps in this infection module are as follows: 
#' \enumerate{
#'  \item Get IDs for current infected and susceptibles given the current disease
#'    status.
#'  \item Call \code{\link{discord.edgelist}} to get the current discordant edgelist 
#'    given step 1.
#'  \item Determine the transmission rates (e.g., as a function of mode).
#'  \item Pull the number of acts per partnership in a time step from the 
#'    \code{act.rate} parameter.
#'  \item Calculate the final transmission probabilities given the transmission 
#'    rates and act rates.
#'  \item Randomly transmit on the discordant edgelist.
#'  \item Conduct bookkeeping for new infections to update status on the nodes
#'    and calculate disease incidence.
#' }
#' 
#' @return
#' The main \code{all} object is returned with updated disease status and summary
#' incidence measures.
#' 
#' @export
#' @keywords epiNetModule internal
#' 
#' @seealso \code{\link{discord.edgelist}} is used within \code{infection} to obtain
#'  discordant edgelist.
#'
infection <- function(all, at) {
  
    trans.rate <- all$trans.rate
    if (is.null(all$trans.rate.m2)) trans.rate.m2 <- NULL
    trans.rate.m2 <- all$trans.rate.m2
    act.rate <- all$act.rate
    if (is.null(all$act.rate)) act.rate <- 1
    
    # Reset new infections vector
    new.inf <- NULL
    all$status.old <- all$status
  
    nw <- all$nw
    tea <- ifelse(any(names(nw$val[[1]]) %in% 'status.active'), TRUE, FALSE)
    
    # Vector of infected and susceptible IDs at last step
    infecteds <- get.stat(all, stat=1, at=at, out='ids')
    susceptible <- get.stat(all, stat=0, at=at, out='ids')
    
    # If some infected AND some susceptible, then proceed
    if (length(infecteds) > 0 && length(infecteds) < node.active(nw, at=at, out='prev')) {
      
      # get discordant edgelist
      df <- discord.edgelist(all, infecteds, susceptible, at=at)
      
      # If some discordant edges, then proceed
      if (!(is.null(df))) {
        
        inf.time <- all$inf.time
        df$inft <- at - inf.time[df$inf]
        df$inft[df$inft == 0] <- 1
        
        # Calculate infection-stage transmission rates
        ltrans.rate <- length(trans.rate)
        if (is.null(trans.rate.m2)) {
          df$trans.rates <- ifelse(df$inft <= ltrans.rate, trans.rate[df$inft], 
                                   trans.rate[ltrans.rate])
        } else {
          df$trans.rates <- ifelse(df$sus <= nw %n% 'bipartite', 
                                   ifelse(df$inft <= ltrans.rate, trans.rate[df$inft], 
                                          trans.rate[ltrans.rate]),
                                   ifelse(df$inft <= ltrans.rate, trans.rate.m2[df$inft], 
                                          trans.rate.m2[ltrans.rate]))
        }
        
        # Calculate infection-stage act/contact rates
        lact.rate <- length(act.rate)
        df$act.rates <- ifelse(df$inft <= lact.rate, act.rate[df$inft], 
                                 act.rate[lact.rate])
        
        # Calculate final transmission probability per timestep
        df$tprob <- 1-(1-df$trans.rates)^df$act.rates
        
        # Randomize transmissions and subset df
        transmit <- runif(nrow(df)) <= df$tprob
        df <- df[transmit, ]
        
        # Set new infections vector
        new.inf <- unique(df$sus)
        if (length(df$sus) == 0) new.inf <- NULL
        
        # Update nw attributes
        if (!is.null(new.inf)) {
          if (tea == TRUE) {
            nw <- set.vertex.attribute(nw, 'inf.time', value=at, v=new.inf)
            nw <- set.vertex.attribute(nw, 'transmitID', value=df$inf, v=new.inf)
            nw <- activate.vertex.attribute(nw, 'status', value=1, onset=at, 
                                      terminus=Inf, v=new.inf)
          }
          all$status[new.inf] <- 1
          all$inf.time[new.inf] <- at
        }
        
        # Substitute PIDs for vital bipartite sims
        if (any(names(nw$gal) %in% 'vertex.pid')) {
          df$sus <- get.vertex.pid(nw, df$sus)
          df$inf <- get.vertex.pid(nw, df$inf)
        }
        
      } # end some discordant edges condition
    } # end some active discordant nodes condition
  
    if (all$save.statmat == TRUE) {
      all$stat.mat[at, ] <- all$status
    }
    
    if (is.null(new.inf)) df <- NULL
    if (at == 2) {
      all$transdf <- df
    } else {
      all$transdf <- rbind(all$transdf, df)
    }
    
    ## Save incidence vector
    if (nw %n% 'bipartite' == FALSE) {
      if (at == 2) {
        all$si.flow.m1 <- c(0, length(new.inf))
      } else {
        all$si.flow.m1[at] <- length(new.inf)
      }
    } else{
      bip.cutoff <- nw %n% 'bipartite'
      if (at == 2) {
        all$si.flow.m1 <- c(0, length(new.inf[new.inf <= bip.cutoff]))
        all$si.flow.m2 <- c(0, length(new.inf[new.inf > bip.cutoff]))
      } else {
        all$si.flow.m1[at] <- length(new.inf[new.inf <= bip.cutoff])
        all$si.flow.m2[at] <- length(new.inf[new.inf > bip.cutoff])
      }
    }
  
   all$nw <- nw
   return(all)
}


#' @title Discordant Edgelist from NetworkDynamic Object
#'
#' @description This function returns a \code{data.frame} with a discordant edgelist, 
#'  defined as the set of edges in which the status of the two partners is one
#'  susceptible and one infected.
#' 
#' @param all a list object containing a \code{networkDynamic} object and other
#'   initialization information passed from \code{\link{epiNet.simTrans}}.
#' @param infecteds vector of IDs for currently infecteds.
#' @param susceptible vector of IDs for currently susceptible.
#' @param at current time step.
#' 
#' @details
#' This internal function works within the parent \code{\link{infection}} function
#' to pull the current edgelist from the dynamic network object, look up the disease
#' status of the head and tails on the edge, and subset the list to those edges
#' with one susceptible and one infected node.
#' 
#' @return
#' This function returns a \code{data.frame} with the following columns:
#' \itemize{
#'  \item \strong{time:} time step queried
#'  \item \strong{inf:} ID number for the infected partner
#'  \item \strong{sus:} ID number for the susceptible partner
#'  \item \strong{eid:} Edge ID number for that unique dyad
#'  \item \strong{infdeg:} Current degree for the infected partner
#'  \item \strong{susdeg:} Current degree for the susceptible partner
#' }
#' The output from this function is added to the transmission \code{data.frame}
#' object that is requested as output in \code{epiNet.simTrans} simulations with
#' the \code{save.trans=TRUE} argument. 
#' 
#' @seealso \code{\link{epiNet.simTrans}}, \code{\link{infection}}
#' 
#' @export
#' @keywords epiNetModule internal
#'
discord.edgelist <- function(all, infecteds, susceptible, at) {
  
  if (length(infecteds) == 0) return(NULL)
  
  # ties that are active at time step
  active <- get.dyads.active(all$nw, at=at)
  active.order <- sample(1:nrow(active))
  ties <- as.data.frame(active[active.order, ])
  colnames(ties) <- c('tail', 'head')
  
  # no ties
  if (nrow(ties) == 0) {
    return(NULL)
  } else {
    # discordant ties
    disc.ties <- ties[sapply(1:nrow(ties), 
                             function(i) xor(ties[i, 'head'] %in% infecteds && 
                                               ties[i, 'tail'] %in% susceptible, 
                                             ties[i, 'tail'] %in% infecteds &&
                                               ties[i, 'head'] %in% susceptible)), ]    
    
    # no discordant ties
    if (nrow(disc.ties) == 0) {
      return(NULL)
    } else {
      inf <- disc.ties$head
      sus <- disc.ties$tail
      
      for (i in seq_along(sus)) {
        if (sus[i] %in% infecteds) {
          sus[i] <- disc.ties[i, 'head']
          inf[i] <- disc.ties[i, 'tail']
        }
      }
      time <- rep(at, length(sus))
      
      infdeg <- sapply(inf, function(k) sum(ties$head == k) + sum(ties$tail == k))
      susdeg <- sapply(sus, function(k) sum(ties$head == k) + sum(ties$tail == k))
      
      de.df <- data.frame(time, 
                          inf, 
                          sus, 
                          infdeg, 
                          susdeg)
      
      return(de.df)
      
    }
  }
}