
##
## Vital dynamics modules for epiNet.simTrans
##

#' @title Deaths for Suceptibles: epiNet.simTrans Module
#'
#' @description This function simulates death for those in the susceptible state
#'   for use in \code{\link{epiNet.simTrans}} simulations.
#' 
#' @param all a list object containing a \code{networkDynamic} object and other
#'   initialization information passed from \code{\link{epiNet.simTrans}}.
#' @param at time step for simulation.
#' 
#' @seealso \code{\link{epiNet.simTrans}}, \code{\link{deaths.inf}}, 
#'  \code{\link{deaths.rec}}
#' 
#' @export
#' @keywords epiNetModule internal
#' 
#'
deaths.sus <- function(all, at) {
  
  nw <- all$nw
  ds.rate <- all$ds.rate
  tea <- ifelse(any(names(nw$val[[1]]) %in% 'status.active'), TRUE, FALSE)
  
  sus.temp <- get.stat(all, stat=0, at=at, out='all')
  if (sus.temp$prev > 0) {
    new.sus.deaths <- rpois(1, ds.rate * sus.temp$prev)
    d.sus <- ssample(sus.temp$ids, new.sus.deaths)
    if (new.sus.deaths > 0) {
      nw <- deactivate.vertices(nw, onset=at, terminus=Inf, v=d.sus,
                                deactivate.edges=TRUE)
    }
  } else {
    d.sus <- NULL
  }
  
  all$status[d.sus] <- 9
  if (all$save.statmat == TRUE) {
    all$stat.mat[at, ] <- all$status
  }
  
  if (nw %n% 'bipartite' == FALSE) {
    if (at == 2) {
      all$ds.flow.m1 <- c(0, length(d.sus))
    } else {
      all$ds.flow.m1[at] <- length(d.sus)
    }
  } else{
    bip.cutoff <- nw %n% 'bipartite'
    if (at == 2) {
      all$ds.flow.m1 <- c(0, length(d.sus[d.sus <= bip.cutoff]))
      all$ds.flow.m2 <- c(0, length(d.sus[d.sus > bip.cutoff]))
    } else {
      all$ds.flow.m1[at] <- length(d.sus[d.sus <= bip.cutoff])
      all$ds.flow.m2[at] <- length(d.sus[d.sus > bip.cutoff])
    }
  }
  
  all$nw <- nw
  return(all)
}


#' @title Deaths for Infecteds: epiNet.simTrans Module
#'
#' @description This function simulates deaths for infecteds
#'   for use in \code{\link{epiNet.simTrans}} simulations.
#' 
#' @param all a list object containing a \code{networkDynamic} object and other
#'   initialization information passed from \code{\link{epiNet.simTrans}}.
#' @param at time step for simulation.
#' 
#' @seealso \code{\link{epiNet.simTrans}}, \code{\link{deaths.sus}}, 
#'  \code{\link{deaths.rec}}
#' 
#' @export
#' @keywords epiNetModule internal
#' 
#'
deaths.inf <- function(all, at) {
  
  nw <- all$nw
  tea <- ifelse(any(names(nw$val[[1]]) %in% 'status.active'), TRUE, FALSE)
  
  inf.temp <- get.stat(all, stat=1, at=at, out='all')
  if (inf.temp$prev > 0) {
      
    di.rate <- all$di.rate
    df <- data.frame(infecteds = inf.temp$ids)
    inf.time <- all$inf.time
    df$inft <- at - inf.time[df$infecteds]
    df$inft[df$inft == 0] <- 1
    ldi.rate <- length(di.rate)
    df$di.rates <- ifelse(df$inft <= ldi.rate, di.rate[df$inft], 
                           di.rate[ldi.rate])
    df$new.inf.deaths <- runif(nrow(df)) <= df$di.rates
    d.inf <- df$infecteds[df$new.inf.deaths]
    
    if (length(d.inf) > 0) {
      nw <- deactivate.vertices(nw, onset=at, terminus=Inf, v=d.inf,
                                deactivate.edges=TRUE)
    }
  } else {
    d.inf <- NULL
  }
  
  all$status[d.inf] <- 9
  if (all$save.statmat == TRUE) {
    all$stat.mat[at, ] <- all$status
  }
  
  if (nw %n% 'bipartite' == FALSE) {
    if (at == 2) {
      all$di.flow.m1 <- c(0, length(d.inf))
    } else {
      all$di.flow.m1[at] <- length(d.inf)
    }
  } else{
    bip.cutoff <- nw %n% 'bipartite'
    if (at == 2) {
      all$di.flow.m1 <- c(0, length(d.inf[d.inf <= bip.cutoff]))
      all$di.flow.m2 <- c(0, length(d.inf[d.inf > bip.cutoff]))
    } else {
      all$di.flow.m1[at] <- length(d.inf[d.inf <= bip.cutoff])
      all$di.flow.m2[at] <- length(d.inf[d.inf > bip.cutoff])
    }
  }
  
  all$nw <- nw
  return(all)
}


#' @title Deaths for Recovered: epiNet.simTrans Module
#'
#' @description This function simulates death for those in the susceptible state
#'   for use in \code{\link{epiNet.simTrans}} simulations.
#' 
#' @param all a list object containing a \code{networkDynamic} object and other
#'   initialization information passed from \code{\link{epiNet.simTrans}}.
#' @param at time step for simulation.
#' 
#' @seealso \code{\link{epiNet.simTrans}}, \code{\link{deaths.sus}}, 
#'  \code{\link{deaths.inf}}
#' 
#' @export
#' @keywords epiNetModule internal
#'
deaths.rec <- function(all, at) {
  
  if (all$type == 'SIR') {
  nw <- all$nw
  tea <- ifelse(any(names(nw$val[[1]]) %in% 'status.active'), TRUE, FALSE)
  
    dr.rate <- all$dr.rate
    rec.temp <- get.stat(all, stat=2, at=at, out='all')
    if (rec.temp$prev > 0) {
      new.rec.deaths <- rpois(1, dr.rate * rec.temp$prev)
      d.rec <- ssample(rec.temp$ids, new.rec.deaths)
      if (new.rec.deaths > 0) {
        nw <- deactivate.vertices(nw, onset=at, terminus=Inf, v=d.rec,
                                  deactivate.edges=TRUE)
      }
    } else {
      d.rec <- NULL
    }
    
    all$status[d.rec] <- 9
    if (all$save.statmat == TRUE) {
      all$stat.mat[at, ] <- all$status
    }
    
    if (nw %n% 'bipartite' == FALSE) {
      if (at == 2) {
        all$dr.flow.m1 <- c(0, length(d.rec))
      } else {
        all$dr.flow.m1[at] <- length(d.rec)
      }
    } else {
      bip.cutoff <- nw %n% 'bipartite'
      if (at == 2) {
        all$dr.flow.m1 <- c(0, length(d.rec[d.rec <= bip.cutoff]))
        all$dr.flow.m2 <- c(0, length(d.rec[d.rec > bip.cutoff]))
      } else {
        all$dr.flow.m1[at] <- length(d.rec[d.rec <= bip.cutoff])
        all$dr.flow.m2[at] <- length(d.rec[d.rec > bip.cutoff])
      }
    }
  all$nw <- nw
  }
  return(all)
}


#' @title Births: epiNet.simTrans Module
#'
#' @description This function simulates new births into the network
#'   for use in \code{\link{epiNet.simTrans}} simulations.
#' 
#' @param all a list object containing a \code{networkDynamic} object and other
#'   initialization information passed from \code{\link{epiNet.simTrans}}.
#' @param at time step for simulation.
#' 
#' @seealso \code{\link{epiNet.simTrans}}
#' 
#' @export
#' @keywords epiNetModule internal
#'
births <- function(all, at) {
  
  b.rate <- all$b.rate
  formation <- all$formation
  nw <- all$nw
  modes <- all$modes
  tea <- ifelse(any(names(nw$val[[1]]) %in% 'status.active'), TRUE, FALSE)
  
  node.stats <- node.active(nw, at=at-1, out='all')
  node.ids <- node.stats$ids$all
  num.temp <- node.stats$prev$all
  modes <- all$modes
  num.m1.temp <- ifelse(modes==1, node.stats$prev$all, node.stats$prev$m1)
  
  
  ## Update Nodes ##
  {
    if (modes == 1) {
      if (num.m1.temp > 0) {
        old.nw.size <- network.size(nw)
        new.births.m1 <- rpois(1, num.m1.temp * b.rate)
      } else {
        new.births.m1 <- 0
      }
      if (new.births.m1 > 0) {
        nw <- add.vertices(nw, nv=new.births.m1)
        new.nodes <- (old.nw.size+1):(old.nw.size+new.births.m1)
        nw <- activate.vertices(nw, onset=at, terminus=Inf, v=new.nodes)
      } else {
        new.nodes <- NULL
      }
    }
    if (modes == 2) {
      if (num.m1.temp > 0) {
        pct.m1.temp <- num.m1.temp / num.temp
        new.births.tot <- rpois(1, num.m1.temp * b.rate)
        new.births.m1 <- rbinom(1, new.births.tot, pct.m1.temp)
        new.births.m2 <- new.births.tot - new.births.m1
        
        old.size.m1 <- length(modeids(nw, 1))
        old.size.m2 <- length(modeids(nw, 2))
        old.size <- old.size.m1 + old.size.m2
        prefixes <- unique(substr(nw %v% 'vertex.names', 1, 1))
  
        if (new.births.m1 > 0) {
          new.pids.m1 <- paste0(prefixes[1], 
                               (old.size.m1 + 1):(old.size.m1 + new.births.m1))
          nw <- add.vertices(nw, new.births.m1, last.mode=FALSE, vertex.pid=new.pids.m1)
          new.nodes.m1 <- (old.size.m1+1):(old.size.m1+new.births.m1)
        } else {
          new.nodes.m1 <- NULL
        }
        if (new.births.m2 > 0) {
          new.pids.m2 <- paste0(prefixes[2], 
                               (old.size.m2 + 1):(old.size.m2 + new.births.m2))
          nw <- add.vertices(nw, new.births.m2, last.mode=TRUE, vertex.pid=new.pids.m2)
          new.size <- network.size(nw)
          new.nodes.m2 <- (new.size-new.births.m2+1):new.size
        } else {
          new.nodes.m2 <- NULL
        }
        new.nodes <- c(new.nodes.m1, new.nodes.m2)
        if (!is.null(new.nodes)) {
          nw <- activate.vertices(nw, onset=at, terminus=Inf, v=new.nodes)
        }
      } else {
        new.births.m1 <- new.births.m2 <- 0
        new.nodes <- NULL
      }
    }
  }
  

  ## Update Nodal Attributes ##
  if (length(new.nodes) > 0) {   
    # Update status
    if (tea == TRUE) {
      nw <- activate.vertex.attribute(nw, prefix='status', value=0, 
                                      onset=at, terminus=Inf, v=new.nodes)
    }
    if (modes == 1) {
      all$status <- c(all$status, rep(0, length(new.nodes)))
      if (all$save.statmat == TRUE) {
        new.mat <- matrix(NA, ncol=length(new.nodes), nrow=all$nsteps)
        all$stat.mat <- cbind(all$stat.mat, new.mat)
        all$stat.mat[at, new.nodes] <- 0
      }
    }
    if (modes == 2) {
      old.status.m1 <- all$status[1:old.size.m1]
        old.status.m2 <- all$status[(old.size.m1+1):old.size]
        new.status.m1 <- c(old.status.m1, rep(0, new.births.m1))
        new.status.m2 <- c(old.status.m2, rep(0, new.births.m2))
        all$status <- c(new.status.m1, new.status.m2)
      old.inf.time.m1 <- all$inf.time[1:old.size.m1]
        old.inf.time.m2 <- all$inf.time[(old.size.m1+1):old.size]
        new.inf.time.m1 <- c(old.inf.time.m1, rep(at, new.births.m1))
        new.inf.time.m2 <- c(old.inf.time.m2, rep(at, new.births.m2))
        all$inf.time <- c(new.inf.time.m1, new.inf.time.m2)
      if (all$save.statmat == TRUE) {
        stat.mat <- all$stat.mat
        stat.mat.m1 <- stat.mat[, 1:old.size.m1]
        stat.mat.m2 <- stat.mat[, (old.size.m1+1):old.size]
        if (new.births.m1 > 0) {
          new.mat.m1 <- matrix(NA, ncol=new.births.m1, nrow=all$nsteps)
          stat.mat.m1 <- cbind(stat.mat.m1, new.mat.m1)
        }
        if (new.births.m2 > 0) {
          new.mat.m2 <- matrix(NA, ncol=new.births.m2, nrow=all$nsteps) 
          stat.mat.m2 <- cbind(stat.mat.m2, new.mat.m2)
        }
        all$stat.mat <- cbind(stat.mat.m1, stat.mat.m2)
        all$stat.mat[at, new.nodes] <- 0
      }
    }
    
    
    # Set nodematch attributes
    nm.attr <- get.nodematch.attributes(nw, formation)
    if (!is.null(nm.attr)) {
      nw <- set.nodematch.attributes(nw, nm.attr, new.nodes)
    }
  }

  
  
  ## Set Summary Stats ##
  if (nw$gal$bipartite == FALSE) {
    if (at == 2) {
      all$b.flow.m1 <- c(0, new.births.m1)
    } else {
      all$b.flow.m1[at] <- new.births.m1
    }
  } else{
    if (at == 2) {
      all$b.flow.m1 <- c(0, new.births.m1)
      all$b.flow.m2 <- c(0, new.births.m2)
    } else {
      all$b.flow.m1[at] <- new.births.m1
      all$b.flow.m2[at] <- new.births.m2
    }
  }
  
  all$nw <- nw
  return(all)
}
