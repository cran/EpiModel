
##
## Recovery and other non-infection state transition modules for epiNet.simTrans
##

#' @title Recovery: epiNet.simTrans Module
#'
#' @description This function simulates recovery from the infected state
#'   either to an distinct recovered state (SIR model type) or back to a
#'   susceptible state (SIS model type), for use in \code{\link{epiNet.simTrans}}.
#' 
#' @param all a list object containing a \code{networkDynamic} object and other
#'   initialization information passed from \code{\link{epiNet.simTrans}}.
#' @param at time step to query. 
#' 
#' @export
#' @keywords internal
#'
recovery <- function(all, at) {
  
  nw <- all$nw
  tea <- ifelse(any(names(nw$val[[1]]) %in% 'status.active'), TRUE, FALSE)
  
  inf.temp <- get.stat(all, stat=1, at=at, out='all')
  if (inf.temp$prev > 0) {
    new.recovs <- rpois(1, all$rec.rate * inf.temp$prev)
    if (new.recovs > inf.temp$prev) new.recovs <- inf.temp$prev
    d.rec <- ssample(inf.temp$ids, new.recovs)
    if (new.recovs > 0) {
      if (all$type == 'SIR') {
        all$status[d.rec] <- 2
        if (tea == TRUE) {
          nw <- activate.vertex.attribute(nw, 'status', value=2, onset=at, 
                                    terminus=Inf, v=d.rec)
        }
        if (all$save.statmat == TRUE) {
          all$stat.mat[at, ] <- all$status
        }
      }
      if (all$type == 'SIS') {
        all$status[d.rec] <- 0
        if (tea == TRUE) {
          nw <- activate.vertex.attribute(nw, 'status', value=0, onset=at, 
                                    terminus=Inf, v=d.rec)
        }
        if (all$save.statmat == TRUE) {
          all$stat.mat[at, ] <- all$status
        }
      }
    }
  } else {
    d.rec <- NULL
  }
  
  ## Save recovered vector
  if (nw %n% 'bipartite' == FALSE) {
    if (at == 2) {
      all$recovs.m1 <- c(0, length(d.rec))
    } else {
      all$recovs.m1[at] <- length(d.rec)
    }
  } else{
    bip.cutoff <- nw %n% 'bipartite'
    if (at == 2) {
      all$recovs.m1 <- c(0, length(d.rec[d.rec <= bip.cutoff]))
      all$recovs.m2 <- c(0, length(d.rec[d.rec > bip.cutoff]))
    } else {
      all$recovs.m1[at] <- length(d.rec[d.rec <= bip.cutoff])
      all$recovs.m2[at] <- length(d.rec[d.rec > bip.cutoff])
    }
  }

  all$nw <- nw
  return(all)
}
