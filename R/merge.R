
#' @title Merge Model Simulations Across epiICM Objects
#'
#' @description This function merges data from two independent \code{epiSPM} 
#' simulations.
#'
#' @param x an \code{EpiModel} object of class \code{\link{epiICM}}.
#' @param y another \code{EpiModel} object of class \code{\link{epiICM}}, with the
#' identical model parameterization as \code{x}.
#' @param ...  additional arguments required by generic merge method, but not
#' used.
#' 
#' @details
#' The purpose of this generic merge function is to facilitate analysis of multiple
#' simulations of \code{\link{epiICM}} class models that may have been simulated 
#' under different function calls, but where the model parameterization between
#' the two calls is exactly the same. Such a situation would occur when one runs
#' sets of simulations in parallel across computing clusters for efficency. 
#' 
#' Note that this merge function does not work the same as the default merge in allowing
#' for a combined object where the structure differs between the individual elements:
#' instead, the function checks that objects are identical in model parameterization in
#' every respect except number of simulations. 
#' 
#' @method merge epiICM
#' @keywords extract
#' @export
#' 
#' @examples
#' x <- epiICM(type="SI", s.num=1000, i.num=100,
#'             trans.rate = 0.2, act.rate = 0.8,
#'             nsteps = 10, nsims = 3)
#' 
#' y <- epiICM(type="SI", s.num=1000, i.num=100,
#'             trans.rate = 0.2, act.rate = 0.8,
#'             nsteps = 10, nsims = 2)
#' 
#' z <- merge(x, y)
#' x$i.num
#' y$i.num
#' z$i.num
#'
merge.epiICM <- function(x, 
                         y,  
                         ...
                         ) {
  
  ## Check structure
  if (length(x) != length(y) || names(x) != names(y))
    stop("x and y have different structure")
  if (x$nsims > 1 & y$nsims > 1 & !all(sapply(x, class) == sapply(y, class)))
    stop("x and y have different structure")
    
  ## Check params
  x.metadata <- intersect(grep("num", names(x), invert=TRUE), 
                      grep("flow", names(x), invert=TRUE))
  x.data <- c(grep("num", names(x), invert=FALSE), 
              grep("flow", names(x), invert=FALSE))
  check.params <- x.metadata[-which(x.metadata %in% which(names(x) %in% c("nsims", "call")))]
  
  for (i in check.params) {
    check <- identical(x[[i]], y[[i]])
    if (check == FALSE) stop("x and y have different parameters")
  }
  
  z <- x
  new.range <- (x$nsims+1):(x$nsims+y$nsims)
  
  # Merge data
  for (i in x.data) {
    if (x$nsims == 1) {
      x[[i]] <- data.frame(x[[i]])
      names(x[[i]]) <- "sim1"
    }
    if (y$nsims == 1) {
      y[[i]] <- data.frame(y[[i]])
      names(y[[i]]) <- "sim1"
    }
    z[[i]] <- cbind(x[[i]], y[[i]])
    names(z[[i]])[new.range] <- paste("sim", new.range, sep="")
  }

  z$nsims <- max(new.range)
  
  return(z)
}

#' @title Merge Model Simulations Across epiNet.simTrans Objects
#'
#' @description This function merges data from two independent 
#' \code{epiNet.simTrans} simulations.
#'
#' @param x an \code{EpiModel} object of class \code{\link{epiNet.simTrans}}.
#' @param y another \code{EpiModel} object of class \code{\link{epiNet.simTrans}}, 
#'   with the identical model parameterization as \code{x}.
#' @param keep.trans if \code{TRUE}, keep in the merged output the transmission 
#'   data frames in \code{trans} and the status matrix in \code{stat.mat} (not only 
#'   not using temporally extended attributes) from the original \code{x} and \code{y}
#'   elements.
#' @param keep.network if \code{TRUE}, keep in the merged output the \code{networkDynamic}
#'   objects from the original \code{x} and \code{y} elements.
#' @param keep.stats if \code{TRUE}, keep in the merged output the network statistics
#'   (as set by the \code{stats.formula} parameter in \code{epiNet.simTrans}) from
#'   the original \code{x} and \code{y} elements.
#' @param ...  additional arguments required by generic merge method, but not
#'   used.
#' 
#' @details
#' The purpose of this generic merge function is to facilitate analysis of multiple
#' simulations of \code{\link{epiNet.simTrans}} class models that may have been simulated 
#' under different function calls, but where the model parameterization between
#' the two calls is exactly the same. Such a situation would occur when one runs
#' sets of simulations in parallel across computing clusters for efficency. 
#' 
#' Note that this merge function does not work the same as the default merge function
#' in generating a combined object where the structure differs between the 
#' individual elements: instead, the function checks that objects are identical 
#' in model parameterization in every respect except number of simulations. 
#' 
#' @method merge epiNet.simTrans
#' @keywords extract
#' @export
#' 
#' @examples
#' nw <- network.initialize(n = 100, directed = FALSE)
#' dissolution <- ~offset(edges)
#' coef.diss <- dissolution.coefs(dissolution, duration=10)
#' est <- epiNet.est(nw,
#'                   formation = ~ edges,
#'                   dissolution = ~offset(edges),
#'                   target.stats = 25,
#'                   coef.diss = coef.diss,
#'                   save.stats = FALSE, verbose = FALSE)
#' nwsims <- epiNet.simNet(est, nsteps = 10, verbose = FALSE)
#' 
#' x <- epiNet.simTrans(nwsims, type = "SI",
#'                      i.num = 1, trans.rate = 0.9,
#'                      sims.per.nw = 2, verbose = FALSE) 
#' y <- epiNet.simTrans(nwsims, type = "SI",
#'                      i.num = 1, trans.rate = 0.9,
#'                      sims.per.nw = 3, verbose = FALSE) 
#'                            
#' z <- merge(x, y)
#' x$i.num
#' y$i.num
#' z$i.num
#'
merge.epiNet.simTrans <- function(x, 
                                  y,                          
                                  keep.trans = TRUE, 
                                  keep.network = TRUE,
                                  keep.stats = TRUE,
                                  ...
                                  ) {
  
  ## Check structure
  if (length(x) != length(y) || names(x) != names(y))
    stop("x and y have different structure")
  x$nsims <- as.integer(x$nsims)
  y$nsims <- as.integer(y$nsims)
  if (x$nsims > 1 & y$nsims > 1 & !all(sapply(x, class) == sapply(y, class)))
    stop("x and y have different structure")
  
  ## Check params
  x.metadata <- intersect(grep("num", names(x), invert=TRUE), 
                          grep("flow", names(x), invert=TRUE))
  x.data <- c(grep("num", names(x), invert=FALSE), 
              grep("flow", names(x), invert=FALSE))
  check.params <- x.metadata[-which(x.metadata %in% 
                                    which(names(x) 
                                          %in% c("nsims", 
                                                 "call", 
                                                 "stat.mat",
                                                 "trans",
                                                 "network",
                                                 "timer",
                                                 "stats")))]
  
  for (i in check.params) {
    check <- identical(x[[i]], y[[i]])
    if (check == FALSE) stop("x and y have different parameters")
  }
  
  z <- x
  new.range <- (x$nsims+1):(x$nsims+y$nsims)
  
  # Merge data
  for (i in x.data) {
    if (x$nsims == 1) {
      x[[i]] <- data.frame(x[[i]])
      names(x[[i]]) <- "sim1"
    }
    if (y$nsims == 1) {
      y[[i]] <- data.frame(y[[i]])
      names(y[[i]]) <- "sim1"
    }
    z[[i]] <- cbind(x[[i]], y[[i]])
    names(z[[i]])[new.range] <- paste0("sim", new.range)
  }
  
  z$nsims <- max(new.range)
  
  if (keep.trans == TRUE) {
    for (i in new.range) {
      if (!is.null(x$stat.mat) & !is.null(y$stat.mat)) {
        z$stat.mat[[i]] <- y$stat.mat[[i-x$nsims]]
        if (!is.null(z$stat.mat))
          names(z$stat.mat)[i] <- paste("sim", i, sep="")
      }
      if (!is.null(x$trans) & !is.null(y$trans)) {
        z$trans[[i]] <- y$trans[[i-x$nsims]]
        if (!is.null(z$trans))
          names(z$trans)[i] <- paste0("sim", i)
      }
    }
  } else {
    z$stat.mat <- NULL
  }
  
  if (keep.network == TRUE & !is.null(x$network) & !is.null(y$network)) {
    for (i in new.range) {
      z$network[[i]] <- y$network[[i-x$nsims]]
      if (!is.null(z$network))
        names(z$network)[i] <- paste0("sim", i)
    }
  } else {
    z$network <- NULL
  }
  
  if (keep.stats == TRUE & !is.null(x$stats) & !is.null(y$stats)) {
    for (i in new.range) {
      z$stats[[i]] <- y$stats[[i-x$nsims]]
      if (!is.null(z$stats))
        names(z$stats)[i] <- paste0("sim", i)
    }
  }
  
  return(z)
}