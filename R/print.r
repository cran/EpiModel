
##' @S3method print epiNet.est
print.epiNet.est <- function(x, digits=3, ...) {
  
  cat("EpiModel Object")
  cat("\n=======================")
  cat("\nModel class:", class(x))
  estmeth <- ifelse(x$edapprox == TRUE, 'ERGM with Edges Approximation',
                                        'Full STERGM Fit')
  cat(paste("\nEsimation Method:", estmeth))
  
  cat("\n\nERGM Model Form")
  cat("\n-----------------------")
  cat("\nFormation: "); print(x$formation)
  cat("Dissolution: "); print(x$dissolution)
  cat("Constraints: "); cat(paste(as.character(x$constraints)[1], 
                                  as.character(x$constraints)[2], sep=''))
  
  if (any(names(x) == 'sim.stats')) {
    cat("\n\nFormation Diagnostics")
    cat("\n----------------------- \n")
    print(round(x$stats.formation, digits=digits))
    
    cat("\nDuration Diagnostics")
    cat("\n----------------------- \n")
    print(round(x$stats.duration, digits=digits))
  }
  cat("\n")
  
  invisible()
}

##' @S3method print epiNet.simNet
print.epiNet.simNet <- function(x, ...) {
  
  cat("EpiModel Object")
  cat("\n=======================")
  cat("\nModel class:", class(x))
  
  cat("\n\nNW Simulation Summary")
  cat("\n-----------------------")
  cat("\nNo. simulations:", x$nsims)
  cat("\nNo. time steps:", x$nsteps)

  cat("\n\nBase Network Object")
  cat("\n----------------------- \n")
  print(x$base.nw)
  cat("\n")
  
  invisible()
}

##' @S3method print epiNet.simTrans
print.epiNet.simTrans <- function(x, ...) {
  
  # model dimensions
  nts <- max(x$time)
  nsims <- x$nsims
  
  if (nsims == 1) {
    simnames <- 'sim1'
  }
  if (nsims == 2) {
    simnames <- 'sim1 sim2'
  }
  if (nsims > 2) {
    simnames <- paste('sim1 ... sim', nsims, sep='')
  } 
  
  cat("EpiModel Object")
  cat("\n=======================")
  cat("\nModel class:", class(x))
  
  cat("\n\nSimulation Summary")
  cat("\n-----------------------")
  cat("\nModel type:", x$type)
  cat("\nNo. simulations:", nsims)
  cat("\nNo. time steps:", nts)
  cat("\nNo. NW modes:", x$modes)
  
  cat("\n\nModel Output")
  cat("\n-----------------------")
  cat("\nCompartments:", names(x)[grep('num', names(x))], fill=60)
  cat("Flows:", names(x)[grep('flow', names(x))], fill=60)
  if (!(is.null(x$network))) cat("Networks:", simnames)
  if (!(is.null(x$stats))) cat("\nStats:", simnames)
  if (!(is.null(x$trans))) cat("\nTransmissions:", simnames)
  cat("")

  invisible()
}


##' @S3method print epiICM
print.epiICM <- function(x, ...) {
  
  # model dimensions
  nts <- max(x$time)
  nsims <- x$nsims
  
  cat("EpiModel Object")
  cat("\n=======================")
  cat("\nModel class:", class(x))
  
  cat("\n\nSimulation Summary")
  cat("\n-----------------------")
  cat("\nModel type:", x$type)
  cat("\nNo. simulations:", nsims)
  cat("\nNo. time steps:", nts)
  cat("\nNo. groups:", x$groups)
  
  cat("\n\nModel Parameters")
  cat("\n-----------------------\n")
  for (i in 1:length(x$params)) {
    cat(names(x$params)[i], "=", x$params[[i]], fill=60)
  }
  
  cat("\nModel Output")
  cat("\n-----------------------")
  cat("\nCompartments:", names(x)[grep('num', names(x))], fill=60)
  cat("Flows:", names(x)[grep('flow', names(x))], fill=60)
  
  invisible()
}


##' @S3method print epiDCM
print.epiDCM <- function(x, ...) {
  
  # model dimensions
  nts <- max(x$time)
  nruns <- x$nruns
  
  if (x$type == 'SI') x$params$rec.rate <- NULL
  if (x$type == 'SI' && x$groups == 2) x$params$rec.rate.g2 <- NULL
  
  cat("EpiModel Object")
  cat("\n=======================")
  cat("\nModel class:", class(x))
  
  cat("\n\nSimulation Summary")
  cat("\n-----------------------")
  cat("\nModel type:", x$type)
  cat("\nNo. runs:", nruns)
  cat("\nNo. time steps:", nts)
  cat("\nNo. groups:", x$groups)
  
  cat("\n\nModel Parameters")
  cat("\n-----------------------\n")
  for (i in 1:length(x$params)) {
    cat(names(x$params)[i], "=", x$params[[i]], fill=60)
  }
  
  cat("\nModel Output")
  cat("\n-----------------------")
  cat("\nCompartments:", names(x)[grep('num', names(x))], fill=60)
  cat("Flows:", names(x)[grep('flow', names(x))], fill=60)
  
  invisible()
}


##' @S3method print dissolution.coefs
print.dissolution.coefs <- function(x, ...) {
  
  cat("Dissolution Coefficients")
  cat("\n=======================")
  cat("\nDissolution Model: "); ; print(x$dissolution)
  cat("Edge Duration:", x$duration)
  cat("\nAdjusted Coefficient:", x$coef.adj)
  cat("\nCrude Coefficient:", x$coef.crude)
  
}