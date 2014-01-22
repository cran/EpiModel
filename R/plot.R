

# Main Exported Methods ---------------------------------------------------

#' @title Plot Values from a Deterministic Compartmental Epidemic Model
#' 
#' @description This function plots values from an deterministic compartment 
#'  epidemic model solved with \code{epiDCM}.
#'
#' @param x an \code{EpiModel} object of class \code{\link{epiDCM}}.
#' @param y output compartments or flows from \code{epiDCM} object to plot.
#' @param popfrac if \code{TRUE}, plot prevalence of values rather than numbers 
#'   (see details).
#' @param run run number to plot for models with multiple runs; with a default 
#'   of run 1.
#' @param col color for lines, either specified as a single color in a standard 
#'   R color format, or alternatively, a color palette from \code{\link{RColorBrewer}} 
#'   (see details).
#' @param lwd line width for output values.
#' @param lty line type for output values.
#' @param alpha transparency level for lines, where 0 = transparent and 1 = opaque 
#'  (see \code{\link{transco}}).
#' @param leg type of legend to plot. Values are "n" for no legend, "full" for 
#'   full legend, and "lim" for limited legend (see details).
#' @param leg.name character string to use for legend, calculated automatically 
#'   if not supplied.
#' @param leg.cex legend scale size, with default of 0.8.
#' @param xlim x-axis scale limits for plot, with default calculated based on 
#'   model time steps.
#' @param ylim y-axis scale limits for plot, with default calculated based on 
#'   range of data.
#' @param main character string for main plot title.
#' @param axs plot axis type (see \code{\link{par}} for details), with default 
#'   to "r".
#' @param add if \code{TRUE}, new plot is not called and lines are added to 
#'   existing plot.
#' @param ... additional arguments to pass to main plot (see \code{\link{plot.default}}).
#' 
#' @details
#' The goal of this plotting function is to facilitate plotting values from a 
#' deterministic compartmental model solved with \code{\link{epiDCM}}. Depending 
#' on the number of model runs (sensitivity analyses) and number of groups, the 
#' default plot without any options set will generally plot the disease prevalence 
#' from the model. The specific compartments or flows to plot may be set using 
#' the \code{y} parameter, and in multiple run models the specific run may also 
#' be set. 
#' 
#' Compartment prevalences are the size of a compartment over some denominator. 
#' To plot the raw numbers from any compartment, use \code{popfrac=FALSE}; this 
#' is the default for any plots of flows. The \code{popfrac} parameter calculates 
#' and plots the denominators of all specified compartments using these rules: 1) 
#' for one-group models, the prevalence of any compartment is the compartment size 
#' divided by the total population size; 2) for two-group models, the prevalence 
#' of any compartment is the compartment size divided by the group population size. 
#' 
#' Since \code{\link{epiDCM}} supports multiple run sensitivity models, plotting 
#' the results of such models must use a complex visual scheme for easy 
#' representation of results. This is accomplished mainly using the 
#' \code{\link{RColorBrewer}} color palettes, in which one may specify a range of 
#' linked colors using named palettes. For \code{plot.epiDCM}, one may either 
#' specify a Brewer color palette listed in \code{\link{display.brewer.all}}, or 
#' alternatively a single or vector of standard R colors (named, hexidecimal, or 
#' positive integers; see \code{\link{col2rgb}}). 
#' 
#' There are three automatic legend types available, and the legend is 
#' automatically added by default for most plots. To shut off the legend, use 
#' \code{leg="n"}. To plot a legend with values for every line, use 
#' \code{leg="full"}. With sensitivity models with many runs, this is generally 
#' not visually appealing Therefore, in those cases, use \code{leg="lim"} to 
#' plot a legend limited to the highest and lowest of the varying parameters in 
#' the model. In cases where the default legend names are not helpful, one may 
#' override those names with the \code{leg.names} argument.
#' 
#' @method plot epiDCM
#' @export 
#' 
#' @keywords plot
#' @seealso \code{\link{epiDCM}}, \code{\link{display.brewer.all}}, 
#'  \code{\link{transco}}
#' 
#' @examples
#' # Deterministic SIR model with varying act rate
#' mod <- epiDCM(type="SIR", s.num=1000, i.num=1, r.num=0, 
#'               trans.rate=0.2, act.rate=1:10, rec.rate=1/3,
#'               b.rate=0.011, ds.rate=0.01, di.rate=0.03, 
#'               dr.rate=0.01, nsteps=500, verbose=TRUE)
#'
#' # Plot disease prevalence by default
#' plot(mod)
#'
#' # Plot prevalence of susceptibles
#' plot(mod, y="s.num", col="Greys")
#' 
#' # Plot number of susceptibles
#' plot(mod, y="s.num", popfrac=FALSE, col="Greys")
#' 
#' # One way to plot multiple runs of multiple compartments together
#' plot(mod, y=c("s.num", "i.num"), run=5, xlim=c(0,200))
#' plot(mod, y=c("s.num", "i.num"), run=10, alpha=0.3, leg="n", add=TRUE)
#' 
plot.epiDCM <- function(x, 
                        y, 
                        popfrac,
                        run, 
                        col, 
                        lwd, 
                        lty, 
                        alpha,
                        leg, 
                        leg.name, 
                        leg.cex,  
                        xlim, 
                        ylim, 
                        main,
                        axs, 
                        add = FALSE, 
                        ...) {
  
  
  ## Set missing flags
  if (missing(y)) noy <- TRUE else noy <- FALSE
  if (missing(run)) norun <- TRUE else norun <- FALSE
  if (missing(col)) nocol <- TRUE else nocol <- FALSE
  if (missing(lwd)) nolwd <- TRUE else nolwd <- FALSE
  if (missing(lty)) nolty <- TRUE else nolty <- FALSE
  if (missing(leg)) noleg <- TRUE else noleg <- FALSE
  if (missing(popfrac)) nopopfrac <- TRUE else nopopfrac <- FALSE
  
  ## Model dimensions
  nts <- max(x$time)
  nvars <- length(x)
  nruns <- x$nruns
  if (norun == FALSE && any(run > nruns)) {
    stop("Run ", max(run), " was requested, but only ", nruns, " runs available in object.")
  }
  
  groups <- x$groups
  type <- x$type
  
  ## Universal defaults
  if (noleg == TRUE) leg <- "n"
  if (missing(alpha)) alpha <- 0.9
  if (missing(leg.cex)) leg.cex <- 0.8
  if (missing(main)) main <- paste("DCM", x$type, "Model")
  
  
  ## Defaults for missing y 
  if (noy == TRUE && nruns == 1) {
    y <- grep(".num", names(x), value=TRUE)
  }
  if (noy == TRUE && nruns > 1) {
    y <- grep("i.num", names(x), value=TRUE)
  }
  if (all(y %in% names(x)) == FALSE) 
    stop("Specified y is unavailable.")
  lcomp <- length(y)
  
  
  ## Prevalence calculations
  if (missing(popfrac)) popfrac <- TRUE
  if (any(grepl(".flow", y))) popfrac <- FALSE
  x <- denom(x, y, popfrac)
  
  
  ## Compartment ymax calculations
  if (popfrac == FALSE) {
    allmax <- sapply(1:lcomp, function(i) max(x[[y[i]]]))
    ymax <- ceiling(max(allmax))
  } else {
    ymax <- 1
  }
  
  
  ## Defaults for ylim, xlim, axs
  if (missing(ylim)) ylim <- c(0, ymax)
  if (missing(xlim)) xlim <- c(0, nts)
  if (missing(axs)) axs <- "r"
  
  
  ## Defaults for lwd
  if (nolwd == FALSE && lcomp > 1 && length(lwd) < lcomp)
    lwd <- rep(lwd, lcomp)
  if (nolwd == FALSE && lcomp == 1 && length(lwd) < nruns)
    lwd <- rep(lwd, nruns)
  if (nolwd == TRUE) {
    lwd <- rep(2.5, lcomp*nruns)
  }
  
  
  ## Defaults for lty
  if (nolty == FALSE && lcomp > 1 && length(lty) < lcomp) 
    lty <- rep(lty, lcomp)
  if (nolty == FALSE && lcomp == 1 && length(lty) < nruns)
    lty <- rep(lty, nruns)
  if (nolty == TRUE) {
    lty <- rep(1, lcomp*nruns)
    if (groups == 2 && noy == TRUE) {
      lty <- rep(1:2, times=lcomp/2)
    }
  }
  
  
  ## Main plot window
  if (add == FALSE) {
    if (popfrac == FALSE) {
      Time <- Number <- 1
      plot(Time, Number, type="n", bty="n", 
           xaxs=axs, yaxs=axs, xlim=xlim, ylim=ylim, 
           main=main, ...)
    } else {
      Time <- Prevalence <- 1
      plot(Time, Prevalence, type="n", bty="n", 
           xaxs=axs, yaxs=axs, xlim=xlim, ylim=ylim, 
           main=main, ...)
    }  
  }
  
  
  ## Default line colors
  pal <- NULL
  # Missing col
  if (nocol == TRUE) {
    if (lcomp == 1) {
      if (nruns == 1) col <- "black"
      if (nruns > 1) col <- "Set1"
      if (nruns > 5) col <- "Spectral"
      if (norun == FALSE && length(run) == 1) col <- "black"
    }
    if (lcomp > 1) {
      col <- "Set1"
    }
  }
  
  # Not missing col, poker chips palette for MM
  if (nocol == FALSE && col[1] == "chips") {
    col <- c("#377EB8", "#E41A1C", "darkgrey")
  } 
  
  # Test if using a RColorBrewer palette
  if (length(col) == 1 && col %in% row.names(brewer.pal.info)) 
    use.brewer <- TRUE else use.brewer <- FALSE 
  
  # Set color palette
  if (is.null(pal)) {
    if (lcomp == 1) {
      if (use.brewer == TRUE) {
        if (nruns < 6) {
          pal <- transco(brewer.pal(5, col)[1:nruns], alpha)
        } else {
          pal <- transco(brewer.ramp(nruns, col), alpha)
        }
      } 
      if (use.brewer == FALSE) {
        pal <- transco(rep(col, nruns), alpha) 
      }
    }
    if (lcomp > 1) {
      if (use.brewer == TRUE) {
        if (lcomp > 4) {
          pal <- transco(brewer.ramp(lcomp, col), alpha)
        } else {
          pal <- transco(brewer.pal(max(c(lcomp, 4)), col), alpha)
          fixpal <- pal
          fixpal[1] <- pal[2]; fixpal[2] <- pal[1]
          pal <- fixpal
        }
        if (groups == 2 && noy == TRUE) {
          pal <- transco(brewer.pal(max(c(lcomp, 4)), col), alpha)
          fixpal <- pal
          fixpal[1] <- pal[2]; fixpal[2] <- pal[1]
          pal <- fixpal
          pal <- rep(pal, each=2)
        }
      } 
      if (use.brewer == FALSE) {
        pal <- transco(rep(col, lcomp), alpha)
        if (groups == 2 && noy == TRUE) {
          pal <- transco(rep(col, each=2), alpha)
        }
      }
    }
  }
  
  
  ## Plot lines
  if (lcomp == 1) {
    if (nruns == 1) {
      lines(x$time, x[[y]][,1], 
            lwd=lwd[1], lty=lty[1], col=pal[1])
    } 
    if (nruns > 1) {
      if (norun == TRUE) {
        for (i in 1:nruns) { 
          lines(x$time, x[[y]][,i], 
                lwd=lwd[i], lty=lty[i], col=pal[i])
        }
      } else {
        if (length(run) == 1) {
          lines(x$time, x[[y]][,run], 
                lwd=lwd[1], lty=lty[1], col=pal[1])
        }
        if (length(run) > 1) {
          for (i in 1:length(run)) {
            lines(x$time, x[[y]][,run[i]], 
                  lwd=lwd[i], lty=lty[i], col=pal[i])
          }
        }
      }
    }
  }
  if (lcomp > 1) {
    if (nruns == 1) {
      for (i in 1:lcomp) {
        lines(x$time, x[[y[i]]][,1], 
              lwd=lwd, lty=lty[i], col=pal[i])
      }
    } 
    if (nruns > 1) {
      if (norun == TRUE) {
        for (i in 1:lcomp) {
          run <- 1
          lines(x$time, x[[y[i]]][,run], 
                lwd=lwd[i], lty=lty[i], col=pal[i])
        }
      } 
      if (norun == FALSE) {
        if (length(run) > 1)
          stop("Plotting multiple runs of multiple outputs (y) is not supported")
        for (i in 1:lcomp) {
          lines(x$time, x[[y[i]]][,run], 
                lwd=lwd[i], lty=lty[i], col=pal[i])
        } 
      }
    }
  }
  
  
  ## Legend
  
  # Default legend type
  if (noleg == TRUE) {
    leg <- "n"
    if (lcomp == 1 & nruns < 3) leg <- "full"
    if (lcomp == 1 & nruns >= 3) leg <- "lim"
    if (lcomp > 1) leg <- "full"
  } else {
    if (leg == "lim" & nruns < 3) leg <- "full"
    if (leg == "lim" & lcomp == 2) leg <- "full"
  }
  
  # Default legend names
  if (missing(leg.name)) {
    if (nruns == 1) leg.names <- y
    if (nruns > 1) {
      if (norun == TRUE & lcomp == 1) leg.names <- names(x[[y[1]]])
      if (norun == FALSE & lcomp == 1) {
        if (length(run) == 1) leg.names <- y
        if (length(run) > 1) leg.names <- names(x[[y[1]]][run])
      }
      if (lcomp > 1) leg.names <- y
    }
  } else {
    if (lcomp == 1) leg.names <- paste(leg.name, 1:nruns, sep=" ")
    if (lcomp > 1) {
      leg.names <- y
      warning("Legend names ignored for multiple outputs (y) plots of multiple run models")
    }
  }
  
  # Legend
  if (leg == "n") { }
  if (norun == TRUE) {
    if (leg == "full") {
      legend("topright", legend=leg.names, 
             bg="white", lty=lty, lwd=lwd,
             col=pal, cex=leg.cex)
    } 
    if (leg == "lim") {
      legend("topright", legend=c(leg.names[1], "...", leg.names[nruns]), 
             bg="white", lty=c(lty[1], lty[nruns]), lwd=lwd+1,
             col=c(pal[1], "white", pal[nruns]), cex=leg.cex) 
    }
  } 
  if (norun == FALSE & leg != "n") {
    if (lcomp == 1) {
      legend("topright", legend=leg.names, 
             bg="white", lty=lty[1:length(run)], lwd=lwd[1:length(run)],
             col=pal[1:length(run)], cex=leg.cex)
    }
    if (lcomp > 1) {
      legend("topright", legend=leg.names, 
             bg="white", lty=lty, lwd=lwd,
             col=pal, cex=leg.cex)
    }
  }
  
  ## Mtext
  if (add == FALSE) {
    if ((noy == TRUE & lcomp == 1) | (nruns > 1 & lcomp > 1)) {
      if (noy == TRUE) mt <- paste("Plotting ", y)
      if (nruns > 1 & lcomp > 1) mt <- paste("Run =", run)
      mtext(mt, 3, line=0, cex=0.75)
    }
  }
  
}


#' @title Plot Simulations from a Individual Contact Epidemic Model
#'
#' @description This function plots values from an individual contact epidemic model 
#'  simulated with \code{epiICM}.
#'
#' @param x an \code{EpiModel} object of class \code{epiICM}.
#' @param y output compartments or flows from \code{epiICM} object to plot.
#' @param popfrac if \code{TRUE}, plot prevalence of values rather than numbers 
#'   (see details).
#' @param sim.lines if \code{TRUE}, plot individual simulation lines.
#' @param sims a vector representing which individual simulation lines to plot, 
#'  with default to plot all lines.
#' @param sim.col a vector of any standard R color format for simulation lines.
#' @param sim.lwd line width for simulation lines.
#' @param sim.alpha transparency level for simulation lines, where 0 = transparent 
#'   and 1 = opaque (see \code{\link{transco}}).
#' @param mean.line if \code{TRUE}, plot mean of simulations across time.
#' @param mean.extinct if \code{TRUE}, include extinct simulations in mean 
#'   calculation (see details).
#' @param mean.col a vector of any standard R color format for mean lines.
#' @param mean.lwd line width for mean lines.
#' @param mean.lty line type for mean lines.
#' @param qnts if numeric, plot polygon of simulation quantiles based on the 
#'   range implied by the argument (see details). If \code{FALSE}, supress polygon 
#'   plot.
#' @param qnts.col a vector of any standard R color format for polygons.
#' @param qnts.alpha transparency level for quantile polygons, where 0 = transparent 
#'   and 1 = opaque (see \code{\link{transco}}).
#' @param leg if \code{TRUE}, plot default legend.
#' @param leg.cex a numeric value to scale the legend size.
#' @param xlim x-axis scale limits for plot, with default calculated based on 
#'   model time steps.
#' @param ylim y-axis scale limits for plot, with default calculated based on 
#'   range of data.
#' @param main character string for main plot title.
#' @param axs plot axis type (see \code{\link{par}} for details), with default 
#'   to \code{"r"}.
#' @param add if \code{TRUE}, new plot is not called and lines are added to 
#'   existing plot.
#' @param ... additional arguments to pass to main plot (see \code{\link{plot.default}}).
#' 
#' @details
#' This plotting function will extract the simulation output values from an 
#' object of class \code{epiICM} and plot the requested time series data of 
#' disease prevalence and other results. The summary statistics that the function 
#' calculates and plots are individual simulation lines, means of the individual 
#' simulation lines, and quantiles of those individual simulation lines. The mean 
#' line, toggled on with \code{mean.line=TRUE} is calculated as the row mean 
#' across simulations at each time step.
#' 
#' Compartment prevalences are the size of a compartment over some denominator. 
#' To plot the raw numbers from any compartment, use \code{popfrac=FALSE}; this 
#' is the default for any plots of flows. The \code{popfrac} parameter calculates 
#' and plots the denominators of all specified compartments using these rules: 1) 
#' for one-group models, the prevalence of any compartment is the compartment size 
#' divided by the total population size; 2) for two-group models, the prevalence 
#' of any compartment is the compartment size divided by the group population size.
#' 
#' The quantiles show the range of outcome values within a certain specified 
#' quantile range. By default, the interquartile range is shown. This is 
#' specified by \code{qnts=0.5}, that is the middle 50% of the data. If one 
#' wanted to show the middle 95% of the date, one would specify \code{qnts=0.95}. 
#' To toggle off the polygons where they are plotted by default, specify 
#' \code{qnts=FALSE}.
#' 
#' @method plot epiICM
#' @export 
#' 
#' @keywords plot
#' @seealso \code{\link{epiICM}}
#' 
#' @examples
#' \dontrun{
#' # Plotting multiple compartment values automatically
#' mod <- epiICM(type="SIR", s.num=500, i.num=1,
#'              trans.rate=0.2, act.rate=0.25, rec.rate=1/50,
#'              nsteps=1000, nsims=10)
#' plot(mod)
#' 
#' # Plot only infecteds and show individual simulation lines
#' mod <- epiICM(type="SI", s.num=500, i.num=1, 
#'               trans.rate=0.2, act.rate=0.2, 
#'               nsteps=500, nsims=10)
#' plot(mod, y="i.num", sim.lines=TRUE)
#' }
#' 
plot.epiICM <- function(x, 
                        y, 
                        popfrac,
                        sim.lines = FALSE, 
                        sims,
                        sim.col, 
                        sim.lwd,
                        sim.alpha,
                        mean.line = TRUE, 
                        mean.extinct = TRUE, 
                        mean.col, 
                        mean.lwd, 
                        mean.lty,
                        qnts, 
                        qnts.col, 
                        qnts.alpha,
                        leg,
                        leg.cex,
                        xlim, 
                        ylim, 
                        main,
                        axs, 
                        add = FALSE, 
                        ...) {
  
  # model dimensions
  nts <- max(x$time)
  nsims <- x$nsims
  if (missing(sims)) sims <- 1:nsims
  if (max(sims) > x$nsims) stop("Set sim to between 1 and ", x$nsims)
  if (missing(mean.lty)) nomeanlty <- TRUE else nomeanlty <- FALSE
  
  type <- x$type
  if (class(x) == "epiICM") modes <- x$groups
  if (class(x) == "epiNet.simTrans") modes <- x$modes
  
  # compartments
  if (missing(y)) nocomp <- TRUE else nocomp <- FALSE
  if (nocomp == TRUE) {
    y <- grep(".num", names(x), value=TRUE)
    if (missing(leg)) leg <- TRUE
  }
  if (nocomp == FALSE) {
    if (any(y %in% names(x) == FALSE))
      stop("Specified y is not available in object.")
  }
  lcomp <- length(y)
  
  # Color palettes
  if (modes == 2 & nocomp == TRUE) {
    pal <- brewer.pal(3, "Set1")
    if (missing(mean.col)) {
      if (type == "SIR") {
        mean.col <- c(pal[2], pal[1], pal[3], pal[2], pal[1], pal[3])
      } else {
        mean.col <- c(pal[2], pal[1], pal[2], pal[1])
      }
    }
    if (missing(qnts.col)) {
      if (type == "SIR") {
        qnts.col <- c(pal[2], pal[1], pal[3], pal[2], pal[1], pal[3])
      } else {
        qnts.col <- c(pal[2], pal[1], pal[2], pal[1])
      }
    }
  } 
  if (missing(sim.lwd)) sim.lwd <- rep(0.75, lcomp)
  if (missing(sim.col)) sim.col <- rep("black", lcomp)
  if (missing(sim.alpha) & nsims == 1) sim.alpha <- 0.9
  if (missing(sim.alpha) & nsims > 1) sim.alpha <- max(c(0.05, 1-log10(nsims)/3))
  sim.pal <- transco(sim.col, sim.alpha)
  bpal <- brewer.pal(3, "Set1")
  bpal <- c(bpal[2], bpal[1], bpal[3])
  if (!(missing(mean.col)) && mean.col == "chips") {
    mean.col <- c(bpal[1], bpal[2], "darkgrey")
    if (modes == 2 & nocomp) {
      mean.col <- rep(mean.col, times=2)
    }
  }
  if (missing(mean.col)) mean.col <- bpal
  mean.pal <- transco(mean.col, 0.9)
  if (!(missing(qnts.col)) && qnts.col == "chips") {
    qnts.col <- c(bpal[1], bpal[2], "darkgrey")
    if (modes == 2 & nocomp) {
      if (type == "SIR") {
        qnts.col <- rep(qnts.col, times=2)
      } else {
        qnts.col <- rep(qnts.col[1:2], times=2) 
      }
    }
  }
  if (missing(qnts.col)) qnts.col <- bpal
  if (missing(qnts.alpha)) qnts.alpha <- 0.4
  qnts.pal <- transco(qnts.col, qnts.alpha)
  
  
  ## Prevalence calculations
  if (missing(popfrac)) popfrac <- TRUE
  if (any(grepl(".flow", y))) popfrac <- FALSE
  x <- denom(x, y, popfrac)
  
  
  # Compartment max
  if (popfrac == FALSE) {
    if (lcomp == 1) {
      max.prev <- max(x[[y]])
    } else {
      max.prev <- max(sapply(y, function(comps) max(x[[comps]])))
    }
  } else {
    max.prev <- 1
  }
  
  # Missing args
  if (missing(xlim)) xlim <- c(0, nts)
  if (missing(ylim)) ylim <- c(0, max.prev)
  
  if (missing(main)) {
    if (class(x) == "epiICM") modclass <- "ICM"
    if (class(x) == "epiNet.simTrans") modclass <- "Network"
    main <- paste(modclass, x$type, "Model")
  }
  
  # Main plot
  if (missing(axs)) axs <- "r"
  if (add == FALSE) {
    if (popfrac == FALSE) {
      Time <- Number <- 1
      plot(Time, Number, type="n", bty="n", 
           xaxs=axs, yaxs=axs, xlim=xlim, ylim=ylim, 
           main=main, ...)
    } else {
      Time <- Prevalence <- 1
      plot(Time, Prevalence, type="n", bty="n", 
           xaxs=axs, yaxs=axs, xlim=xlim, ylim=ylim, 
           main=main, ...)
    }
  }
  
  # Quantiles
  if (missing(qnts)) disp.qnts <- FALSE else disp.qnts <- TRUE
  if (nsims == 1) disp.qnts <- FALSE
  if (modes == 1 & nocomp == TRUE & missing(qnts)) {
    disp.qnts <- TRUE
    qnts <- 0.5
  }
  if (disp.qnts) {
    if (qnts > 1 | qnts < 0) stop("qnts must be between 0 and 1")
    for (j in seq_len(lcomp)) {
      quants <- c((1-qnts)/2, 1-((1-qnts)/2))
      qnt.prev <- apply(x[[y[j]]], 1, 
                        function(x) quantile(x, c(quants[1], quants[2])))
      xx <- c(1:nts, nts:1)
      yy <- c(qnt.prev[1,], rev(qnt.prev[2,]))
      polygon(xx, yy, col=qnts.pal[j], border=NA)
    }
  }
  
  # Individual sim lines
  if (sim.lines) {
    if (nsims == 1) {
      for (j in seq_len(lcomp)) {
        lines(1:nts, x[[y[j]]][,1], lwd=sim.lwd, col=sim.pal[j])
      }
    }
    if (nsims > 1) {
      for (j in seq_len(lcomp)) {
        for (i in sims) {
          lines(1:nts, x[[y[j]]][,i], lwd=sim.lwd, col=sim.pal[j])
        }
      }
    } 
  }
  
  # Mean of sim lines
  if (!missing(mean.lwd) && length(mean.lwd) < lcomp) 
    mean.lwd <- rep(mean.lwd, lcomp)
  if (missing(mean.lwd)) mean.lwd <- rep(2.5, lcomp)
  if (nomeanlty == TRUE) {
    if (nocomp == TRUE) {
      if (modes == 1) {
        mean.lty <- rep(1, lcomp)
      }
      if (modes == 2) {
        mean.lty <- rep(1:2, each=lcomp/2)
      } 
    }
    if (nocomp == FALSE) {
      mean.lty <- rep(1, lcomp)
    }
  }
  if (nomeanlty == FALSE & length(mean.lty) < lcomp) {
    mean.lty <- rep(mean.lty, lcomp)
  } 
  if (mean.line == TRUE) {
    if (nsims == 1) {
      for (j in seq_len(lcomp)) {
        mean.prev <- x[[y[j]]][,1]
        lines(1:nts, mean.prev, lwd=mean.lwd[j], col=mean.pal[j], lty=mean.lty[j])
        
      }
    }
    if (nsims > 1) {
      for (j in seq_len(lcomp)) {
        if (mean.extinct == TRUE) {
          mean.prev <- apply(x[[y[j]]], 1, mean)
        } else {
          non.extinct <- as.vector(which(apply(x$si.flow, 2, max) > 0))
          mean.prev <- apply(x[[y[j]]][,non.extinct], 1, mean)
        }
        lines(1:nts, mean.prev, lwd=mean.lwd[j], col=mean.pal[j], lty=mean.lty[j])
      }
    } 
  }
  
  if (missing(leg)) leg <- FALSE
  if (missing(leg.cex)) leg.cex <- 0.9
  if (leg == TRUE) {
    if (modes == 2 & nocomp) {
      legend("topright", legend=y, lty=mean.lty, lwd=3, 
             col=mean.pal, cex=leg.cex, bg="white")
    } else {
      legend("topright", legend=y, lty=1, lwd=3, 
             col=mean.pal, cex=leg.cex, bg="white")
    }
  }
  
}


#' @title Plot Diagnostics from an epiNet.est Object
#' 
#' @description This function plots values from diagnostic simulations
#'  in \code{epiNet.est}.
#'
#' @param x an \code{EpiModel} object of class \code{\link{epiNet.est}}.
#' @param type plot type, with options of \code{type="formation"} for 
#'   partnership formation statistics or other network statistics specified in 
#'   \code{epiNet.est}, or \code{type="duration"} to plot the mean ages of 
#'   partnerships over time.
#' @param dx.start start time for diagnostic plots. This must be a positive integer.
#' @param dx.end end time for diagnostic plots. This must be less than or equal 
#'   to the number of time steps simulated in the \code{epiNet.est} diagnostics.
#' @param dx.leg if \code{TRUE}, show legend (only if \code{plots.joined=TRUE})
#' @param plots.joined if \code{TRUE}, combine all target statistics in one plot,
#'  versus one plot per target statistic if \code{FALSE}
#' @param ... additional arguments to pass to main plot window
#' 
#' @details
#' The plot function for \code{epiNet.est} objects will generate plots of two
#' types of model diagnostic statistics that were run as part of the diagnostic
#' tools within that estimation function. The \code{formation} plot shows the 
#' summary statistics requested in \code{stats.form}, which defaults to those in 
#' the formation formula. The \code{duration} plot shows the average age of all 
#' partnerships at each time step, up until the maximum time step requested. This 
#' is calculated with the \code{\link{edgelist.meanage}} function.
#' 
#' The \code{plots.joined} argument will control whether the target statistics 
#' in the \code{formation} plot are joined in one plot or plotted separately. 
#' The default is based on the number of network statistics requested. The 
#' layout of the separate plots within the larger plot window is also based on 
#' the number of statistics.
#' 
#' Required for these plots is that the estimation diagnostics are run in 
#' \code{epiNet.est}. This happens by default, and is set with the \code{stats=TRUE} 
#' argument in that function. Since these diagnostics only generate one network 
#' simulation, the plots here will only show one line per network statistic. If the 
#' range of possible statistics is of interest, then the fitted STERGM should be 
#' simulated multiple times, which is possible with \code{\link{epiNet.simNet}}, 
#' which has its own plotting function, \code{\link{plot.epiNet.simNet}}.
#' 
#' @method plot epiNet.est
#' @export
#' 
#' @keywords plot
#' @seealso \code{\link{epiNet.est}}
#' 
#' @examples
#' \dontrun{
#' # Initialize bipartite network
#' nw <- network.initialize(500, bipartite=250, directed=FALSE)
#' 
#' # Fit model with limited degree terms, but monitor wider range
#' est <- epiNet.est(nw, 
#'                   formation = ~ edges + b1degree(0:1) + b2degree(0:1), 
#'                   dissolution = ~ offset(edges), 
#'                   target.stats = c(165, 100, 137.5, 120, 102.5), 
#'                   coef.diss = dissolution.coefs(~ offset(edges), duration=25),  
#'                   stats.form = ~ edges + b1degree(0:5) + b2degree(0:5))
#' 
#' # Formation plot with plots split given the 11 statistics
#' plot(est)
#' 
#' # Formation plot with plots joined, for time steps 400 to 500
#' plot(est, plots.joined=TRUE, dx.start=400, dx.end=500)
#' 
#' # Partnership duration plot
#' plot(est, type="duration")
#' 
#' # Truncate plot to start after time 200, given the age ramp up period
#' plot(est, type="duration", dx.start=200)
#' }
#' 
plot.epiNet.est <- function(x,     
                            type = "formation",
                            dx.start,
                            dx.end,
                            dx.leg = TRUE,
                            plots.joined,
                            ...) {
  
  if (is.null(x$sim.stats)) {
    if (class(x) == "epiNet.est") {
      stop("\nepiNet.est object missing diagnostics. Rerun epiNet.est function with stats=TRUE")
    }
  }
  
  stats <- x$sim.stats
  
  if (missing(dx.start)) dx.start <- 1
  if (missing(dx.end)) dx.end <- nrow(stats)
  if (dx.end > nrow(stats)) 
    stop("epiNet.est object diagnostics only extend to time ", nrow(stats))
  
  if (type == "formation") {
    
    if (missing(plots.joined)) plots.joined <- ifelse(ncol(stats) > 5, FALSE, TRUE)
    
    plot.stats <- data.frame(stats[dx.start:dx.end, ])
    names(plot.stats) <- attributes(stats)$dimnames[[2]]
    plot.stats.mn <- data.frame(names=names(plot.stats),
                                stats=colMeans(plot.stats))
    plot.stats.mn$sorder <- 1:nrow(plot.stats.mn)
    stats.table <- merge(x$target.stats, plot.stats.mn, all=TRUE)
    stats.table$names <- as.character(stats.table$names)
    stats.table <- stats.table[do.call("order",
                                       stats.table[, "sorder", drop = FALSE]), , drop = FALSE]
    targs <- which(!is.na(stats.table$targets))
    sts <- which(!is.na(stats.table$stats))
    nstats <- ncol(plot.stats)
    
    if (plots.joined == TRUE) {
      
      pal <- brewer.ramp(nrow(stats.table), "Dark2")
      plot(1,1, xlim=c(dx.start, dx.end), ylim=c(min(plot.stats), max(plot.stats)), 
           type="n", xlab="Timestep", ylab="Statistic", ...)
      for (i in sts) {
        lines(x=dx.start:dx.end, y=plot.stats[,stats.table$names[i]], col=pal[i], lwd=2)
        if (any(targs == i)) {
          abline(h=stats.table$targets[i], lty=2, lwd=1, col=pal[i])
        }
      }
      if (dx.leg == TRUE) {
        legend("topleft", legend=attributes(stats)$dimnames[[2]], lwd=3, 
               col=pal[1:nstats], cex=0.75, bg="white") 
      }
    }
    
    if (plots.joined == FALSE) {
      
      if (nstats == 1) dimens <- c(1,1)
      if (nstats == 2) dimens <- c(1,2)
      if (nstats == 3) dimens <- c(1,3)
      if (nstats == 4) dimens <- c(2,2)
      if (nstats == 5) dimens <- c(2,3)
      if (nstats == 6) dimens <- c(2,3)
      if (nstats %in% 7:9) dimens <- c(3,3)
      if (nstats %in% 10:12) dimens <- c(4,3)
      if (nstats %in% 13:16) dimens <- c(4,4)
      if (nstats > 16) dimens <- rep(ceiling(sqrt(nstats)),2)
      
      ops <- list(mar=par()$mar, mfrow=par()$mfrow, mgp=par()$mgp)
      par(mar=c(2.5,2.5,2,1), mgp=c(2,1,0), mfrow=dimens)
      
      for (i in sts) {
        plot(x=dx.start:dx.end, y=plot.stats[,stats.table$names[i]], 
             xlim=c(dx.start, dx.end), type="l", 
             lwd=2, col="darkgrey", main=stats.table$names[i],
             xlab="", ylab="", ...)
        abline(h=stats.table$targets[i], lty=2, lwd=1)
      }
      
      # Reset graphical parameters
      on.exit(par(ops))
    } 
  }
  
  if (type == "duration") {
    
    # Edges only dissolution model
    if (x$dissolution == ~offset(edges)) {
      
      pages <- x$pages
      plot(x=dx.start:dx.end, y=pages[dx.start:dx.end], type="l", lwd=2, 
           col="darkgrey", xlab="Timestep", ylab="Partnership Age", ...)
      duration.expected <- exp(x$coef.diss[1]) + 1
      abline(h=duration.expected, lty=2, lwd=1.5)
      
    }
    
    # Edges + nodematch dissolution model
    if (length(all.names(x$dissolution)) > 3 && 
          any(all.names(x$dissolution) == "nodematch")) {
      
      pages <- x$pages
      pal <- c("steelblue", "firebrick")
      ymax <- max(c(max(pages[[1]], na.rm=TRUE)), max(pages[[2]], na.rm=TRUE))
      plot(x=dx.start:dx.end, y=pages$mean.page.match[dx.start:dx.end], 
           type="l", lwd=2, col=pal[1],
           ylim=c(0, ymax), xlab="Timestep", ylab="Partnership Age", ...)
      abline(h=(exp(x$coef.diss[1] + x$coef.diss[2]) + 1), col=pal[1], lty=2, lwd=1.5)
      lines(x=dx.start:dx.end, y=pages$mean.page.nomatch[dx.start:dx.end], 
            lwd=2, col=pal[2])
      abline(h=(exp(x$coef.diss[1]) + 1), col=pal[2], lty=2, lwd=1.5)
      legend("bottomright", legend=c("match", "no match"), lty=1, 
             col=pal, lwd=3, cex=0.75, bty="n")
    }
    
  }
  
}



#' @title Plot Diagnostics from an epiNet.simNet Object
#' 
#' @description This function plots values from diagnostic simulations
#'  in \code{epiNet.simNet}.
#'
#' @param x an \code{EpiModel} object of class \code{\link{epiNet.simNet}}.
#' @param type plot type, with options either of \code{type="formation"} for 
#'   partnership formation statistics or other network statistics specified in 
#'   \code{epiNet.simNet}, or \code{type="duration"} to plot the mean ages of 
#'   partnerships over time.
#' @param sim network simulation number to plot, with default to plot all 
#'   simulations in \code{epiNet.simNet} object.
#' @param dx.start start time for diagnostic plots. This must be a positive 
#'   integer.
#' @param dx.end end time for diagnostic plots. This must be less than or equal 
#'   to the number of time steps simulated in the \code{epiNet.est} diagnostics.
#' @param dx.leg if \code{TRUE}, show legend (only if \code{plots.joined=TRUE})
#' @param plots.joined if \code{TRUE}, combine all target statistics in one plot,
#'   versus one plot per target statistic if \code{FALSE}.
#' @param alpha transparency level for lines, where 0 = transparent and 1 = opaque.
#' @param lwd line width for output values.
#' @param ... additional arguments to pass to main plot.
#' 
#' @details
#' The plot function for \code{epiNet.simNet} objects will generate plots of two
#' types of model diagnostic statistics that were run as part of the diagnostic
#' tools within that simulation function. The \code{formation} plot shows the 
#' summary statistics requested in \code{stats.form}, which defaults to those in 
#' the formation formula. The \code{duration} plot shows the average age of all 
#' partnerships at each time step, up until the maximum time step requested. 
#' This is estimated with the \code{\link{edgelist.meanage}} function.
#' 
#' The \code{plots.joined} argument will control whether the target statistics 
#' in the \code{formation} plot are joined in one plot or plotted separately. The 
#' default is based on the number of network statistics requested. The layout of 
#' the separate plots within the larger plot window is also based on the number 
#' of statistics.
#' 
#' Required for these plots is that the estimation diagnostics are run in 
#' \code{epiNet.simNet}. This happens by default, and is set with the 
#' \code{stats=TRUE} argument in that function.
#' 
#' @method plot epiNet.simNet
#' @export
#' 
#' @keywords plot
#' @seealso \code{\link{epiNet.simNet}}
#' 
#' @examples
#' \dontrun{
#' # See epiNet.est help to understand these steps
#' nw <- network.initialize(n = 500, directed = FALSE)  
#' nw <- set.vertex.attribute(nw, "race", value = rep(0:1, each = 250))
#' formation <- ~ edges + nodematch("race") + degree(0) + concurrent
#' target.stats <- c(225, 187, 180, 90)
#' dissolution <- ~ offset(edges)
#' coef.diss <- dissolution.coefs(dissolution, durations = 20)
#' 
#' # Model estimation using the direct STERGM fit, then simulation
#' est <- epiNet.est(
#'   nw, 
#'   formation, 
#'   dissolution, 
#'   target.stats, 
#'   coef.diss,
#'   edapprox = FALSE)
#' nwsims <- epiNet.simNet(est, nsteps = 250, nsims = 10)
#' 
#' # Default will plot all simulations for target statistics formation plot
#' plot(nwsims)
#' 
#' # If that is too busy, try adjusting the line width or transparency with alpha
#' plot(nwsims, lwd=0.5)
#' plot(nwsims, lwd=0.75, alpha=0.5)
#' 
#' # If still too busy, then split plots
#' plot(nwsims, plots.joined=FALSE)
#' 
#' # Or plot just a subset of simulations
#' plot(nwsims, sim=1)
#' plot(nwsims, sim=4:6, lwd=0.6)
#' 
#' # Or plot just a subset of time
#' plot(nwsims, dx.start=100, dx.end=150, alpha=0.4)
#' 
#' # Duration plot shows the average age of partnerships along each time step
#' # Since a nodematch term was used in the dissolution, it will plot 
#' #   partnerships matched and unmatched separately
#' plot(nwsims, type="duration", alpha=0.5)
#' 
#' # Truncate plot to start after time 100, given the age ramp up period
#' plot(nwsims, type="duration", dx.start=100, lwd=0.5)
#' }
#' 
plot.epiNet.simNet <- function(x,  
                               type = "formation",
                               sim,
                               dx.start,
                               dx.end,
                               dx.leg = TRUE,
                               plots.joined,
                               alpha,
                               lwd,
                               ...) {
  
  # General warnings
  if (is.null(x$stats)) {
    stop("\nepiNet.simNet object missing diagnostics. Rerun epiNet.simNet function with stats=TRUE")
  }
  
  # Pull existing graphical parameters
  ops <- list(mar=par()$mar, mfrow=par()$mfrow, mgp=par()$mgp)
  
  # Pull data from object
  stats <- x$stats
  nts <- nrow(stats[[1]])
  if (missing(sim)) sim <- 1:x$nsims
  if (missing(dx.start)) dx.start <- 1
  if (missing(dx.end)) dx.end <- nts
  if (dx.end > nts) 
    stop("epiNet.simNet object diagnostics only extend to time ", nts)
  if (max(sim) > x$nsims) stop("Set sim to between 1 and ", x$nsims)
  
  # Set default alpha and lwd
  if (missing(alpha)) {
    if (x$nsims > 1) {
      alpha <- 0.8
    } else {
      alpha <- 1
    }
  }
  if (missing(lwd)) {
    if (x$nsims > 1) {
      lwd <- 1
    } else {
      lwd <- 2
    }
  }
  
  ## Formation target statistics plot ##
  if (type == "formation") {
    
    # Calculate formation stats summaries and targets to plot against
    plot.stats <- list()
    for (ns in 1:x$nsims) {
      plot.stats[[ns]] <- data.frame(stats[[ns]][dx.start:dx.end, ])
      names(plot.stats[[ns]]) <- attributes(stats[[1]])$names
    }
    
    plot.stats.mn <- data.frame(names=names(plot.stats[[1]]),
                                stats=colMeans(plot.stats[[1]]))
    plot.stats.mn$sorder <- 1:nrow(plot.stats.mn)
    stats.table <- merge(x$target.stats, plot.stats.mn, all=TRUE)
    stats.table$names <- as.character(stats.table$names)
    stats.table <- stats.table[do.call("order",
                                       stats.table[, "sorder", drop = FALSE]), , drop = FALSE]
    targs <- which(!is.na(stats.table$targets))
    sts <- which(!is.na(stats.table$stats))
    nstats <- ncol(plot.stats[[1]])
    
    if (missing(plots.joined)) {
      plots.joined <- ifelse(ncol(stats[[1]]) > 5, FALSE, TRUE)
    }
    if (plots.joined == TRUE) {
      
      if (nstats <= 9) {
        pal <- transco(brewer.pal(9, "Set1"), alpha)
      } else {
        pal <- transco(brewer.ramp(nrow(stats.table), "Dark2"), alpha)
      }
      
      ymin <- min(sapply(plot.stats, min))
      ymax <- max(sapply(plot.stats, max))
      plot(1,1, type="n",
           xlim=c(dx.start, dx.end), ylim=c(ymin, ymax), 
           xlab="", ylab="", ...)
      for (i in sts) {
        for (ns in sim) {
          lines(x=dx.start:dx.end, y=plot.stats[[ns]][,stats.table$names[i]], 
                col=pal[i], lwd=lwd)
        }
        if (any(targs == i)) {
          abline(h=stats.table$targets[i], lty=2, lwd=2.5, col=pal[i])
        }
      }
      if (dx.leg == TRUE) {
        legend("topleft", legend=attributes(stats[[1]])$names, lwd=3, 
               col=pal[1:nstats], cex=0.75, bg="white") 
      }
    }
    if (plots.joined == FALSE) {
      
      if (nstats == 1) dimens <- c(1,1)
      if (nstats == 2) dimens <- c(1,2)
      if (nstats == 3) dimens <- c(1,3)
      if (nstats == 4) dimens <- c(2,2)
      if (nstats == 5) dimens <- c(2,3)
      if (nstats == 6) dimens <- c(2,3)
      if (nstats %in% 7:9) dimens <- c(3,3)
      if (nstats %in% 10:12) dimens <- c(4,3)
      if (nstats %in% 13:16) dimens <- c(4,4)
      if (nstats > 16) dimens <- rep(ceiling(sqrt(nstats)),2)
      
      par(mar=c(2.5,2.5,2,1), mgp=c(2,1,0), mfrow=dimens)
      
      if (x$nsims > 1) {
        pal <- transco(brewer.ramp(x$nsims, "Set1"), alpha)
      } else {
        pal <- transco("darkgrey", alpha)
      }
      
      for (i in sts) {
        all <- do.call("cbind", plot.stats)
        stsvars <- which(names(all) == stats.table$names[i])
        allsts <- all[, stsvars]
        if (class(allsts) == "numeric") allsts <- data.frame(allsts)
        ymin <- min(allsts)
        ymax <- max(allsts)
        plot(1, 1, type="n", 
             xlim=c(dx.start, dx.end), ylim=c(ymin,ymax),
             main=stats.table$names[i],
             xlab="", ylab="", ...)
        for (ns in sim) {
          lines(x=dx.start:dx.end, y=allsts[,ns], 
                col=pal[ns], lwd=lwd)
        }
        abline(h=stats.table$targets[i], lty=2.5, lwd=1)
      } 
    }
  }
  
  ## Duration/Age Plot ##
  if (type == "duration") {
    
    # Edges only dissolution model
    if (x$dissolution == ~ offset(edges)) {
      
      if (x$nsims > 1) {
        pal <- transco(brewer.ramp(x$nsims, "Set1"), alpha)
      } else {
        pal <- transco("darkgrey", alpha)
      }
      
      pages <- x$pages
      ymax <- max(sapply(pages, max))
      plot(1,1, type="n", 
           xlab="Timestep", ylab="Partnership Age",
           xlim=c(dx.start, dx.end), ylim=c(0, ymax), ...)
      for (i in sim) {
        lines(x=dx.start:dx.end, y=pages[[i]][dx.start:dx.end], 
              col=pal[i], lwd=lwd)
      }
      abline(h=(exp(x$coef.diss[1]) + 1), lty=2, lwd=1.5)
    }
    
    # Edges + nodematch dissolution model
    if (length(all.names(x$dissolution)) > 3 && any(all.names(x$dissolution) == "nodematch")) {
      
      pages <- x$pages
      
      if (x$nsims > 1) {
        bpal <- transco(brewer.ramp(x$nsims, "Blues"), alpha)
        rpal <- transco(brewer.ramp(x$nsims, "Reds"), alpha)
      } else {
        bpal <- transco("steelblue", alpha)
        rpal <- transco("firebrick", alpha)
      }
      
      allmax <- vector()
      for (i in 1:length(pages)) allmax[i] <- max(sapply(pages[[i]], max))
      ymax <- max(allmax)
      plot(1,1, type="n", 
           xlab="Timestep", ylab="Partnership Age",
           xlim=c(dx.start, dx.end), ylim=c(0, ymax), ...)
      for (i in sim) {
        lines(x=dx.start:dx.end, y=pages[[i]][[1]][dx.start:dx.end], lwd=lwd, col=bpal[i])
        lines(x=dx.start:dx.end, y=pages[[i]][[2]][dx.start:dx.end], lwd=lwd, col=rpal[i])
      }
      abline(h=(exp(x$coef.diss[1] + x$coef.diss[2]) + 1), col=bpal[1], lty=2, lwd=2.5)
      abline(h=(exp(x$coef.diss[1]) + 1), col=rpal[1], lty=2, lwd=2.5)
      legend("bottomright", legend=c("match", "no match"), lty=1, 
             col=pal, lwd=3, cex=0.75, bty="n")
    }
    
  }
  
  # Reset graphical parameters
  on.exit(par(ops))
}


#' @title Plot Simulations from a Stochastic Network Epidemic Model
#'
#' @description This function plots three types of output from a stochastic
#' network model simulated through \code{\link{epiNet.simTrans}}.
#'
#' @param x an \code{EpiModel} model object of class \code{epiNet.simTrans}.
#' @param type type of plot: \code{type="sim"} for epidemic model results, 
#'   \code{type="network"} for a static network plot (\code{plot.network}),
#'   or \code{type="formation"} for network formation statistics.
#' @param sim if \code{type="network"}, simulation number for network graph.
#' @param at if \code{type="network"}, time step for network graph.
#' @param col.inf if \code{TRUE} and \code{type="network"}, automatic disease 
#'   status colors (blue = susceptible, red = infected, , green = recovered).
#' @param shp.bip if \code{type="network"}, specify shapes in bipartite network 
#'   for the second mode vertices, with acceptable inputs of "triangle" and 
#'   "square".
#' @param zeromarg if \code{TRUE} and \code{type="network"}, automatically 
#'   sets plot margins to 0 on all sides.
#' @param alpha if \code{type="formation"}, transparency level for lines, 
#'   where 0 = transparent and 1 = opaque.
#' @param lwd if \code{type="formation"}, line width for output values.
#' @param plots.joined if \code{TRUE} and \code{type="formation"}, combine all 
#'   target statistics in one plot, versus one plot per target statistic if \code{FALSE}.
#' @param ... additional arguments to pass to either plot type.
#' 
#' @details
#' This generic plot function can produce three types of plots given a stochastic
#' network model simulated through \code{\link{epiNet.simTrans}}:
#' \enumerate{
#'  \item If \code{type="sim"}, epidemic model results (e.g., disease prevalence
#'    and incidence) may be plotted. In this case, this plotting function wraps the
#'    \code{\link{plot.epiICM}} function, as the stochastic epidemiological results
#'    are in the same data structure. Consult the help page for \code{plot.epiICM}
#'    for all the plotting parameters.
#'  \item If \code{type="network"}, a static network plot will be generated. A 
#'    static network plot of a dynamic network of the sort simulated in 
#'    \code{EpiModel} is a cross-sectional extraction of that dynamic network at
#'    a specific time point. This plotting function wraps the 
#'    \code{\link{plot.network}} function in the \code{network} package. Consult 
#'    the help page for \code{plot.network} for all the plotting parameters. In 
#'    addition, five plotting parameters specific to \code{EpiModel} plots are 
#'    available: \code{sim}, \code{at}, \code{col.inf}, \code{shp.bip}, and 
#'    \code{zeromarg}.
#'  \item If \code{type="formation"}, summary network statistics related to partnership
#'    formation will be plotted. These formation plots are similar to the formation
#'    plots for both \code{\link{epiNet.est}} and \code{\link{epiNet.simNet}} objects.
#'    When running a \code{epiNet.simTrans} simulation, one must specify there that 
#'    \code{save.stats=TRUE}; the plot here will then show the network statistics
#'    requested explicitly in \code{stat.formula}, or use the formation formula
#'    set in \code{epiNet.est} otherwise. Note, that these network statistics are
#'    not saved for independent \code{epiNet.simTrans} simulations, since all of 
#'    the network data, including these statistics, is simulated and saved in the
#'    call to \code{epiNet.simNet}; one should plot the \code{epiNet.simNet} object
#'    in that situation.
#' }
#' 
#' @method plot epiNet.simTrans
#' @export 
#' 
#' @keywords plot
#' @seealso \code{\link{plot.epiICM}}, \code{\link{plot.network}}
#' 
#' @examples
#' ## See EpiModel Tutorial vignette for examples ##
#' 
plot.epiNet.simTrans <- function(x, 
                                 type = "sim", 
                                 sim, 
                                 at = 1, 
                                 col.inf = FALSE, 
                                 shp.bip = NULL,
                                 zeromarg = TRUE,
                                 alpha,
                                 lwd,
                                 plots.joined,
                                 ...) {
  
  if (type == "network") {
    
    if (at > max(x$time)) 
      stop("Specify a time step between 1 and ", max(x$time), sep="")
    
    if (missing(sim)) sim <- 1
    if (sim > x$nsims)
      stop("Specify sim between 1 and ", x$nsims)
    obj <- network.extract(x$network[[sim]], at=at)
    tea <- ifelse(any(names(obj$val[[1]]) %in% "status.active"), TRUE, FALSE)
    
    ops <- list(mar=par()$mar)
    if (zeromarg == TRUE) par(mar=c(0,0,0,0))
    
    if (!is.null(shp.bip)) {
      if (all(shp.bip != c("square", "triangle")))
        stop("shp.bip accepts inputs of either square or triangle")
      
      if (is.numeric(obj$gal$bipartite)) {
        mids <- idmode(obj)
        if (shp.bip == "square") {
          vertex.sides <- ifelse(mids == 1, 50, 4)
          vertex.rot <- 45
        }
        if (shp.bip == "triangle") {
          vertex.sides <- ifelse(mids == 1, 50, 3)
          vertex.rot <- 90
        }
        
      } else {
        warning("shp.bip only applies to bipartite networks, so ignoring")
        vertex.sides <- 50
        vertex.rot <- 0
      }
    } else {
      vertex.sides <- 50
      vertex.rot <- 0
    }
    if (col.inf == TRUE) {
      pal <- transco(c("firebrick", "steelblue", "seagreen"), 0.75)
      if (tea == TRUE) {
        cols <- ifelse(get.vertex.attribute.active(obj, "status", at=at) == 1,
                       pal[1], pal[2])
        cols <- ifelse(get.vertex.attribute.active(obj, "status", at=at) == 2,
                       pal[3], cols)
      } else {
        if (is.null(x$stat.mat)) {
          stop("Network plots of epiNet.simTrans objects require either tea=TRUE or save.statmat=TRUE")
        } else {
          if (is.numeric(obj$gal$bipartite)) {
            status <- x$stat.mat[[sim]][at, node.active(nw=x$network[[sim]], at=at, out="vec")]
          } else {
            status <- x$stat.mat[[sim]][at,]
          }
          cols <- ifelse(status == 1, pal[1], pal[2])
          cols <- ifelse(status == 2, pal[3], cols) 
        }
      }
      plot.network(obj, 
                   vertex.col = cols, 
                   vertex.border = "lightgrey", 
                   edge.col = "darkgrey", 
                   vertex.sides = vertex.sides,
                   vertex.rot = vertex.rot,
                   ...)
    } else {
      plot.network(obj, 
                   vertex.sides = vertex.sides, 
                   vertex.rot = vertex.rot,
                   ...)
    }
    
    on.exit(par(ops))
  }
  
  if (type == "sim") {
    plot.epiICM(x, ...)
  }
  
  ## Formation statistics plot ##
  if (type == "formation") {
    
    if (x$vital == FALSE)
      stop('For formation plot when vital=FALSE, plot the epiNet.simNet object')
    if (any(names(x) == 'stats') == FALSE)
      stop('For formation plot, set epiNet.simTrans parameter save.stats=TRUE')
    
    # Pull existing graphical parameters
    ops <- list(mar=par()$mar, mfrow=par()$mfrow, mgp=par()$mgp)
    
    # Set default alpha and lwd
    if (missing(alpha)) {
      if (x$nsims > 1) {
        alpha <- 0.8
      } else {
        alpha <- 1
      }
    }
    if (missing(lwd)) {
      if (x$nsims > 1) {
        lwd <- 1
      } else {
        lwd <- 2
      }
    }
    
    stats <- x$stats
    
    # Calculate formation stats summaries and targets to plot against
    plot.stats <- list()
    for (ns in 1:x$nsims) {
      plot.stats[[ns]] <- data.frame(stats[[ns]])
      names(plot.stats[[ns]]) <- attributes(stats[[1]])$names
    }
    
    plot.stats.mn <- data.frame(names=names(plot.stats[[1]]),
                                stats=colMeans(plot.stats[[1]]))
    nstats <- ncol(plot.stats[[1]])
    
    if (missing(plots.joined)) {
      plots.joined <- ifelse(ncol(stats[[1]]) > 5, FALSE, TRUE)
    }
    if (plots.joined == TRUE) {
      
      if (nstats <= 9) {
        pal <- transco(brewer.pal(9, "Set1"), alpha)
      } else {
        pal <- transco(brewer.ramp(nstats, "Dark2"), alpha)
      }
      
      ymin <- min(sapply(plot.stats, min))
      ymax <- max(sapply(plot.stats, max))
      plot(1,1, type="n",
           xlim=c(1, max(x$time)), ylim=c(ymin, ymax), 
           xlab="", ylab="", ...)
      for (i in 1:nstats) {
        for (j in 1:x$nsims) {
          lines(plot.stats[[j]][,i], col=pal[i], lwd=lwd)
        }
      }
      legend("topleft", legend=attributes(stats[[1]])$names, lwd=3, 
               col=pal[1:nstats], cex=0.75, bg="white") 
    }
    if (plots.joined == FALSE) {
      
      if (nstats == 1) dimens <- c(1,1)
      if (nstats == 2) dimens <- c(1,2)
      if (nstats == 3) dimens <- c(1,3)
      if (nstats == 4) dimens <- c(2,2)
      if (nstats == 5) dimens <- c(2,3)
      if (nstats == 6) dimens <- c(2,3)
      if (nstats %in% 7:9) dimens <- c(3,3)
      if (nstats %in% 10:12) dimens <- c(4,3)
      if (nstats %in% 13:16) dimens <- c(4,4)
      if (nstats > 16) dimens <- rep(ceiling(sqrt(nstats)),2)
      
      par(mar=c(2.5,2.5,2,1), mgp=c(2,1,0), mfrow=dimens)
      
      if (x$nsims > 1) {
        pal <- transco(brewer.ramp(x$nsims, "Set1"), alpha)
      } else {
        pal <- transco("darkgrey", alpha)
      }
      
      for (i in 1:nstats) {
        all <- do.call("cbind", plot.stats)
        stsvars <- which(names(all) == names(plot.stats[[1]])[i])
        allsts <- all[, stsvars, drop = FALSE]
        ymin <- min(allsts)
        ymax <- max(allsts)
        plot(1, 1, type="n", 
             xlim=c(1, max(x$time)), ylim=c(ymin,ymax),
             main=names(allsts)[1],
             xlab="", ylab="", ...)
        for (j in 1:x$nsims) {
          lines(allsts[, j], col=pal[j], lwd=lwd)
        }
      } 
    }
  
    # Reset graphical parameters
    on.exit(par(ops))
  }
  
  
}


#' @title Plot Compartment Diagram from Deterministic or Stochastic Model
#' 
#' @description This function plots a compartment diagram for all three 
#'   classes of epidemic models in EpiModel: \code{epiDCM}, \code{epiICM}, 
#'   and \code{epiNet.simTrans}.
#' 
#' @param x an \code{EpiModel} object of class \code{epiDCM}, \code{epiICM}, or 
#'   \code{epiNet.simTrans}
#' @param time timestep at which summary data is extracted
#' @param run model run number for plotting, only for \code{epiDCM} class
#'   models with multiple runs (sensitivity analyses)
#' @param digits number of significant digits to print
#' @param ... additional arguments passed to plot (not currently used)
#' 
#' @details
#' The \code{comp.plot} function is a visual summary of an EpiModel at a
#' specific point in time. The information contained in \code{comp.plot} 
#' is essentially the same as in the \code{summary} functions for a model, 
#' but presented in a way commonly seen in the epidemiological literature.
#'   
#' For \code{epiDCM} class plots, one must specify the specific run number to 
#' be plotted, in case the object contains multiple runs from a sensitivity
#' analysis. For \code{epiICM} and \code{epiNet.simTrans} class plots, the 
#' \code{run} argument is not used because the plots show the means and 
#' standard deviations for each state size and flow at the specified timestep.
#'   
#' The functionality of these plots is currently limited to one-group or one-
#' mode models in each of the three model classes, and that functionality 
#' will be expanded in future EpiModel releases.
#' 
#' @export
#' @keywords plot
#' 
#' @examples
#' # Example for epiDCM: SIR model with varying act.rate
#' mod <- epiDCM(type="SIR", s.num = 1000, i.num = 1, 
#'               trans.rate = 0.2, act.rate = 5:7, rec.rate = 1/3,
#'               b.rate = 1/90, ds.rate = 1/100, di.rate = 1/35, 
#'               dr.rate = 1/100, nsteps = 25)
#' comp.plot(mod, time=25, run=3)
#' 
#' # Example for epiICM: SIR model with 3 simulations
#' modICM <- epiICM(type = "SIR", s.num = 500, i.num = 1, 
#'                  trans.rate = 0.2, act.rate = 3,
#'                  rec.rate = 1/50, b.rate = 1/100,
#'                  ds.rate = 1/100, di.rate = 1/90, dr.rate = 1/100,
#'                  nsteps = 25, nsims = 3)
#' comp.plot(modICM, time=25)
#' 
comp.plot <- function(x, time, run, digits, ...) {
  UseMethod("comp.plot")
}


#' @method comp.plot epiDCM
#' @rdname comp.plot
#' @export
comp.plot.epiDCM <- function(x, 
                             time,
                             run = 1, 
                             digits = 3,
                             ...
                             ) {
  
  
  ## Number of model runs
  nruns <- x$nruns
  
  ## Errors
  if (x$groups != 1)
    stop("Only 1-group epiDCM models currently supported")
  if (run > nruns) 
    stop("Specify a run between 1 and ", nruns)
  
  ## Time
  if (missing(time) || (time > max(x$time) | time < min(x$time))) 
    stop("Specify a timestep between 1 and ", max(x$time))
  intime <- time
  time <- which(x$time == intime)
  
  ## Dataframe subsets
  df <- as.data.frame(x, run=run)
  df <- round(df[time, ], digits)
  
  ## Change graphical parameters
  ops <- list(mar=par()$mar, mfrow=par()$mfrow, mgp=par()$mgp)
  par(mar=c(0,0,2,0))
  options(scipen = 10)
  
  ## Main Plot
  plot(0:100, 0:100, type="n", axes=FALSE)
  title(main=paste(x$type, "Model Diagram"))
  mtext(paste("time=", intime, "  |  run=", run, sep=""), side=3, cex=0.8, line=-1)
  
  ## 1. SI Model
  if (x$type == "SI" && x$groups == 1) {
    mbox(22, 40, "Susceptible", df$s.num)
    mbox(57, 40, "Infected", df$i.num)
    harrow(22, 40, "si.flow", df$si.flow, dir="right")
    if (x$vital == TRUE) {
      varrow(22, 40, "ds.flow", df$ds.flow, dir="out")
      varrow(57, 40, "di.flow", df$di.flow, dir="out")
      varrow(22, 40, "b.flow", df$b.flow, dir="in")
    }
  }
  
  ## 2. SIR Model
  if (x$type == "SIR" && x$groups == 1) {
    mbox(5, 40, "Susceptible", df$s.num)
    mbox(40, 40, "Infected", df$i.num)
    mbox(75, 40, "Recovered", df$r.num)
    harrow(5, 40, "si.flow", df$si.flow, dir="right")   
    harrow(40, 40, "ir.flow", df$ir.flow, dir="right")
    if (x$vital == TRUE) {
      varrow(5, 40, "ds.flow", df$ds.flow, dir="out")
      varrow(40, 40, "di.flow", df$di.flow, dir="out")
      varrow(75, 40, "dr.flow", df$dr.flow, dir="out")
      varrow(5, 40, "b.flow", df$b.flow, dir="in")
    }
  } 
  
  ## 3. SIS Model
  if (x$type == "SIS" && x$groups == 1) {
    mbox(22, 40, "Susceptible", df$s.num)
    mbox(57, 40, "Infected", df$i.num)    
    harrow(22, 40, "si.flow", df$si.flow, dir="right")    
    harrow(22, 40, "is.flow", df$is.flow, dir="left")
    if (x$vital == TRUE) {
      varrow(22, 40, "ds.flow", df$ds.flow, dir="out")
      varrow(57, 40, "di.flow", df$di.flow, dir="out")
      varrow(22, 40, "b.flow", df$b.flow, dir="in")
    }
  }
  
  # Reset graphical parameters
  on.exit(par(ops))
}


#' @method comp.plot epiICM
#' @rdname comp.plot
#' @export
comp.plot.epiICM <- function(x, 
                             time = 1,
                             run,
                             digits = 3,
                             ...
                             ) {
  
  # Number of model runs
  nsims <- x$nsims
  
  # Errors
  if (class(x) == "epiICM") groups <- x$groups
  if (class(x) == "epiNet.simTrans") groups <- x$modes
  if (groups != 1) stop("Only 1-group/mode models currently supported")
  
  # Time
  if (missing(time) || (time > max(x$time) | time < min(x$time))) 
    stop("Specify a timestep between 1 and ", max(x$time))
  
  ## Dataframe subsets for plots
  df.mn <- as.data.frame(x, out="mean")
  df.mn <- round(df.mn[time == df.mn$time, ], digits)
  df.sd <- as.data.frame(x, out="sd")
  df.sd <- round(df.sd[time == df.sd$time, ], digits)
  
  ## Change graphical parameters
  ops <- list(mar=par()$mar, mfrow=par()$mfrow, mgp=par()$mgp)
  par(mar=c(0,0,2,0))
  options(scipen = 10)
  
  ## Main Plot
  plot(0:100, 0:100, type="n", axes=FALSE)
  title(main=paste(x$type, "Model Diagram"))
  mtext(paste("Simulation means(sd) | time=", time, sep=""), side=3, cex=0.8, line=-1)
  
  ## 1. SI Model
  if (x$type == "SI" && groups == 1) {
    mbox(22, 40, "Susceptible", paste(df.mn$s.num, "(", df.sd$s.num, ")", sep=""))
    mbox(57, 40, "Infected", paste(df.mn$i.num, "(", df.sd$i.num, ")", sep=""))
    harrow(22, 40, "si.flow", df.mn$si.flow, dir="right")
    if (x$vital == TRUE) {
      varrow(22, 40, "ds.flow", df.mn$ds.flow, dir="out")
      varrow(57, 40, "di.flow", df.mn$di.flow, dir="out")
      varrow(22, 40, "b.flow", df.mn$b.flow, dir="in")
    }
  }
  
  ## 2. SIR Model
  if (x$type == "SIR" && groups == 1) {
    mbox(5, 40, "Susceptible", paste(df.mn$s.num, "(", df.sd$s.num, ")", sep=""))
    mbox(40, 40, "Infected", paste(df.mn$i.num, "(", df.sd$i.num, ")", sep=""))
    mbox(75, 40, "Recovered", paste(df.mn$r.num, "(", df.sd$r.num, ")", sep=""))
    harrow(5, 40, "si.flow", df.mn$si.flow, dir="right")   
    harrow(40, 40, "ir.flow", df.mn$ir.flow, dir="right")
    if (x$vital == TRUE) {
      varrow(5, 40, "ds.flow", df.mn$ds.flow, dir="out")
      varrow(40, 40, "di.flow", df.mn$di.flow, dir="out")
      varrow(75, 40, "dr.flow", df.mn$dr.flow, dir="out")
      varrow(5, 40, "b.flow", df.mn$b.flow, dir="in")
    }
  } 
  
  ## 3. SIS Model
  if (x$type == "SIS" && groups == 1) {
    mbox(22, 40, "Susceptible", paste(df.mn$s.num, "(", df.sd$s.num, ")", sep=""))
    mbox(57, 40, "Infected", paste(df.mn$i.num, "(", df.sd$i.num, ")", sep=""))    
    harrow(22, 40, "si.flow", df.mn$si.flow, dir="right")    
    harrow(22, 40, "is.flow", df.mn$is.flow, dir="left")
    if (x$vital == TRUE) {
      varrow(22, 40, "ds.flow", df.mn$ds.flow, dir="out")
      varrow(57, 40, "di.flow", df.mn$di.flow, dir="out")
      varrow(22, 40, "b.flow", df.mn$b.flow, dir="in")
    }
  }
  
  # Reset graphical parameters
  on.exit(par(ops))
  
}

#' @method comp.plot epiNet.simTrans
#' @rdname comp.plot
#' @export
comp.plot.epiNet.simTrans <- function(x, 
                                      time = 1,
                                      run,
                                      digits = 3,
                                      ...
                                      ) {
  
  comp.plot.epiICM(x = x,
                   time = time,
                   run = run,
                   digits = digits,
                   ...)
  
}



# Helper Functions --------------------------------------------------------


# Calculate denominators
denom <- function(x, y, popfrac) {
  
  if (class(x) == "epiDCM") {
    if (popfrac == TRUE) {
      den <- data.frame(den=rep(NA, length(x$time)))
      if (x$groups == 1) {
        for (i in 1:x$nruns) {
          den[,i] <- as.data.frame(x, run=i)$num
        }
        for (i in 1:length(y)) {
          x[[y[i]]] <- x[[y[i]]] / den
        }
      }
      if (x$groups == 2) {
        den <- list(den.g1 = den, den.g2 = den)
        for (i in 1:x$nruns) {
          den[[1]][, i] <- as.data.frame(x, run=i)$num
          den[[2]][, i] <- as.data.frame(x, run=i)$num.g2
        }
        y.group.num <- ifelse(grepl("num.g2$", y), 2, 1)
        for (j in 1:length(y)) {
          x[[y[j]]] <- x[[y[j]]] / den[[y.group.num[j]]]
        }
      }  
    }
    if (popfrac == FALSE && x$nruns == 1) {
      for (j in 1:length(y)) {
        x[[y[j]]] <- data.frame(x[[y[j]]])
      }
    }
  }
  
  if (class(x) == "epiICM") {
    if (popfrac == TRUE) {
      den <- data.frame(den=rep(NA, length(x$time)))
      if (x$groups == 1) {
        for (i in 1:x$nsims) {
          den[,i] <- as.data.frame(x, sim=i, out="vals")$num
        }
        for (i in 1:length(y)) {
          x[[y[i]]] <- x[[y[i]]] / den
        }
      }
      if (x$groups == 2) {
        den <- list(den.g1 = den, den.g2 = den)
        for (i in 1:x$nsims) {
          den[[1]][, i] <- as.data.frame(x, sim=i, out="vals")$num
          den[[2]][, i] <- as.data.frame(x, sim=i, out="vals")$num.g2
        }
        y.group.num <- ifelse(grepl("num.g2$", y), 2, 1)
        for (j in 1:length(y)) {
          x[[y[j]]] <- x[[y[j]]] / den[[y.group.num[j]]]
        }
      } 
    }
    if (popfrac == FALSE && x$nsims == 1) {
      for (j in 1:length(y)) {
        x[[y[j]]] <- data.frame(x[[y[j]]])
      }
    }
  }
  
  if (class(x) == "epiNet.simTrans") {
    if (popfrac == TRUE) {
      den <- data.frame(den=rep(NA, length(x$time)))
      if (x$modes == 1) {
        for (i in 1:x$nsims) {
          den[,i] <- as.data.frame(x, sim=i, out="vals")$num
        }
        for (i in 1:length(y)) {
          x[[y[i]]] <- x[[y[i]]] / den
        }
      }
      if (x$modes == 2) {
        den <- list(den.g1 = den, den.g2 = den)
        for (i in 1:x$nsims) {
          den[[1]][, i] <- as.data.frame(x, sim=i, out="vals")$num
          den[[2]][, i] <- as.data.frame(x, sim=i, out="vals")$num.m2
        }
        y.group.num <- ifelse(grepl("num.m2$", y), 2, 1)
        for (j in 1:length(y)) {
          x[[y[j]]] <- x[[y[j]]] / den[[y.group.num[j]]]
        }
      }  
    }
    if (popfrac == FALSE && x$nsims == 1) {
      for (j in 1:length(y)) {
        x[[y[j]]] <- data.frame(x[[y[j]]])
      }
    }
  }
  
  return(x)
}

## comp.plot helper utilities ##
#  Text box
mbox <- function(x, y, title, val) {
  polygon(c(x, x+20, x+20, x), c(y, y, y+20, y+20))
  text(x+10, y+10, paste(title, "\n n=", val, sep=""), cex=0.9)
}
#  Horizontal arrow
harrow <- function(xbox, ybox, title, val, dir) {
  if (dir == "right") {
    arrows(xbox+20, ybox+12, xbox+35, lwd=2, length=0.15)
    text(xbox+27.5, ybox+17, paste(title, val, sep="="), cex=0.8) 
  }
  if (dir == "left") {
    arrows(xbox+20+15, ybox+5, xbox+20, lwd=2, length=0.15)
    text(xbox+27.5, ybox+2, paste(title, val, sep="="), cex=0.8)
  }
}
#  Vertical arrow
varrow <- function(xbox, ybox, title, val, dir) {
  if (dir == "out") {
    arrows(xbox+10, ybox, xbox+10, ybox-25, lwd=2, length=0.15)
    text(xbox+10, ybox-12.5, paste(title, val, sep="="), cex=0.8, pos=4)
  }
  if (dir == "in") {
    arrows(xbox+10, ybox+45, xbox+10, ybox+20, lwd=2, length=0.15)
    text(xbox+10, ybox+32.5, paste(title, val, sep="="), cex=0.8, pos=4)
  }
}
