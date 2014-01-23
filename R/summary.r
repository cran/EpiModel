
#' @title Summary Model Statistics
#'
#' @description This function extracts and prints model statistics solved 
#'  with the \code{epiDCM} function.
#'
#' @param object an \code{EpiModel} object of class \code{epiDCM}.
#' @param time time step for model statistics.
#' @param run run number for model, if sensitivity analyses were conducted
#'   through \code{epiDCM}.
#' @param digits number of significant digits to print.
#' @param comp.plot if \code{TRUE}, plot compartments and flows in summary. This
#'   can be called separately using the \code{\link{comp.plot}} function.
#' @param ... additional summary function arguments.
#' 
#' @details
#' Summary statistics for the main epidemiological outcomes (state and transition
#' size and prevalence) from an \code{epiDCM} model. Time-specific summary measures
#' are provided, so it is necessary to input a time of interest. For multiple-run
#' models (e.g., sensitivity analyses), one may also input a model run number. See
#' the examples below.
#' 
#' @seealso \code{\link{epiDCM}}
#' 
#' @method summary epiDCM
#' @keywords extract
#' @export
#' 
#' @examples
#' \dontrun{
#' # Deterministic SIR model with varying act.rate
#' mod <- epiDCM(type="SIR", s.num=1000, i.num=1, r.num=0,
#'               trans.rate=0.2, act.rate=2:8, rec.rate=1/3,
#'               b.rate=0.011, ds.rate=0.01, di.rate=0.03, 
#'               dr.rate=0.01, nsteps=500, verbose=TRUE)
#' summary(mod, time=25, run=1, comp.plot=TRUE)
#' summary(mod, time=25, run=7, comp.plot=TRUE)
#' summary(mod, time=26, run=7, comp.plot=TRUE)
#' }
#' 
summary.epiDCM <- function(object, 
                           time, 
                           run = 1, 
                           digits = 3, 
                           comp.plot = FALSE,
                           ...) {
  
  nvars <- length(object)
  nruns <- object$nruns
  type <- object$type
  groups <- object$groups
  vital <- object$vital
  
  df <- as.data.frame(object, run=run)
  
  if (missing(time) || (time > max(object$time) | time < min(object$time))) 
    stop("Specify a time between 1 and ", max(object$time))
  df <- df[df$time == time, ]
  
  ## Some calculations
  df$s.prev <- df$s.num / df$num
  df$i.prev <- df$i.num / df$num
  if (type == "SIR") {
    df$r.prev <- df$r.num / df$num
  }
  if (groups == 2) {
    df$s.prev.g2 <- df$s.num.g2 / df$num.g2
    df$i.prev.g2 <- df$i.num.g2 / df$num.g2
    if (type == "SIR") {
      df$r.prev.g2 <- df$r.num.g2 / df$num.g2
    }
  }
  
  if (type == "SI") {
    stats <- with(df, c(s.num, s.prev, 
                        i.num, i.prev,
                        num, 1,
                        si.flow, NA))
    mat <- matrix(stats, byrow=T, nrow=length(stats)/2)
    rownames(mat) <- c("Suscept.","Infect.","Total","S -> I")
    if (vital == TRUE) {
      stats <- with(df, c(b.flow, NA,
                          ds.flow, NA,
                          di.flow, NA))
      mat <- rbind(mat, matrix(stats, byrow=T, nrow=length(stats)/2))
      rownames(mat)[rownames(mat)==""] <- c("Birth ->", "S Death ->", "I Death ->")
    }
    if (groups == 2) {
      stats <- with(df, c(s.num.g2, s.prev.g2, 
                          i.num.g2, i.prev.g2,
                          num.g2, 1,
                          si.flow.g2, NA))
      mat.g2 <- matrix(stats, byrow=T, nrow=length(stats)/2)
      if (vital == TRUE) {
        stats <- with(df, c(b.flow.g2, NA,
                            ds.flow.g2, NA,
                            di.flow.g2, NA))
        mat.g2 <- rbind(mat.g2, matrix(stats, byrow=T, nrow=length(stats)/2))
      }
      mat <- cbind(mat, mat.g2)
    }
  }
  if (type == "SIR") {
    stats <- with(df, c(s.num, s.prev, 
                        i.num, i.prev, 
                        r.num, r.prev,
                        num, 1,
                        si.flow, NA,
                        ir.flow, NA))
    mat <- matrix(stats, byrow=T, nrow=length(stats)/2)
    rownames(mat) <- c("Suscept.","Infect.","Recov.","Total",
                       "S -> I", "I -> R")
    if (vital == TRUE) {
      stats <- c(df$b.flow, NA,
                 df$ds.flow, NA,
                 df$di.flow, NA,
                 df$dr.flow, NA)
      mat <- rbind(mat, matrix(stats, byrow=T, nrow=length(stats)/2))
      rownames(mat)[rownames(mat)==""] <- c("Birth ->", "S Death ->", 
                                            "I Death ->", "R Death ->")
    }
    if (groups == 2) {
      stats <- with(df, c(s.num.g2, s.prev.g2, 
                          i.num.g2, i.prev.g2,
                          r.num.g2, r.prev.g2,
                          num.g2, 1,
                          si.flow.g2, NA,
                          ir.flow.g2, NA))
      mat.g2 <- matrix(stats, byrow=T, nrow=length(stats)/2)
      if (vital == TRUE) {
        stats <- with(df, c(b.flow.g2, NA,
                            ds.flow.g2, NA,
                            di.flow.g2, NA,
                            dr.flow.g2, NA))
        mat.g2 <- rbind(mat.g2, matrix(stats, byrow=T, nrow=length(stats)/2))
      }
      mat <- cbind(mat, mat.g2)
    }
  }
  if (type == "SIS") {
    stats <- with(df, c(s.num, s.prev, 
                        i.num, i.prev, 
                        num, 1,
                        si.flow, NA,
                        is.flow, NA))
    mat <- matrix(stats, byrow=T, nrow=length(stats)/2)
    rownames(mat) <- c("Suscept.","Infect.","Total",
                       "S -> I","I -> S") 
    if (vital == TRUE) {
      stats <- c(df$b.flow, NA,
                 df$ds.flow, NA,
                 df$di.flow, NA)
      mat <- rbind(mat, matrix(stats, byrow=T, nrow=length(stats)/2))
      rownames(mat)[rownames(mat)==""] <- c("Birth ->", "S Death ->", "I Death ->")
    }
    if (groups == 2) {
      stats <- with(df, c(s.num.g2, s.prev.g2, 
                          i.num.g2, i.prev.g2,
                          num.g2, 1,
                          si.flow.g2, NA,
                          is.flow.g2, NA))
      mat.g2 <- matrix(stats, byrow=T, nrow=length(stats)/2)
      if (vital == TRUE) {
        stats <- with(df, c(b.flow.g2, NA,
                            ds.flow.g2, NA,
                            di.flow.g2, NA))
        mat.g2 <- rbind(mat.g2, matrix(stats, byrow=T, nrow=length(stats)/2))
      }
      mat <- cbind(mat, mat.g2)
    }
  }
  
  if (groups == 1) colnames(mat) <- c("n","pct")
  if (groups == 2) colnames(mat) <- c("n:g1", "pct:g1", "n:g2", "pct:g2")
  
  mat <- round(mat, digits)
  
  # print it
  cat("EpiModel Summary")
  cat("\n=======================")
  cat("\nModel class:", class(object))
  
  cat("\n\nSimulation Summary")
  cat("\n-----------------------")
  cat("\nModel type:", type)
  cat("\nNo. runs:", nruns)
  cat("\nNo. time steps:", max(object$time))
  cat("\nNo. groups:", groups)
  
  if (groups == 1) statsep <- "---------------------------------"
  if (groups == 2) statsep <- "--------------------------------------------------------"
  
  cat("\n\nModel Statistics\n")
  cat(statsep)
  cat("\nTime:", time)
  cat("\t Run:", run, "\n")
  cat(statsep, "\n")
  print(mat, print.gap=4)
  cat(statsep, "\n")
  
  if (comp.plot == TRUE) {
    comp.plot(object, 
              time=time,
              run=run, 
              digits=digits)
  }
}


#' @title Summary Model Statistics
#'
#' @description This function extracts and prints model statistics 
#'   solved with the \code{epiICM} function.
#'
#' @param object an \code{EpiModel} object of class \code{epiICM}.
#' @param time time step for model statistics.
#' @param digits number of significant digits to print.
#' @param comp.plot if \code{TRUE}, plot compartments and flows in summary. This
#'  can be called separately using the \code{\link{comp.plot}} function.
#' @param ... additional summary function arguments.
#' 
#' @details
#' Summary statistics for the main epidemiological outcomes (state and transition
#' size and prevalence) from an \code{epiICM} model. Time-specific summary measures
#' are provided, so it is necessary to input a time of interest. One may 
#' simultaneously obtain console-based statistics and a compartment plot via
#' the \code{comp.plot} function through the \code{comp.plot} argument.
#' 
#' @seealso \code{\link{epiICM}}, \code{\link{comp.plot}}
#' 
#' @method summary epiICM
#' @keywords extract
#' @export
#' 
#' @examples
#' \dontrun{
#' # Stochastic SI model with 10 simulations
#' mod <- epiICM(type = "SI", s.num = 500, i.num = 1, 
#'               trans.rate = 0.2, act.rate = 0.25, 
#'               nsteps = 500, nsims = 10)
#' summary(mod, time = 25)
#' summary(mod, time = 50, comp.plot = TRUE)
#' }
#' 
summary.epiICM <- function(object, 
                           time, 
                           digits = 3, 
                           comp.plot = FALSE,
                           ...) {
  
  nvars <- length(object)
  nsims <- object$nsims
  type <- object$type
  groups <- object$groups
  vital <- object$vital
  nts <- max(object$time)
  
  if (missing(time) || (time > max(object$time) | time < 1))
    stop("Specify a timestep between 1 and ", max(object$time))
  
  df.mn <- as.data.frame(object, time=time)
    df.mn <- df.mn[df.mn$time == time, ]
  df.sd <- as.data.frame(object, time=time, out="sd")
    df.sd <- df.sd[df.sd$time == time, ]
  
  if (type == "SI") {
    
    ## Prevalence calcs
    s.prev <- df.mn$s.num/df.mn$num
    i.prev <- df.mn$i.num/df.mn$num
    if (groups == 2) {
      s.prev.g2 <- df.mn$s.num.g2/df.mn$num.g2
      i.prev.g2 <- df.mn$i.num.g2/df.mn$num.g2
    }
    
    ## Group 1 stats
    stats <- c(df.mn$s.num, df.sd$s.num, s.prev,
               df.mn$i.num, df.sd$i.num, i.prev,
               df.mn$num, df.sd$num, 1, 
               df.mn$si.flow, df.sd$si.flow, NA)
    mat <- matrix(stats, byrow=T, nrow=length(stats)/3)
    rownames(mat) <- c("Suscept.",
                       "Infect.",
                       "Total",
                       "S -> I")
    if (vital == TRUE) {
      stats <- c(df.mn$b.flow, df.sd$b.flow, NA,
                 df.mn$ds.flow, df.sd$ds.flow, NA,
                 df.mn$di.flow, df.sd$di.flow, NA)
      mat <- rbind(mat, matrix(stats, byrow=T, nrow=length(stats)/3))
      rownames(mat)[rownames(mat)==""] <- c("Birth ->", "S Death ->", "I Death ->")
    }
    
    ## Group 2 stats
    if (groups == 2) {
      stats <- c(df.mn$s.num.g2, df.sd$s.num.g2, s.prev.g2,
                 df.mn$i.num.g2, df.sd$i.num.g2, i.prev.g2,
                 df.mn$num.g2, df.sd$num.g2, 1, 
                 df.mn$si.flow.g2, df.sd$si.flow.g2, NA)
      mat.g2 <- matrix(stats, byrow=T, nrow=length(stats)/3)
      if (vital == TRUE) {
        stats <- c(df.mn$b.flow.g2, df.sd$b.flow.g2, NA,
                   df.mn$ds.flow.g2, df.sd$ds.flow.g2, NA,
                   df.mn$di.flow.g2, df.sd$di.flow.g2, NA)
        mat.g2 <- rbind(mat.g2, matrix(stats, byrow=T, nrow=length(stats)/3))
      }
      mat <- cbind(mat, mat.g2)
    }
    
  } # end SI summary
  
  
  if (type == "SIR") {
    
    ## Prevalence calcs
    s.prev <- df.mn$s.num/df.mn$num
    i.prev <- df.mn$i.num/df.mn$num
    r.prev <- df.mn$r.num/df.mn$num
    if (groups == 2) {
      s.prev.g2 <- df.mn$s.num.g2/df.mn$num.g2
      i.prev.g2 <- df.mn$i.num.g2/df.mn$num.g2
      r.prev.g2 <- df.mn$r.num.g2/df.mn$num.g2
    }
    
    ## Group 1 stats
    stats <- c(df.mn$s.num, df.sd$s.num, s.prev,
               df.mn$i.num, df.sd$i.num, i.prev,
               df.mn$r.num, df.sd$r.num, r.prev,
               df.mn$num, df.sd$num, 1,
               df.mn$si.flow, df.sd$si.flow, NA,
               df.mn$ir.flow, df.sd$ir.flow, NA)
    mat <- matrix(stats, byrow=T, nrow=length(stats)/3)
    rownames(mat) <- c("Suscept.",
                       "Infect.",
                       "Recov.",
                       "Total",
                       "S -> I",
                       "I -> R")
    if (vital == TRUE) {
      stats <- c(df.mn$b.flow, df.sd$b.flow, NA,
                 df.mn$ds.flow, df.sd$ds.flow, NA,
                 df.mn$di.flow, df.sd$di.flow, NA,
                 df.mn$dr.flow, df.sd$dr.flow, NA)
      mat <- rbind(mat, matrix(stats, byrow=T, nrow=length(stats)/3))
      rownames(mat)[rownames(mat)==""] <- c("Birth ->", 
                                            "S Death ->", 
                                            "I Death ->",
                                            "R Death ->")
    }
    
    ## Group 2 stats
    if (groups == 2) {
      stats <- c(df.mn$s.num.g2, df.sd$s.num.g2, s.prev.g2,
                 df.mn$i.num.g2, df.sd$i.num.g2, i.prev.g2,
                 df.mn$r.num.g2, df.sd$r.num.g2, r.prev.g2,
                 df.mn$num.g2, df.sd$num.g2, 1,
                 df.mn$si.flow.g2, df.sd$si.flow.g2, NA,
                 df.mn$ir.flow.g2, df.sd$ir.flow.g2, NA)
      mat.g2 <- matrix(stats, byrow=T, nrow=length(stats)/3)
      if (vital == TRUE) {
        stats <- c(df.mn$b.flow.g2, df.sd$b.flow.g2, NA,
                   df.mn$ds.flow.g2, df.sd$ds.flow.g2, NA,
                   df.mn$di.flow.g2, df.sd$di.flow.g2, NA,
                   df.mn$dr.flow.g2, df.sd$dr.flow.g2, NA)
        mat.g2 <- rbind(mat.g2, matrix(stats, byrow=T, nrow=length(stats)/3))
      }
      mat <- cbind(mat, mat.g2)
    }
  } # end SIR summary
  
  
  if (type == "SIS") {
    
    ## Prevalence calcs
    s.prev <- df.mn$s.num/df.mn$num
    i.prev <- df.mn$i.num/df.mn$num
    if (groups == 2) {
      s.prev.g2 <- df.mn$s.num.g2/df.mn$num.g2
      i.prev.g2 <- df.mn$i.num.g2/df.mn$num.g2
    }
    
    ## Group 1 stats
    stats <- c(df.mn$s.num, df.sd$s.num, s.prev,
               df.mn$i.num, df.sd$i.num, i.prev,
               df.mn$num, df.sd$num, 1,
               df.mn$si.flow, df.sd$si.flow, NA,
               df.mn$is.flow, df.sd$is.flow, NA)
    mat <- matrix(stats, byrow=T, nrow=length(stats)/3)
    rownames(mat) <- c("Suscept.",
                       "Infect.",
                       "Total",
                       "S -> I",
                       "I -> S")
    if (vital == TRUE) {
      stats <- c(df.mn$b.flow, df.sd$b.flow, NA,
                 df.mn$ds.flow, df.sd$ds.flow, NA,
                 df.mn$di.flow, df.sd$di.flow, NA)
      mat <- rbind(mat, matrix(stats, byrow=T, nrow=length(stats)/3))
      rownames(mat)[rownames(mat)==""] <- c("Birth ->", 
                                            "S Death ->", 
                                            "I Death ->")
    }
    
    ## Group 2 stats
    if (groups == 2) {
      stats <- c(df.mn$s.num.g2, df.sd$s.num.g2, s.prev.g2,
                 df.mn$i.num.g2, df.sd$i.num.g2, i.prev.g2,
                 df.mn$num.g2, df.sd$num.g2, 1,
                 df.mn$si.flow.g2, df.sd$si.flow.g2, NA,
                 df.mn$is.flow.g2, df.sd$is.flow.g2, NA)
      mat.g2 <- matrix(stats, byrow=T, nrow=length(stats)/3)
      if (vital == TRUE) {
        stats <- c(df.mn$b.flow.g2, df.sd$b.flow.g2, NA,
                   df.mn$ds.flow.g2, df.sd$ds.flow.g2, NA,
                   df.mn$di.flow.g2, df.sd$di.flow.g2, NA)
        mat.g2 <- rbind(mat.g2, matrix(stats, byrow=T, nrow=length(stats)/3))
      }
      mat <- cbind(mat, mat.g2)
    }
  } # end SIS summary
  
  if (groups == 1) colnames(mat) <- c("mean","sd", "perc")
  if (groups == 2) colnames(mat) <- c("mean:g1", "sd:g1", "perc:g1",
                                      "mean:g2", "sd:g2", "perc:g2")
  mat <- round(mat, digits)
  
  ## Print it
  cat("\nEpiModel Summary")
  cat("\n=======================")
  cat("\nModel class:", class(object))
  
  cat("\n\nSimulation Details")
  cat("\n-----------------------")
  cat("\nModel type:", type)
  cat("\nNo. simulations:", nsims)
  cat("\nNo. time steps:", nts)
  cat("\nNo. groups:", groups)
  
  if (groups == 1) statsep <- "---------------------------------"
  if (groups == 2) statsep <- "------------------------------------------------------------------------"
  
  cat("\n\nModel Statistics\n")
  cat(statsep)
  cat("\nTime:", time, "\n")
  cat(statsep, "\n")
  print(mat, print.gap=4)
  cat(statsep, "\n")
  
  if (comp.plot == TRUE) {
    comp.plot(object, 
              time=time,
              digits=digits)
  }  
}


#' @title Summary Model Statistics
#'
#' @description This function extracts and prints model statistics solved 
#'  with the \code{epiNet.simTrans} function.
#'
#' @param object an \code{EpiModel} object of class \code{epiNet.simTrans}.
#' @param time timestep for model statistics.
#' @param digits number of significant digits to print.
#' @param comp.plot if \code{TRUE}, plot compartments and flows in summary. This
#'  can be called separately using the \code{\link{comp.plot}} function.
#' @param ... additional summary function arguments.
#' 
#' @details
#' Summary statistics for the main epidemiological outcomes (state and transition
#' size and prevalence) from an \code{epiNet.simTrans} model. Time-specific 
#' summary measures are provided, so it is necessary to input a time of interest. 
#' One may simultaneously obtain console-based statistics and a compartment plot via
#' the \code{comp.plot} function through the \code{comp.plot} argument.
#' 
#' @seealso \code{\link{epiNet.simTrans}}, \code{\link{comp.plot}}
#' 
#' @method summary epiNet.simTrans
#' @keywords extract
#' @export
#' 
#' @examples
#' ## See EpiModel Tutorial vignette ##
#' 
summary.epiNet.simTrans <- function(object, 
                                    time, 
                                    digits = 3, 
                                    comp.plot = FALSE,
                                    ...) {
  nvars <- length(object)
  nsims <- object$nsims
  type <- object$type
  modes <- object$modes
  vital <- object$vital
  nts <- max(object$time)
  
  if (missing(time) || (time > max(object$time) | time < 1))
    stop("Specify a timestep between 1 and ", max(object$time))
  
  df.mn <- as.data.frame(object, time=time)
    df.mn <- df.mn[df.mn$time == time, ]
  df.sd <- as.data.frame(object, time=time, out="sd")
    df.sd <- df.sd[df.sd$time == time, ]
  
  if (type == "SI") {
    
    ## Prevalence calcs
    s.prev <- df.mn$s.num/df.mn$num
    i.prev <- df.mn$i.num/df.mn$num
    if (modes == 2) {
      s.prev.m2 <- df.mn$s.num.m2/df.mn$num.m2
      i.prev.m2 <- df.mn$i.num.m2/df.mn$num.m2
    }
    
    ## Group 1 stats
    stats <- c(df.mn$s.num, df.sd$s.num, s.prev,
               df.mn$i.num, df.sd$i.num, i.prev,
               df.mn$num, df.sd$num, 1, 
               df.mn$si.flow, df.sd$si.flow, NA)
    mat <- matrix(stats, byrow=T, nrow=length(stats)/3)
    rownames(mat) <- c("Suscept.",
                       "Infect.",
                       "Total",
                       "S -> I")
    if (vital == TRUE) {
      stats <- c(df.mn$b.flow, df.sd$b.flow, NA,
                 df.mn$ds.flow, df.sd$ds.flow, NA,
                 df.mn$di.flow, df.sd$di.flow, NA)
      mat <- rbind(mat, matrix(stats, byrow=T, nrow=length(stats)/3))
      rownames(mat)[rownames(mat)==""] <- c("Birth ->", "S Death ->", "I Death ->")
    }
    
    ## Group 2 stats
    if (modes == 2) {
      stats <- c(df.mn$s.num.m2, df.sd$s.num.m2, s.prev.m2,
                 df.mn$i.num.m2, df.sd$i.num.m2, i.prev.m2,
                 df.mn$num.m2, df.sd$num.m2, 1, 
                 df.mn$si.flow.m2, df.sd$si.flow.m2, NA)
      mat.m2 <- matrix(stats, byrow=T, nrow=4)
      if (vital == TRUE) {
        stats <- c(df.mn$b.flow.m2, df.sd$b.flow.m2, NA,
                   df.mn$ds.flow.m2, df.sd$ds.flow.m2, NA,
                   df.mn$di.flow.m2, df.sd$di.flow.m2, NA)
        mat.m2 <- rbind(mat.m2, matrix(stats, byrow=T, nrow=length(stats)/3))
      }
      mat <- cbind(mat, mat.m2)
    }
    
  } ## end SI summary
  
  
  if (type == "SIR") {
    
    ## Prevalence calcs
    s.prev <- df.mn$s.num/df.mn$num
    i.prev <- df.mn$i.num/df.mn$num
    r.prev <- df.mn$r.num/df.mn$num
    if (modes == 2) {
      s.prev.m2 <- df.mn$s.num.m2/df.mn$num.m2
      i.prev.m2 <- df.mn$i.num.m2/df.mn$num.m2
      r.prev.m2 <- df.mn$r.num.m2/df.mn$num.m2
    }
    
    ## Group 1 stats
    stats <- c(df.mn$s.num, df.sd$s.num, s.prev,
               df.mn$i.num, df.sd$i.num, i.prev,
               df.mn$r.num, df.sd$r.num, r.prev,
               df.mn$num, df.sd$num, 1,
               df.mn$si.flow, df.sd$si.flow, NA,
               df.mn$ir.flow, df.sd$ir.flow, NA)
    mat <- matrix(stats, byrow=T, nrow=length(stats)/3)
    rownames(mat) <- c("Suscept.",
                       "Infect.",
                       "Recov.",
                       "Total",
                       "S -> I",
                       "I -> R")
    if (vital == TRUE) {
      stats <- c(df.mn$b.flow, df.sd$b.flow, NA,
                 df.mn$ds.flow, df.sd$ds.flow, NA,
                 df.mn$di.flow, df.sd$di.flow, NA,
                 df.mn$dr.flow, df.sd$dr.flow, NA)
      mat <- rbind(mat, matrix(stats, byrow=T, nrow=length(stats)/3))
      rownames(mat)[rownames(mat)==""] <- c("Birth ->", 
                                            "S Death ->", 
                                            "I Death ->",
                                            "R Death ->")
    }
    
    ## Group 2 stats
    if (modes == 2) {
      stats <- c(df.mn$s.num.m2, df.sd$s.num.m2, s.prev.m2,
                 df.mn$i.num.m2, df.sd$i.num.m2, i.prev.m2,
                 df.mn$r.num.m2, df.sd$r.num.m2, r.prev.m2,
                 df.mn$num.m2, df.sd$num.m2, 1,
                 df.mn$si.flow.m2, df.sd$si.flow.m2, NA,
                 df.mn$ir.flow.m2, df.sd$ir.flow.m2, NA)
      mat.m2 <- matrix(stats, byrow=T, nrow=length(stats)/3)
      if (vital == TRUE) {
        stats <- c(df.mn$b.flow.m2, df.sd$b.flow.m2, NA,
                   df.mn$ds.flow.m2, df.sd$ds.flow.m2, NA,
                   df.mn$di.flow.m2, df.sd$di.flow.m2, NA,
                   df.mn$dr.flow.m2, df.sd$dr.flow.m2, NA)
        mat.m2 <- rbind(mat.m2, matrix(stats, byrow=T, nrow=length(stats)/3))
      }
      mat <- cbind(mat, mat.m2)
    }
  }
  
  
  if (type == "SIS") {
    
    ## Prevalence calcs
    s.prev <- df.mn$s.num/df.mn$num
    i.prev <- df.mn$i.num/df.mn$num
    if (modes == 2) {
      s.prev.m2 <- df.mn$s.num.m2/df.mn$num.m2
      i.prev.m2 <- df.mn$i.num.m2/df.mn$num.m2
    }
    
    ## Group 1 stats
    stats <- c(df.mn$s.num, df.sd$s.num, s.prev,
               df.mn$i.num, df.sd$i.num, i.prev,
               df.mn$num, df.sd$num, 1,
               df.mn$si.flow, df.sd$si.flow, NA,
               df.mn$is.flow, df.sd$is.flow, NA)
    mat <- matrix(stats, byrow=T, nrow=length(stats)/3)
    rownames(mat) <- c("Suscept.",
                       "Infect.",
                       "Total",
                       "S -> I",
                       "I -> S")
    if (vital == TRUE) {
      stats <- c(df.mn$b.flow, df.sd$b.flow, NA,
                 df.mn$ds.flow, df.sd$ds.flow, NA,
                 df.mn$di.flow, df.sd$di.flow, NA)
      mat <- rbind(mat, matrix(stats, byrow=T, nrow=length(stats)/3))
      rownames(mat)[rownames(mat)==""] <- c("Birth ->", 
                                            "S Death ->", 
                                            "I Death ->")
    }
    
    ## Group 2 stats
    if (modes == 2) {
      stats <- c(df.mn$s.num.m2, df.sd$s.num.m2, s.prev.m2,
                 df.mn$i.num.m2, df.sd$i.num.m2, i.prev.m2,
                 df.mn$num.m2, df.sd$num.m2, 1,
                 df.mn$si.flow.m2, df.sd$si.flow.m2, NA,
                 df.mn$is.flow.m2, df.sd$is.flow.m2, NA)
      mat.m2 <- matrix(stats, byrow=T, nrow=length(stats)/3)
      if (vital == TRUE) {
        stats <- c(df.mn$b.flow.m2, df.sd$b.flow.m2, NA,
                   df.mn$ds.flow.m2, df.sd$ds.flow.m2, NA,
                   df.mn$di.flow.m2, df.sd$di.flow.m2, NA)
        mat.m2 <- rbind(mat.m2, matrix(stats, byrow=T, nrow=length(stats)/3))
      }
      mat <- cbind(mat, mat.m2)
    }
  }
  
  if (modes == 1) colnames(mat) <- c("mean","sd", "perc")
  if (modes == 2) colnames(mat) <- c("mean:m1", "sd:m1", "perc:m1",
                                      "mean:m2", "sd:m2", "perc:m2")
  mat <- round(mat, digits)
  
  ## Print it
  cat("\nEpiModel Summary")
  cat("\n=======================")
  cat("\nModel class:", class(object))
  
  cat("\n\nSimulation Details")
  cat("\n-----------------------")
  cat("\nModel type:", type)
  cat("\nNo. simulations:", nsims)
  cat("\nNo. time steps:", nts)
  cat("\nNo. NW modes:", modes)
  
  if (modes == 1) statsep <- "---------------------------------"
  if (modes == 2) statsep <- "------------------------------------------------------------------------"
  
  cat("\n\nModel Statistics\n")
  cat(statsep)
  cat("\nTime:", time, "\n")
  cat(statsep, "\n")
  print(mat, print.gap=4)
  cat(statsep, "\n")
  
  if (comp.plot == TRUE) {
    comp.plot(object, 
              time=time,
              digits=digits)
  }
  
}



#' @title Summary for Network Model Fit
#'
#' @description This function prints the summary model fit statistics for an
#'   ERGM or STERGM fit.
#'
#' @param object an \code{EpiModel} object of class \code{epiNet.est}.
#' @param ... additional summary function arguments.
#' 
#' @method summary epiNet.est
#' @keywords extract
#' @export
#' 
#' @details
#' This function is simply a wrapper function for \code{summary.ergm} and 
#'   \code{summary.stergm}.
#' 
summary.epiNet.est <- function(object, 
                               ...) {
  
  summary(object$fit, ...)

  
}


