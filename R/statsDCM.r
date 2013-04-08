#' 
#' @title Extract Model Statistics
#'
#' @description This function extracts and prints model statistics solved 
#' with the \code{mod.} class functions.
#'
#' @param out array containing model results
#' @param time timestep for model statistics
#' @param run run number for model if sensitivity analyses are conducted
#' @param digits number of significant digits to print
#' @param console print to R console
#' 
#' @author Samuel M. Jenness <sjenness@@uw.edu>
#' @keywords model
#' @export
#' 
#' @examples
#' out <- epiDCM(type='SIR', s.num=1000, i.num=1, r.num=0,
#'                 beta=0.2, cont=2:8, nu=1/3,
#'                 b=0.011, ms=0.01, mi=0.03, mr=0.01,
#'                 dt=500, verbose=TRUE)
#' statsDCM(out, time=25, run=1)
#' statsDCM(out, time=25, run=7)
#' statsDCM(out, time=26, run=7)
#' 
statsDCM <- function(out, time, run=1, digits=3, console=TRUE){
    
  nts <- max(out$time)
  nvars <- length(out)-1
  if (class(out$s.num)=='numeric') nruns <- 1
  if (class(out$s.num)=='data.frame') nruns <- dim(out$s.num)[2]
  
  df <- getDCM(out, run)
  df <- df[df$time==time,]

  mytab <- matrix(c(df$s.num, df$s.prev,
                    df$i.num, df$i.prev,
                    df$r.num, df$r.prev,
                    df$i.incid, NA,
                    df$R0, NA,
                    df$Rn, NA), byrow=T, nrow=6)
  rownames(mytab) <- c('S','I','R','Inc.Rate','R0','Rn')
  colnames(mytab) <- c('n','pct')
  
  if (nruns > 1) {
    mytab <- rbind(mytab, c(df[,nvars], NA))
    rownames(mytab)[7] <-  names(df)[nvars]
  }
  
  mytab <- round(mytab, digits)
  
  # print it
  if (console) {
    print(mytab, print.gap=3)
    cat('\n')
  } else {
    return(mytab)
  }
  
}