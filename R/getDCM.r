#' 
#' @title Extract Model Run
#'
#' @description This function extracts a specific model run for epiDCM 
#'
#' @param out array containing model results
#' @param run run number for model if sensitivity analyses are conducted
#' 
#' @author Samuel M. Jenness <sjenness@@uw.edu>
#' @keywords model
#' @export
#' 
#' @examples
#' out <- epiDCM(type='SIR', s.num=1000, i.num=1, r.num=0,
#'                beta=0.2, cont=1:4, nu=1/3,
#'                b=0.011, ms=0.01, mi=0.03, mr=0.01,
#'                dt=500, verbose=TRUE)
#' mod.r5 <- getDCM(out, run=3)
#' head(mod.r5)
#' 
getDCM <- function(out, run=1){
  
  if (class(out$s.num) == 'numeric') {
    if (run > 1) stop('\n Max run for model is 1')
    new.out <- as.data.frame(out)
  }
  if (class(out$s.num) == 'data.frame') {
    if (run > ncol(out$s.num)) stop(paste('\n Max run for model is', ncol(out$s.num)))
    new.out <- data.frame(time=out$time)
    for (i in 3:length(out)) new.out[,(i-1)] <- out[[i]][run]
    names(new.out) <- names(out)[2:length(out)]
  } 
  
  return(new.out)
}