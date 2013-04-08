#'
#' @title Drop Variables from a Data Frame
#'
#' @description This function drops named variables from a data frame
#'
#' @param df data frame with multiple variables
#' @param vars vector of variables in character format
#' 
#' @author Samuel M. Jenness <sjenness@@uw.edu>
#' @export
#' 
#' @examples
#' w <- x <- y <- z <- 1:10
#' df <- data.frame(w,x,y,z)
#' ( newdf <- dropdf(df, c('x', 'y')) )
#' 
dropdf <- function(df, vars){
  df2 <- df[,!(names(df) %in% vars)]
}

#'
#' @title Expand a Frequency-Weighted Data Frame
#'
#' @description This function expands a frequency-weighted data frame by the specified weights.
#'
#' @param df data frame with one row per combination of variables
#' @param count variable on the data frame with the count variable
#' 
#' @author Samuel M. Jenness <sjenness@@uw.edu>
#' @export
#' 
#' @examples
#' df <- expand.grid(a=0:1, b=0:1, c=0:1)
#' df$count <- rpois(nrow(df), 5)
#' ( df2 <- expanddf(df, df$count) )
#' 
expanddf <- function(df, count){
  df2 <- as.data.frame(lapply(df, function(x) rep(x, count)))
}

#' 
#' @title Summarizing Regression Model Fits
#'
#' @description Provides summary statistics from a regression model fit, 
#' including exponentiated coefficients, robust standard errors, and
#' confidence intervals.
#'
#' @param fit fitted regression model
#' @param exp exponentiate the summary coefficients
#' @param robust robust standard errors
#' @param digits digits to provide in output
#' 
#' @author Samuel M. Jenness <sjenness@@uw.edu>
#' @export
#' 
#' @examples
#' y <- rbinom(100, 1, 0.5)
#' x <- ifelse(y == 1, rnorm(length(y[y==1]), 100, 20), 
#'             rnorm(length(y[y==0]), 50, 20))
#' mod <- glm(y ~ x, family=binomial)
#' summary(mod)
#' mysummary(mod, exp=TRUE, robust=TRUE)
#' 
mysummary <- function(fit, exp=FALSE, robust=FALSE, digits=3) {
  require(sandwich)
  c <- coef(fit)
  if (robust == TRUE) {
    s <- sqrt(diag(vcovHC(fit, 'HC1')))
  }
  if (robust == FALSE) {
    s <- sqrt(diag(vcov(fit)))
  }
  w <- c/s
  df <- length(fitted(fit))-length(c)
  p <- pt(abs(w), df, lower.tail=FALSE)*2
  i <- cbind(c,c) + t(qt(c(0.025, 0.975), df) %o% s)
  if (exp==F) {
    o <- cbind('Coef'=c, 'SE'=s, 't'=w, 'p'=p, 
               '95% CI'=i[,1], ' '=i[,2])
  }
  if (exp==T) {
    o <- cbind('Coef'=exp(c), 'SE'=s, 't'=w, 'p'=p, 
               '95% CI'=exp(i[,1]), ' '=exp(i[,2]))
  }
  print(round(o, digits))
  invisible(o)
}

#' 
#' @title Sample Rows of a Data Frame
#'
#' @description This function samples rows from a data frame.
#'
#' @param df data frame to be sampled
#' @param size number of rows to sample
#' @param replace sample with replacement
#' @param sortby sort the data frame by a variable in the data frame
#' 
#' @author Samuel M. Jenness <sjenness@@uw.edu>
#' @export
#' 
#' @examples
#' df <- expand.grid(a=0:1, b=0:1, c=0:1)
#' df$count <- rpois(nrow(df), 5)
#' df2 <- expanddf(df, df$count)
#' ( df3 <- sampledf(df2, size=5, replace=FALSE, sortby='c') )
#' 
sampledf <- function(df, size, replace=FALSE, sortby=''){
  require(reshape)
  df2 <- df[sample(nrow(df), size=size, replace=replace),]
  if (sortby!='') {
    df2 <- sort_df(df2, sortby)
    #df2 <- df2[do.call("order", df2[, sortby, drop = FALSE]), , drop = FALSE]
  }
  return(df2)
}