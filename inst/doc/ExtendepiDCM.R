
## ----, include=FALSE-----------------------------------------------------
library(knitr)
opts_chunk$set(cache=TRUE, comment=NA, tidy=FALSE, fig.width=9.5, fig.height=6)


## ------------------------------------------------------------------------
library(deSolve)


## ------------------------------------------------------------------------
SI <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {    
    
    ## Dynamic Calculations
    N <- S+I
    lambda <- rho*c*I/N
    
    ## Differential Equations
    dS <- -lambda*S
    dI <- lambda*S 
    
    ## Output
    list(c(dS, dI), 
         N = N, 
         lambda = lambda)
  })
}


## ------------------------------------------------------------------------
params <- c(c = 4,
            rho = 0.2)
t0 <- c(S = 999, 
        I = 1)
dt <- seq(from = 0, 
          to = 25, 
          by = 0.25)


## ------------------------------------------------------------------------
df <- data.frame(
  ode(y = t0, times = dt, func = SI, parms = params, method = 'rk4')
)


## ------------------------------------------------------------------------
head(df, 20)
df$I[1:100]


## ------------------------------------------------------------------------
par(mar=c(3.2,3.2,2,1), mgp=c(2,1,0))
plot(x=df$time, y=df$S, type='l', col='steelblue', lwd=3, xlab='Time', ylab='Prevalence')
lines(x=df$time, y=df$I, col='firebrick', lwd=3)
legend(x=22, y=500, legend=c('Sus', 'Inf'), lty=1, lwd=3, col=c('steelblue', 'firebrick'), cex=0.9, bty='n')


## ------------------------------------------------------------------------
Qmod <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {    
    
    ## Dynamic Calculations ##
    
    # Popsize
    N.high <- S.high + I.high
    N.low <- S.low + I.low
    N <- N.high + N.low
    prev <- (I.high+I.low)/N
    
    # Contact rates
    c.high <- (c.mean*N - c.low*N.low)/N.high
    
    # mixing matrix calculations based on Q
    g.hh <- ((c.high*N.high) + (Q*c.low*N.low)) / ((c.high*N.high) + (c.low*N.low))
    g.lh <- 1 - g.hh 
    g.hl <- (1 - g.hh) * ((c.high*N.high) / (c.low*N.low))
    g.ll <- 1 - g.hl
    
    # prob that p is infected
    p.high <- (g.hh*I.high/N.high)+(g.lh*I.low/N.low)
    p.low <- (g.ll*I.low/N.low)+(g.hl*I.high/N.high)
    
    # lambda
    lambda.high <- rho * c.high * p.high
    lambda.low <- rho * c.low * p.low
    
    
    ## Differential Equations ##
    dS.high <- -lambda.high*S.high + nu*I.high
    dI.high <- lambda.high*S.high - nu*I.high
    
    dS.low <- -lambda.low*S.low + nu*I.low
    dI.low <- lambda.low*S.low - nu*I.low
    
    
    ## Output ##
    list(c(dS.high, dI.high,
           dS.low, dI.low), 
           N = N,
           prev = prev)
  })
}


## ------------------------------------------------------------------------
params <- c(
  c.mean = 2,
  c.low = 1.4,
  rho = 0.75,
  nu = 6
)
p1 <- c(params, Q = -0.45)
p2 <- c(params, Q = -0.33)
p3 <- c(params, Q = 0)
p4 <- c(params, Q = 0.5)
p5 <- c(params, Q = 1)


## ------------------------------------------------------------------------
N.tot <- 20000000
prop.high <- 0.02
prop.low <- 0.98

t0 <- c(
  S.high = N.tot*prop.high - 1,
  I.high = 1,
  S.low = N.tot*prop.low - 1,
  I.low = 1
)
dt <- seq(1, 25, 0.02)


## ------------------------------------------------------------------------
df1 <- data.frame(ode(y=t0, times=dt, func=Qmod, parms=p1, method='rk4'))
df2 <- data.frame(ode(y=t0, times=dt, func=Qmod, parms=p2, method='rk4'))
df3 <- data.frame(ode(y=t0, times=dt, func=Qmod, parms=p3, method='rk4'))
df4 <- data.frame(ode(y=t0, times=dt, func=Qmod, parms=p4, method='rk4'))
df5 <- data.frame(ode(y=t0, times=dt, func=Qmod, parms=p5, method='rk4'))


## ----, eval=FALSE--------------------------------------------------------
## Qsens <- function(params, Q, t0, Qmod) {
##   temp.params <- c(params, Q = Q)
##   df <- data.frame(ode(y=t0, times=dt, func=Qmod, parms=temp.params, method='rk4'))
##   return(df)
## }
## df1 <- Qsens(params, Q = -0.45, t0, Qmod)
## df2 <- Qsens(params, Q = -0.33, t0, Qmod)
## df3 <- Qsens(params, Q = 0, t0, Qmod)
## df4 <- Qsens(params, Q = 0.5, t0, Qmod)
## df5 <- Qsens(params, Q = 1, t0, Qmod)


## ------------------------------------------------------------------------
library(RColorBrewer)
pal <- brewer.pal(5, 'Set1')
par(mar=c(3.2,3.2,2,1), mgp=c(2,1,0))


## ------------------------------------------------------------------------
plot(df1$time, df1$prev, type='l', ylim=c(0,0.04), col=pal[1], lwd=3, xlab='Time', ylab='Prevalence')
lines(df2$time, df2$prev, col=pal[2], lwd=3)
lines(df3$time, df3$prev, col=pal[3], lwd=3)
lines(df4$time, df4$prev, col=pal[4], lwd=3)
lines(df5$time, df5$prev, col=pal[5], lwd=3)
legend('topright', legend=paste('Q =', c(-0.45, -0.33, 0, 0.5, 1)), lty=1, lwd=3, col=pal, cex=0.9, bg='white')


## ------------------------------------------------------------------------
## Set a basic SI model
SI <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {    
    
    ## Dynamic Calculations
    # Population size
    N <- S+I
    
    # Intervention start time: if time > start, 
    #   then multiply lambda by relative hazard
    if (t < start) {
      lambda <- rho*c*I/N
    } else {
      lambda <- (rho*c*I/N) * rel.haz
    }
    
    ## Differential Equations
    dS <- -lambda*S
    dI <- lambda*S 
    
    ## Output
    list(c(dS, dI), 
         N = N, 
         lambda = lambda)
  })
}

## Try it with an SIS model
SIS <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {    
    
    ## Dynamic Calculations
    
    # Population size
    N <- S+I
    
    # Intervention start time: if time > start, 
    #   then multiply lambda by relative hazard
    if (t < start) {
      lambda <- rho*c*I/N
    } else {
      lambda <- (rho*c*I/N) * rel.haz
    }
    
    ## Differential Equations
    dS <- -lambda*S + nu*I
    dI <- lambda*S - nu*I
    
    ## Output
    list(c(dS, dI), 
         N = N, 
         lambda = lambda)
  })
}

## This function wraps calcuations over a vector of start.times
time.sens <- function(mod, params, t0, dt, start.times, rel.haz) {
  
  # Create an empty matrix for storing results
  out <- matrix(rep(NA, max(dt)*length(start.times)), ncol=length(start.times))
  
  # The base parameters are those fixed for all start times (which varies)
  base.params <- c(params, rel.haz=rel.haz)
  
  # Loop over all start times
  for (i in 1:length(start.times)) {
    temp.params <- c(base.params, start = start.times[i])
    df <- data.frame(ode(y = t0, times = dt, func = mod, parms = temp.params, method = 'rk4'))
    out[,i] <- df$I
  }
  
  # Set up a plot for each of the start times
  pal <- rainbow(length(start.times))
  plot(out[,1], type='n', col=pal[1], lwd=3, ylim=c(0,max(out)),
       xlab='Time', ylab='Number Infected')
  for (i in 1:ncol(out)) {
    lines(x=dt, y=out[,i], col=pal[i], lwd=3)
  }
  legend('bottomright', paste('Start=', start.times, sep=''), 
         lty=1, col=pal, lwd=4, cex=0.8)
  
  # Return the model output to a matrix
  return(out)
}

## SI model tests
params <- c(c = 0.5,
            rho = 0.1)
t0 <- c(S = 999, 
        I = 1)
dt <- 1:1000
start.times <- seq(50, 200, 50)

# Baseline model (intervention has no effect)
out <- time.sens(mod=SI, params, t0, dt, start.times, rel.haz=1)

# Counterfactual models (intervention reduces risk by 50%)
out <- time.sens(mod=SI, params, t0, dt, start.times, rel.haz=0.5)

## SIS model tests
params <- c(c = 0.5,
            rho = 0.2,
            nu = 0.05)
t0 <- c(S = 999, 
        I = 1)
dt <- 1:1000
start.times <- seq(50, 200, 50)

# Baseline model (intervention has no effect)
out <- time.sens(mod=SIS, params, t0, dt, start.times, rel.haz=1)

# Counterfacutal model (intervention reduces risk by 25%)
out <- time.sens(mod=SIS, params, t0, dt, start.times, rel.haz=0.75)



