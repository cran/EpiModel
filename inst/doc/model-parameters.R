## ---- echo = FALSE, include = FALSE-------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo = FALSE, include = FALSE-------------------------------------
library(EpiModel)

## ----model--------------------------------------------------------------------
nw <- network_initialize(n = 50)

est <- netest(
  nw, formation = ~edges,
  target.stats = 25,
  coef.diss = dissolution_coefs(~offset(edges), 10, 0),
  verbose = FALSE
)

param <- param.net(
  inf.prob = 0.3,
  act.rate = 0.5,
  dummy.param = 4,
  dummy.strat.param = c(0, 1)
)

init <- init.net(i.num = 10)
control <- control.net(type = "SI", nsims = 1, nsteps = 5, verbose = FALSE)
mod <- netsim(est, param, init, control)
mod

## ----generators---------------------------------------------------------------
my.randoms <- list(
  act.rate = param_random(c(0.25, 0.5, 0.75)),
  dummy.param = function() rbeta(1, 1, 2),
  dummy.strat.param = function() c(
    rnorm(1, 0.05, 0.01),
    rnorm(1, 0.15, 0.03)
  )
)

param <- param.net(
  inf.prob = 0.3,
  random.params = my.randoms
)

param

## ----generators_run-----------------------------------------------------------
control <- control.net(type = "SI", nsims = 3, nsteps = 5, verbose = FALSE)
mod <- netsim(est, param, init, control)

mod

## ----generators_inspect-------------------------------------------------------
str(mod$param$random.params.values)

## ----set_df-------------------------------------------------------------------
n <- 5

related.param <- data.frame(
  dummy.param = rbeta(n, 1, 2)
)

related.param$dummy.strat.param_1 <- related.param$dummy.param + rnorm(n)
related.param$dummy.strat.param_2 <- related.param$dummy.param * 2 + rnorm(n)

related.param

## ----set_param----------------------------------------------------------------
my.randoms <- list(
  act.rate = param_random(c(0.25, 0.5, 0.75)),
  param_random_set = related.param

)

param <- param.net(
  inf.prob = 0.3,
  random.params = my.randoms
)

param

## ----set_run------------------------------------------------------------------
control <- control.net(type = "SI", nsims = 3, nsteps = 5, verbose = FALSE)
mod <- netsim(est, param, init, control)

mod

## ----updater-example, eval = FALSE--------------------------------------------
#  list(
#    at = 10,
#    param = list(
#      inf.prob = 0.3,
#      act.rate = 0.5
#    )
#  )

## ----updater-list-------------------------------------------------------------
# Create a `list.of.updaters`
list.of.updaters <- list(
  # this is one updater
  list(
    at = 100,
    param = list(
      inf.prob = 0.3,
      act.rate = 0.3
    )
  ),
  # this is another updater
  list(
    at = 125,
    param = list(
      inf.prob = 0.01
    )
  )
)

 # The `list.of.updaters` goes into `param.net` under `param.updater.list`
 param <- param.net(
   inf.prob = 0.1,
   act.rate = 0.1,
   param.updater.list = list.of.updaters
 )

## ----updater-module, fig.align = "center", fig.height = 4, fig.width = 6------
 control <- control.net(
   type = NULL, # must be NULL as we use a custom module
   nsims = 1,
   nsteps = 200,
   verbose = FALSE,
   updater.FUN = updater.net,
   infection.FUN = infection.net
 )

nw <- network_initialize(n = 50)
nw <- set_vertex_attribute(nw, "race", rbinom(50, 1, 0.5))
est <- netest(
  nw,
  formation = ~edges,
  target.stats = 25,
  coef.diss = dissolution_coefs(~offset(edges), 10, 0),
  verbose = FALSE
)

init <- init.net(i.num = 10)
mod <- netsim(est, param, init, control)

plot(mod)

## ----updater-verbose, echo = FALSE--------------------------------------------
list(
  at = 10,
  param = list(
    inf.prob = 0.3,
    act.rate = 0.5
  ),
  verbose = TRUE
)

## ----updater-function, echo = FALSE-------------------------------------------
list(
  at = 10,
  param = list(
    inf.prob = function(x) plogis(qlogis(x) + log(2)),
    act.rate = 0.5
  )
)

