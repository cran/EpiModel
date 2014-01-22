
##
## Testing Series for epiNet
##


run.tests.epiNet <- function(plots=FALSE, seed=FALSE) {

require(EpiModel)
require(testthat)

if (plots == TRUE) pdf('EpiModel/epiNet.pdf', h=5, w=10)
par(mar=c(5,5,3,1), mgp=c(3,2,1))

cat('\n===============================')
cat('\nTesting Series for epiNet')
cat('\n===============================')


# Network Estimation ------------------------------------------------------

test <- 'Estimation: edges only, edapprox'
EpiModel:::mcat(test)

nw <- network.initialize(n = 100, directed = FALSE)
est <- epiNet.est(
  nw,
  formation = ~ edges,
  dissolution = ~offset(edges),
  target.stats = 25,
  coef.diss = dissolution.coefs(~offset(edges), 10, 0),
  save.stats = FALSE,
  verbose = FALSE)

expect_is(est, 'epiNet.est')

capture.output(print(est), 
               summary(est),
               file='NULL')
invisible(file.remove('NULL'))
cat('passed')

############

test <- 'Estimation: edges + degree0, edapprox'
EpiModel:::mcat(test)

nw <- network.initialize(n = 100, directed = FALSE)
est2 <- epiNet.est(
  nw,
  formation = ~ edges + degree(1),
  dissolution = ~offset(edges),
  target.stats = c(25, 25),
  coef.diss = dissolution.coefs(~offset(edges), 10, 0),
  save.stats = FALSE, 
  verbose = FALSE)

expect_is(est2, 'epiNet.est')

capture.output(print(est2),
               summary(est2),
               file='NULL')
invisible(file.remove('NULL'))
cat('passed')

############

test <- 'Estimation: edges + nodematch, edapprox'
EpiModel:::mcat(test)

nw <- network.initialize(n = 100, directed = FALSE)
nw %v% 'sex' <- rep(c('F', 'M'), each=50)
est3 <- epiNet.est(
  nw,
  formation = ~ edges + nodematch('sex'),
  dissolution = ~offset(edges),
  target.stats = c(25, 10),
  coef.diss = dissolution.coefs(~offset(edges), 10, 0),
  save.stats = FALSE, 
  verbose = FALSE)

expect_is(est3, 'epiNet.est')

capture.output(print(est3), 
               summary(est3), 
               file='NULL')
invisible(file.remove('NULL'))
cat('passed')

############

test <- 'Estimation: edges only, edapprox=F'
EpiModel:::mcat(test)

nw <- network.initialize(n = 100, directed = FALSE)
est4 <- epiNet.est(
  nw,
  formation = ~ edges,
  dissolution = ~offset(edges),
  target.stats = 25,
  coef.diss = dissolution.coefs(~offset(edges), 10, 0),
  save.stats = FALSE,
  edapprox = FALSE,
  verbose = FALSE)

expect_is(est4, 'epiNet.est')

capture.output(print(est4), 
               summary(est4),
               file='NULL')
invisible(file.remove('NULL'))
cat('passed')

############

test <- 'Estimation: bipartite, edapprox=T'
EpiModel:::mcat(test)

nw <- network.initialize(n = 100, bipartite=50, directed = FALSE)
est5 <- epiNet.est(
  nw,
  formation = ~ edges,
  dissolution = ~offset(edges),
  target.stats = 25,
  coef.diss = dissolution.coefs(~offset(edges), 10, 0),
  save.stats = FALSE,
  edapprox = TRUE,
  verbose = FALSE)

expect_is(est5, 'epiNet.est')

capture.output(print(est5), 
               summary(est5), 
               file='NULL')
invisible(file.remove('NULL'))
cat('passed')

############

test <- 'Estimation: bipartite, full STERGM'
EpiModel:::mcat(test)

nw <- network.initialize(n = 100, bipartite=50, directed = FALSE)
est6 <- epiNet.est(
  nw,
  formation = ~ edges,
  dissolution = ~offset(edges),
  target.stats = 25,
  coef.diss = dissolution.coefs(~offset(edges), 10, 0),
  save.stats = FALSE,
  edapprox = FALSE,
  verbose = FALSE)

expect_is(est6, 'epiNet.est')

capture.output(print(est6), 
               summary(est6), 
               file='NULL')
invisible(file.remove('NULL'))
cat('passed')



# Network simulation ------------------------------------------------------

test <- 'NW simulation: 25 steps, 1 sim'
EpiModel:::mcat(test)
  
nwsims <- epiNet.simNet(
  est, 
  nsteps = 25, 
  verbose = FALSE)

expect_is(nwsims, 'epiNet.simNet')

capture.output(print(nwsims), file='NULL')
invisible(file.remove('NULL'))
cat('passed')

############

test <- 'NW simulation: 25 steps, 5 sim'
EpiModel:::mcat(test)

nwsims2 <- epiNet.simNet(
  est, 
  nsteps = 25,
  nsim = 5,
  verbose = FALSE)

expect_is(nwsims2, 'epiNet.simNet')

capture.output(print(nwsims2), file='NULL')
invisible(file.remove('NULL'))
cat('passed')

############

test <- 'NW simulation: 25 steps, 5 sim, start from edapprox=F'
EpiModel:::mcat(test)

nwsims3 <- epiNet.simNet(
  est4, 
  nsteps = 25,
  nsim = 5,
  verbose = FALSE)

expect_is(nwsims3, 'epiNet.simNet')

capture.output(print(nwsims3), file='NULL')
invisible(file.remove('NULL'))
cat('passed')

############

test <- 'NW simulation: 25 steps, 5 sim, bipartite'
EpiModel:::mcat(test)

nwsims4 <- epiNet.simNet(
  est5, 
  nsteps = 25,
  nsim = 5,
  verbose = FALSE)

expect_is(nwsims4, 'epiNet.simNet')

capture.output(print(nwsims4), file='NULL')
invisible(file.remove('NULL'))
cat('passed')


# SI Models ---------------------------------------------------------------


test <- 'SI model, 1 mode, no vital, 1 sim'
EpiModel:::mcat(test)
  
x <- epiNet.simTrans(
  nwsims, 
  type = 'SI',
  i.num = 1, 
  trans.rate = 0.5,
  sims.per.nw = 1,
  verbose = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_true(max(x$i.num) >= 1)
expect_true(max(x$i.num) <= 100)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)

############

test <- 'SI model, 1 mode, no vital, 1 sim, trans.rate=0'
EpiModel:::mcat(test)
  
x <- epiNet.simTrans(
  nwsims, 
  type = 'SI',
  i.num = 1, 
  trans.rate = 0,
  sims.per.nw = 1,
  i.rand = FALSE,
  verbose = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_true(max(x$i.num) == 1)
expect_true(max(x$si.flow) == 0)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)

############

test <- 'SI model, 1 mode, no vital, 5 sim'
EpiModel:::mcat(test)
  
x <- epiNet.simTrans(
  nwsims, 
  type = 'SI',
  i.num = 1, 
  trans.rate = 0.5,
  sims.per.nw = 5,
  verbose = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_true(max(x$i.num) >= 1)
expect_true(max(x$i.num) <= 100)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)

############

test <- 'SI model, 1 mode, no vital, 5 sim, use TEAs'
EpiModel:::mcat(test)
  
x <- epiNet.simTrans(
  nwsims, 
  type = 'SI',
  i.num = 1, 
  trans.rate = 1,
  sims.per.nw = 5,
  verbose = FALSE,
  i.rand = FALSE,
  tea = TRUE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_true(max(x$i.num) > 1)
expect_true(max(x$i.num) <= 100)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)

############

test <- 'SI model, 1 mode, no vital, 5 sim, set i.ids'
EpiModel:::mcat(test)
  
x <- epiNet.simTrans(
  nwsims, 
  type = 'SI',
  i.ids = 1:10, 
  trans.rate = 0.5,
  sims.per.nw = 5,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')

# this tests that starting infected is same across sims
expect_true(all(x$stat.mat[[1]][1, ] == x$stat.mat[[2]][1, ]))
expect_true(min(x$i.num) == 10)
expect_true(max(x$i.num) <= 100)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)

############

test <- 'SI model, 1 mode, no vital, 2 sim, use full STERGM fit model'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  nwsims3, 
  type = 'SI',
  i.num = 10, 
  trans.rate = 0.5,
  sims.per.nw = 2,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_true(max(x$i.num) >= 1)
expect_true(max(x$i.num) <= 100)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)


############

test <- 'SI model, bipartite, no vital, 2 sim'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  nwsims4, 
  type = 'SI',
  i.num = 10, 
  i.rand = FALSE,
  trans.rate = 0.5,
  trans.rate.m2 = 0.1,
  sims.per.nw = 2,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_equal(x$modes, 2)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)


############

test <- 'SI model, vital, 1 sim'
EpiModel:::mcat(test)

nw <- network.initialize(n = 100, directed = FALSE)
est.vit <- epiNet.est(
  nw,
  formation = ~ edges,
  dissolution = ~offset(edges),
  target.stats = 25,
  coef.diss = dissolution.coefs(~offset(edges), 10, 0.01),
  save.stats = FALSE,
  verbose = FALSE)
x <- epiNet.simTrans(
  est.vit, 
  type = 'SI',
  vital = TRUE,
  i.num = 10, 
  trans.rate = 0.5,
  act.rate = 2,
  b.rate = 0.01,
  ds.rate = 0.01,
  di.rate = 0.01,
  nsteps = 25,
  sims.per.nw = 1,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_equal(x$vital, TRUE)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)


############

test <- 'SI model, vital, 3 sim'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  est.vit, 
  type = 'SI',
  vital = TRUE,
  i.num = 10, 
  trans.rate = 0.5,
  act.rate = 2,
  b.rate = 0.01,
  ds.rate = 0.01,
  di.rate = 0.01,
  nsteps = 25,
  sims.per.nw = 3,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_equal(x$vital, TRUE)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)

############

test <- 'SI model, vital, 1 sim, bipartite'
EpiModel:::mcat(test)

nw <- network.initialize(n = 100, bipartite=50, directed = FALSE)
est5.vit <- epiNet.est(
  nw,
  formation = ~ edges,
  dissolution = ~offset(edges),
  target.stats = 25,
  coef.diss = dissolution.coefs(~offset(edges), 10, 0.01),
  save.stats = FALSE,
  edapprox = TRUE,
  verbose = FALSE)
x <- epiNet.simTrans(
  est5.vit, 
  type = 'SI',
  vital = TRUE,
  i.num = 10, 
  trans.rate = 0.5,
  trans.rate.m2 = 0.1,
  act.rate = 2,
  b.rate = 0.01,
  ds.rate = 0.01,
  di.rate = 0.01,
  nsteps = 25,
  sims.per.nw = 1,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_equal(x$vital, TRUE)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)


############

test <- 'SI model, vital, 3 sim, bipartite'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  est5.vit, 
  type = 'SI',
  vital = TRUE,
  i.num = 10, 
  trans.rate = 0.5,
  trans.rate.m2 = 0.1,
  act.rate = 2,
  b.rate = 0.01,
  ds.rate = 0.01,
  di.rate = 0.01,
  nsteps = 25,
  sims.per.nw = 3,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_equal(x$vital, TRUE)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)


# SIR models --------------------------------------------------------------

test <- 'SIR model, 1 mode, no vital, 1 sim'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  nwsims, 
  type = 'SIR',
  i.num = 1, 
  trans.rate = 0.5,
  rec.rate = 0.1,
  sims.per.nw = 1,
  verbose = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_true(max(x$i.num) >= 1)
expect_true(max(x$i.num) <= 100)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)

############

test <- 'SIR model, 1 mode, no vital, 1 sim, trans.rate=0'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  nwsims, 
  type = 'SIR',
  i.num = 1, 
  trans.rate = 0,
  rec.rate = 0.1,
  sims.per.nw = 1,
  i.rand = FALSE,
  verbose = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_true(max(x$i.num) == 1)
expect_true(max(x$si.flow) == 0)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)

############

test <- 'SIR model, 1 mode, no vital, 5 sim'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  nwsims, 
  type = 'SIR',
  i.num = 1, 
  trans.rate = 0.5,
  rec.rate = 0.1,
  sims.per.nw = 5,
  verbose = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_true(max(x$i.num) >= 1)
expect_true(max(x$i.num) <= 100)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)

############

test <- 'SIR model, 1 mode, no vital, 5 sim, use TEAs'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  nwsims, 
  type = 'SIR',
  i.num = 1, 
  trans.rate = 1,
  rec.rate = 0.1,
  sims.per.nw = 5,
  verbose = FALSE,
  i.rand = FALSE,
  tea = TRUE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_true(max(x$i.num) >= 1)
expect_true(max(x$i.num) <= 100)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)

############

test <- 'SIR model, 1 mode, no vital, 5 sim, set i.ids'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  nwsims, 
  type = 'SIR',
  i.ids = 1:10, 
  trans.rate = 0.5,
  rec.rate = 0.1,
  sims.per.nw = 5,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')

# this tests that starting infected is same across sims
expect_true(all(x$stat.mat[[1]][1, ] == x$stat.mat[[2]][1, ]))
expect_true(max(x$i.num) <= 100)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)

############

test <- 'SIR model, 1 mode, no vital, 2 sim, use full STERGM fit model'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  nwsims3, 
  type = 'SIR',
  i.num = 10, 
  trans.rate = 0.5,
  rec.rate = 0.1,
  sims.per.nw = 2,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_true(max(x$i.num) >= 1)
expect_true(max(x$i.num) <= 100)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)


############

test <- 'SIR model, bipartite, no vital, 2 sim'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  nwsims4, 
  type = 'SIR',
  i.num = 10, 
  i.rand = FALSE,
  trans.rate = 0.5,
  trans.rate.m2 = 0.1,
  rec.rate = 0.1,
  sims.per.nw = 2,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_equal(x$modes, 2)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)


############

test <- 'SIR model, vital, 1 sim'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  est.vit, 
  type = 'SIR',
  vital = TRUE,
  i.num = 10, 
  trans.rate = 0.5,
  rec.rate = 0.1,
  act.rate = 2,
  b.rate = 0.01,
  ds.rate = 0.01,
  di.rate = 0.01,
  dr.rate = 0.01,
  nsteps = 25,
  sims.per.nw = 1,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_equal(x$vital, TRUE)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)


############

test <- 'SIR model, vital, 3 sim'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  est.vit, 
  type = 'SIR',
  vital = TRUE,
  i.num = 10, 
  trans.rate = 0.5,
  rec.rate = 0.1,
  act.rate = 2,
  b.rate = 0.01,
  ds.rate = 0.01,
  di.rate = 0.01,
  dr.rate = 0.01,
  nsteps = 25,
  sims.per.nw = 3,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_equal(x$vital, TRUE)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)

############

test <- 'SIR model, vital, 1 sim, bipartite'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  est5.vit, 
  type = 'SIR',
  vital = TRUE,
  i.num = 10, 
  trans.rate = 0.5,
  trans.rate.m2 = 0.1,
  rec.rate = 0.1,
  act.rate = 2,
  b.rate = 0.01,
  ds.rate = 0.01,
  di.rate = 0.01,
  dr.rate = 0.01,
  nsteps = 25,
  sims.per.nw = 1,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_equal(x$vital, TRUE)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)


############

test <- 'SIR model, vital, 3 sim, bipartite'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  est5.vit, 
  type = 'SIR',
  vital = TRUE,
  i.num = 10, 
  trans.rate = 0.5,
  trans.rate.m2 = 0.1,
  rec.rate = 0.1,
  act.rate = 2,
  b.rate = 0.01,
  ds.rate = 0.01,
  di.rate = 0.01,
  dr.rate = 0.01,
  nsteps = 25,
  sims.per.nw = 3,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_equal(x$vital, TRUE)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)


# SIS models --------------------------------------------------------------

test <- 'SIS model, 1 mode, no vital, 1 sim'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  nwsims, 
  type = 'SIS',
  i.num = 1, 
  trans.rate = 0.9,
  rec.rate = 0.01,
  sims.per.nw = 1,
  verbose = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_true(max(x$i.num) >= 1)
expect_true(max(x$i.num) <= 100)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)

############

test <- 'SIS model, 1 mode, no vital, 1 sim, trans.rate=0'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  nwsims, 
  type = 'SIS',
  i.num = 1, 
  trans.rate = 0,
  rec.rate = 0.01,
  sims.per.nw = 1,
  i.rand = FALSE,
  verbose = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_true(max(x$i.num) == 1)
expect_true(max(x$si.flow) == 0)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)

############

test <- 'SIS model, 1 mode, no vital, 5 sim'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  nwsims, 
  type = 'SIS',
  i.num = 1, 
  trans.rate = 0.5,
  rec.rate = 0.01,
  sims.per.nw = 5,
  verbose = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_true(max(x$i.num) >= 1)
expect_true(max(x$i.num) <= 100)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)

############

test <- 'SIS model, 1 mode, no vital, 5 sim, use TEAs'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  nwsims, 
  type = 'SIS',
  i.num = 1, 
  trans.rate = 1,
  rec.rate = 0.01,
  sims.per.nw = 5,
  verbose = FALSE,
  i.rand = FALSE,
  tea = TRUE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_true(max(x$i.num) > 1)
expect_true(max(x$i.num) <= 100)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)

############

test <- 'SIS model, 1 mode, no vital, 5 sim, set i.ids'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  nwsims, 
  type = 'SIS',
  i.ids = 1:10, 
  trans.rate = 0.5,
  rec.rate = 0.01,
  sims.per.nw = 5,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')

# this tests that starting infected is same across sims
expect_true(all(x$stat.mat[[1]][1, ] == x$stat.mat[[2]][1, ]))
expect_true(max(x$i.num) <= 100)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)

############

test <- 'SIS model, 1 mode, no vital, 2 sim, use full STERGM fit model'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  nwsims3, 
  type = 'SIS',
  i.num = 10, 
  trans.rate = 0.5,
  rec.rate = 0.01,
  sims.per.nw = 2,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_true(max(x$i.num) >= 1)
expect_true(max(x$i.num) <= 100)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)


############

test <- 'SIS model, bipartite, no vital, 2 sim'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  nwsims4, 
  type = 'SIS',
  i.num = 10, 
  i.rand = FALSE,
  trans.rate = 0.5,
  trans.rate.m2 = 0.25,
  rec.rate = 0.01,
  sims.per.nw = 2,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_equal(x$modes, 2)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)


############

test <- 'SIS model, vital, 1 sim'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  est.vit, 
  type = 'SIS',
  vital = TRUE,
  i.num = 10, 
  trans.rate = 0.5,
  rec.rate = 0.01,
  act.rate = 2,
  b.rate = 0.01,
  ds.rate = 0.01,
  di.rate = 0.01,
  nsteps = 25,
  sims.per.nw = 1,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_equal(x$vital, TRUE)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)


############

test <- 'SIS model, vital, 3 sim'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  est.vit, 
  type = 'SIS',
  vital = TRUE,
  i.num = 10, 
  trans.rate = 0.5,
  rec.rate = 0.01,
  act.rate = 2,
  b.rate = 0.01,
  ds.rate = 0.01,
  di.rate = 0.01,
  nsteps = 25,
  sims.per.nw = 3,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_equal(x$vital, TRUE)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)

############

test <- 'SIS model, vital, 1 sim, bipartite'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  est5.vit, 
  type = 'SIS',
  vital = TRUE,
  i.num = 10, 
  trans.rate = 0.5,
  trans.rate.m2 = 0.25,
  rec.rate = 0.01,
  act.rate = 2,
  b.rate = 0.01,
  ds.rate = 0.01,
  di.rate = 0.01,
  nsteps = 25,
  sims.per.nw = 1,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_equal(x$vital, TRUE)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)


############

test <- 'SIS model, vital, 3 sim, bipartite'
EpiModel:::mcat(test)

x <- epiNet.simTrans(
  est5.vit, 
  type = 'SIS',
  vital = TRUE,
  i.num = 10, 
  trans.rate = 0.5,
  trans.rate.m2 = 0.25,
  rec.rate = 0.01,
  act.rate = 2,
  b.rate = 0.01,
  ds.rate = 0.01,
  di.rate = 0.01,
  nsteps = 25,
  sims.per.nw = 3,
  verbose = FALSE,
  tea = FALSE) 

expect_is(x, 'epiNet.simTrans')
expect_is(as.data.frame(x), 'data.frame')
expect_equal(x$vital, TRUE)

capture.output(print(x),
               summary(x, time=10),
               file='NULL')
invisible(file.remove('NULL'))
if (plots == TRUE) {
  plot(x, sub=test)
}
cat('passed')
EpiModel:::test.epiNet(x)



# Cases -------------------------------------------------------------------

test <- 'simTrans: di.rate = 0 gives warning'
EpiModel:::mcat(test)

expect_error(
x <- epiNet.simTrans(
  est.vit, 
  type = 'SI',
  vital = TRUE,
  i.num = 10, 
  trans.rate = 0.5,
  act.rate = 2,
  b.rate = 0.01,
  ds.rate = 0.00,
  di.rate = 0.00,
  nsteps = 25,
  sims.per.nw = 1,
  verbose = FALSE,
  tea = FALSE) 
)
cat('passed')
EpiModel:::test.epiNet(x)

#####

test <- 'merge: works for vital sims saving nw stats'
EpiModel:::mcat(test)
nw <- network.initialize(n = 100, directed = FALSE)
est <- epiNet.est(nw,
                  formation = ~ edges,
                  dissolution = ~offset(edges),
                  target.stats = 25,
                  coef.diss = dissolution.coefs(~offset(edges), 10, 0),
                  save.stats = FALSE, verbose = FALSE)
x <- epiNet.simTrans(est, type = "SI",
                     i.num = 1, 
                     trans.rate = 0.9,
                     sims.per.nw = 2,
                     vital = TRUE,
                     b.rate = 0.01,
                     ds.rate = 0.01,
                     di.rate = 0.01,
                     nsteps = 10,
                     save.stats = TRUE,
                     stats.formula = ~edges + meandeg + degree(0) + concurrent,
                     verbose = FALSE) 
y <- epiNet.simTrans(est, type = "SI",
                     i.num = 1, 
                     trans.rate = 0.9,
                     sims.per.nw = 3, 
                     vital = TRUE,
                     b.rate = 0.01,
                     ds.rate = 0.01,
                     di.rate = 0.01,
                     nsteps = 10,
                     save.stats = TRUE,
                     stats.formula = ~edges + meandeg + degree(0) + concurrent,
                     verbose = FALSE) 
z <- merge(x, y)
expect_equal(length(z$stats), 5)
expect_true(all(sapply(z$stats, dim)[1,] == 10) & all(sapply(z$stats, dim)[2,] == 4))
cat('passed')

if (plots == TRUE) invisible(dev.off())
}
