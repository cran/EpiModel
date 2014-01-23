
##
## External testing cases for epiNet functions
##

library(EpiModel)
library(testthat)

# One mode SI, no vital ------------------------------------------

rm(list=ls())
num <- 500

nw <- network.initialize(num, directed=F)  
nw %v% 'race' <- rep(0:1, each=250)
nw

formation <- ~ edges + nodematch('race') + degree(0) + concurrent
target.stats <- c(225, 187, 180, 90)
dissolution <- ~ offset(edges) + offset(nodematch('race'))
dissolution <- ~ offset(edges)

durations <- c(20, 10)
coef.diss <- dissolution.coefs(dissolution, durations)
coef.diss
coef.diss <- coef.diss[1]

est <- epiNet.est(
  nw, 
  formation, 
  dissolution, 
  target.stats, 
  coef.diss, 
  edapprox = TRUE
)
est

par(mar=c(3,3,1,1), mgp=c(2,1,0))
plot(est, plots.joined=F)
plot(est, type='duration')

nwsims <- epiNet.simNet(est, nsteps=300, nsims=5)
nwsims
plot(nwsims, plots.joined=F)
plot(nwsims, type='duration')


x <- epiNet.simTrans(
  nwsims, 
  type = 'SI',
  vital = FALSE,
  trans.rate = 1, 
  i.num = 50,
  sims.per.nw = 1,
  verbose = TRUE,
  tea = FALSE
)

x
head(as.data.frame(x), 25)
head(as.data.frame(x, out='vals'), 25)
summary(x, time=20)
comp.plot(x, time=20, digits=0)
head(x$trans$sim1, 25)
EpiModel:::test.epiNet(x)

par(mar=c(3,3,1,1))
plot(x, y='i.num', popfrac=T) 
plot(x, sim.lines=T, sim.col=c('steelblue', 'firebrick'))

plot(x, type='network', at=10, col.inf=T)
plot(x, type='network', at=50, col.inf=T)
plot(x, type='network', at=50, sim=5, col.inf=T)



# One mode SIR, no vital ------------------------------------------


x
head(as.data.frame(x), 25)
head(as.data.frame(x, out='vals'), 25)
summary(x, time=20)
comp.plot(x, time=20)
head(x$trans$sim1, 25)
test.epiNet(x)

plot(x)
plot(x, popfrac=F, qnts=1)
plot(x, sim.lines=T, sim.col=c('steelblue', 'firebrick', 'seagreen'))
plot(x, type='network', at=10, col.inf=T)
plot(x, type='network', at=50, col.inf=T)



# One mode SIS, no vital ------------------------------------------

x <- epiNet.simTrans(
  nwsims, 
  type = 'SIS',
  vital = FALSE,
  trans.rate = 0.1, 
  rec.rate = 0.02,
  i.num = 10,
  i.rand = FALSE,
  nsteps = 10,
  sims.per.nw = 10
)

x
head(as.data.frame(x), 25)
head(as.data.frame(x, out='vals'), 25)
summary(x, time=20, comp.plot=T)
comp.plot(x, time=20, digits=1)
head(x$trans$sim1, 25)
test.epiNet(x)


plot(x)
plot(x, sim.lines=T, sim.col=c('steelblue', 'firebrick'))
plot(x, type='network', at=10, col.inf=T)
plot(x, type='network', at=50, col.inf=T)



# Two mode SI, no vital ------------------------------------------

rm(list=ls())

num.m1 <- 100 
num.m2 <- 100
nw <- network.initialize(num.m1 + num.m2, 
                         bipartite=num.m1, 
                         directed=FALSE)

deg.dist.m1 <- c(0.40, 0.55, 0.04, 0.01)
deg.dist.m2 <- c(0.48, 0.41, 0.08, 0.03)
bip.degdist.check(num.m1, num.m2, 
                  deg.dist.m1, deg.dist.m2)

formation <- ~ edges + b1degree(0:1) + b2degree(0:1)
target.stats <- c(330, 200, 275, 240, 205)/5
dissolution <- ~ offset(edges)

duration <- 25
coef.diss <- dissolution.coefs(dissolution, duration)
coef.diss

dx.stats <- ~ edges + b1degree(0:5) + b2degree(0:5)
est <- epiNet.est(
  nw, 
  formation, 
  dissolution, 
  target.stats, 
  coef.diss, 
  edapprox = TRUE, 
  save.stats = TRUE, 
  stats.formula = dx.stats)
est
plot(est, plots.joined = FALSE)
plot(est, type='duration')

nwsims <- epiNet.simNet(est, nsteps = 300, nsims = 1)
nwsims
plot(nwsims)
plot(nwsims, type='duration')

x <- epiNet.simTrans(
  nwsims, 
  type = 'SI',
  vital = FALSE,
  trans.rate = 0.2, 
  trans.rate.m2 = 0.1, 
  i.num = 10,
  sims.per.nw = 1,
  #nsteps = 10,
  verbose = TRUE
)

x
head(as.data.frame(x), 25)
head(as.data.frame(x, out='vals'), 25)
summary(x, time=20)
comp.plot(x, time=20, digits=0)
head(x$trans$sim1, 25)
test.epiNet(x)

plot(x)
plot(x, y='i.num', leg=T)

plot(x, type='network', at=2, col.inf=T, shp.bip='triangle')
plot(x, type='network', at=150, col.inf=T, shp.bip='square')
plot(x, type='network', at=50, sim=2, col.inf=T)



# Two mode SIR, no vital --------------------------------------------------

x <- epiNet.simTrans(
  nwsims, 
  type = 'SIR',
  vital = FALSE,
  trans.rate = 0.5, 
  trans.rate.m2 = 0.25, 
  rec.rate = 0.01,
  i.num = 10,
  sims.per.nw = 1,
  prog.plot = FALSE,
  tea = FALSE
)

x
head(as.data.frame(x), 25)
head(as.data.frame(x, out='vals',sim=9), 50)
summary(x, time=20)
comp.plot(x, time=20, digits=0)
head(x$trans$sim1, 25)
test.epiNet(x)

plot(x)
plot(x, popfrac=F, ylim=c(0, 300))
plot(x, y='i.num', sim.lines=T)
plot(x, y='i.num', leg=T)

plot(x, type='network', at=10, col.inf=T, shp.bip='triangle')
plot(x, type='network', at=300, col.inf=T, shp.bip='square')
plot(x, type='network', at=50, sim=2, col.inf=T)



# Two mode SIS, no vital --------------------------------------------------

x <- epiNet.simTrans(
  nwsims, 
  type = 'SIS',
  vital = FALSE,
  trans.rate = 0.2, 
  trans.rate.m2 = 0.1, 
  rec.rate = 0.02,
  i.num = 50,
  i.rand = FALSE,
  nsteps = 200,
  sims.per.nw = 20
)

x
head(as.data.frame(x), 25)
head(round(as.data.frame(x, out='sd'), 2), 25)
head(as.data.frame(x, out='vals'), 25)
summary(x, time=20)
comp.plot(x, time=20, digits=0)
head(x$trans$sim1, 25)
test.epiNet(x)

plot(x)
plot(x, y='i.num', sim.lines=T, sim.col='steelblue')
plot(x, y='i.num.m2', add=T, sim.lines=T, 
     sim.col='firebrick', mean.col='firebrick')

plot(x, type='network', at=10, col.inf=T, shp.bip='triangle')
plot(x, type='network', at=200, col.inf=T)
plot(x, type='network', at=200, sim=2, col.inf=T, shp.bip='square')



# One mode SI, vital ------------------------------------------------------

rm(list=ls())

num <- 500
nw <- network.initialize(num, directed=F)  
nw %v% 'race' <- rep(0:1, each=250)
nw

formation <- ~ edges + nodematch('race') + degree(0) + concurrent
target.stats <- c(225, 187, 180, 90)
dissolution <- ~ offset(edges) 

durations <- 25
coef.diss <- dissolution.coefs(dissolution, durations, d.rate = 0.0025)
coef.diss

est <- epiNet.est(
  nw, 
  formation, 
  dissolution, 
  target.stats, 
  coef.diss, 
  edapprox=TRUE
)
est
plot(est)
plot(est, type='duration')

dx.stats <- ~ edges + degree(0:5)
x <- epiNet.simTrans(
  est,
  type = 'SI',
  vital = TRUE,
  trans.rate = 0.2,
  i.num = 50,
  b.rate = 0.0025,
  ds.rate = 0.0025,
  di.rate = 0.0025,
  save.stats = TRUE, 
  stats.formula = dx.stats,
  nsteps = 100,
  sims.per.nw = 3,
  verbose = T,
  save.statmat = T
)

x
head(as.data.frame(x), 25)
head(as.data.frame(x, out='vals'), 25)
summary(x, time=20)
comp.plot(x, time=25, digits=2)
head(x$trans$sim1, 25)
EpiModel:::test.epiNet(x)

plot(x)
plot(x, type='network', at=10, col.inf=T)
plot(x, type='network', at=25, col.inf=T)
plot(x, type='network', at=100, sim=4, col.inf=T)



# One mode SIR, vital -----------------------------------------------------

trans.rate <- 0.15
b.rate <- 0.0025
ds.rate <- 0.0025
di.rate <- 0.0025
dr.rate <- 0.0025
rec.rate <- 0.1

dx.stats <- ~ edges + degree(0:5)
x <- epiNet.simTrans(
  est,
  type = 'SIR',
  vital = TRUE,
  nsteps = 25,
  trans.rate = trans.rate,
  rec.rate = rec.rate,
  i.num = 100,
  b.rate = b.rate,
  ds.rate = ds.rate,
  di.rate = di.rate,
  dr.rate = dr.rate,
  save.stats = TRUE, 
  stats.formula = dx.stats,
  sims.per.nw = 3,
  save.statmat = F,
  verbose = T)

x
head(as.data.frame(x), 25)
head(as.data.frame(x, out='vals'), 25)
summary(x, time=20)
summary(x, time=100)
comp.plot(x, time=50, digits=2)
head(x$trans$sim1, 25)
test.epiNet(x)

plot(x)
plot(x, popfrac=F)
plot(x, type='network', at=10, col.inf=T)
plot(x, type='network', at=100, col.inf=T)
plot(x, type='network', at=100, sim=4, col.inf=T)


# One mode SIS, vital -----------------------------------------------------

trans.rate <- 0.5
b.rate <- 0.0025
ds.rate <- 0.0025
di.rate <- 0.0025
rec.rate <- 0.04

x <- epiNet.simTrans(
  est,
  type = 'SIS',
  vital = TRUE,
  nsteps = 25,
  trans.rate = trans.rate,
  rec.rate = rec.rate,
  i.num = 1,
  b.rate = b.rate,
  ds.rate = ds.rate,
  di.rate = di.rate,
  sims.per.nw = 3,
  tea = FALSE,
  verbose = TRUE)

x
head(as.data.frame(x), 25)
head(as.data.frame(x, out='vals'), 25)
summary(x, time=20)
summary(x, time=100)
comp.plot(x, time=50, digits=2)
head(x$trans$sim1, 25)
EpiModel:::test.epiNet(x)

plot(x)
plot(x, popfrac=F)
plot(x, y='is.flow', popfrac=F, leg=T)
plot(x, type='network', at=10, col.inf=T)
plot(x, type='network', at=100, col.inf=T)
plot(x, type='network', at=100, sim=4, col.inf=T)


# Two mode SI, vital ------------------------------------------------------

rm(list=ls())

num.m1 <- 250
num.m2 <- 250
nw <- network.initialize(num.m1 + num.m2, bipartite=num.m1, directed=F)

deg.dist.m1 <- c(0.40, 0.55, 0.04, 0.01)
deg.dist.m2 <- c(0.48, 0.41, 0.08, 0.03)
bip.degdist.check(num.m1, num.m2, 
                  deg.dist.m1, deg.dist.m2)

formation <- ~ edges + b1degree(0:1) + b2degree(0:1)
target.stats <- c(330, 200, 275, 240, 205)/2
dissolution <- ~ offset(edges)

duration <- 25
coef.diss <- dissolution.coefs(dissolution, duration, d.rate = 0.01)
coef.diss

dx.stats <- ~ edges + b1degree(0:5) + b2degree(0:5)
est <- epiNet.est(
  nw, 
  formation, 
  dissolution, 
  target.stats, 
  coef.diss, 
  edapprox=TRUE, 
  save.stats = TRUE, 
  stats.formula = dx.stats)
est
plot(est, plots.joined=F)
plot(est, type='duration')

trans.rate <- 0.1
trans.rate.m2 <- trans.rate/2
b.rate <- 0.01
ds.rate <- 0.01
di.rate <- 0.01

x <- epiNet.simTrans(
  est,
  type = 'SI',
  vital = TRUE,
  i.num = 10,
  trans.rate = trans.rate,
  trans.rate.m2 = trans.rate.m2,
  b.rate = b.rate*2,
  ds.rate = ds.rate,
  di.rate = di.rate,
  nsteps = 200,
  sims.per.nw = 1,
  tea = F,
  verbose = TRUE
)
x
head(as.data.frame(x), 25)
head(as.data.frame(x, out='vals'), 25)
summary(x, time=20)
comp.plot(x, time=50, digits=2)
head(x$trans$sim1, 25)
test.epiNet(x)

plot(x)
plot(x, type='network', at=1, col.inf=T)
plot(x, type='network', at=1, col.inf=T, shp.bip='square')
plot(x, type='network', at=100, sim=4, col.inf=T)


# Two mode SIR, vital -----------------------------------------------------

trans.rate <- 0.1
trans.rate.m2 <- trans.rate/2
b.rate <- 0.01
ds.rate <- 0.01
di.rate <- 0.01
rec.rate <- 1/30

x <- epiNet.simTrans(
  est,
  type = 'SIR',
  vital = TRUE,
  i.num = 100,
  trans.rate = trans.rate,
  trans.rate.m2 = trans.rate.m2,
  rec.rate = rec.rate,
  b.rate = b.rate,
  ds.rate = ds.rate,
  di.rate = di.rate,
  dr.rate = ds.rate,
  nsteps = 100,
  sims.per.nw = 1,
  verbose = TRUE
)
x
head(as.data.frame(x), 25)
head(as.data.frame(x, out='vals'), 25)
summary(x, time=20)
comp.plot(x, time=50, digits=2)
head(x$trans$sim1, 25)
EpiModel:::test.epiNet(x)

plot(x)
plot(x, type='network', at=10, col.inf=T)
plot(x, type='network', at=100, col.inf=T)
plot(x, type='network', at=100, sim=4, col.inf=T)


# Two mode SIS, vital -----------------------------------------------------

trans.rate <- 0.1
trans.rate.m2 <- trans.rate/2
b.rate <- 0.01
ds.rate <- 0.01
di.rate <- 0.01
rec.rate <- 1/30

x <- epiNet.simTrans(
  est,
  type = 'SIS',
  vital = TRUE,
  i.num = 10,
  trans.rate = trans.rate,
  trans.rate.m2 = trans.rate.m2,
  rec.rate = rec.rate,
  b.rate = b.rate,
  ds.rate = ds.rate,
  di.rate = di.rate,
  nsteps = 10,
  sims.per.nw = 5,
  verbose = TRUE
)
x
head(as.data.frame(x), 25)
head(as.data.frame(x, out='vals'), 25)
summary(x, time=20)
comp.plot(x, time=50, digits=2)
head(x$trans$sim1, 25)
test.epiNet(x)

plot(x)
plot(x, type='network', at=10, col.inf=T, shp.bip='triangle')
plot(x, type='network', at=100, col.inf=T)
plot(x, type='network', at=100, sim=4, col.inf=T)



# Two mode SI, no vital, set IDs ------------------------------------------

rm(list=ls())

num.m1 <- 250
num.m2 <- 250
nw <- network.initialize(num.m1 + num.m2, bipartite=num.m1, directed=F)

deg.dist.m1 <- c(0.40, 0.55, 0.04, 0.01)
deg.dist.m2 <- c(0.48, 0.41, 0.08, 0.03)
bip.degdist.check(num.m1, num.m2, 
                  deg.dist.m1, deg.dist.m2)

formation <- ~ edges + b1degree(0:1) + b2degree(0:1)
target.stats <- c(330, 200, 275, 240, 205)/2
dissolution <- ~ offset(edges)
constraints <- ~ bd(maxout=3)

( coef.diss <- dissolution.coefs(dissolution, duration = 25) )

dx.stats <- ~ edges + b1degree(0:5) + b2degree(0:5)
est <- epiNet.est(
  nw, 
  formation, 
  dissolution, 
  target.stats, 
  coef.diss, 
  constraints, 
  stats.formula = dx.stats)
est
plot(est, plots.joined = FALSE)
plot(est, type='duration', dx.start=100)

trans.rate <- 0.1
trans.rate.m2 <- trans.rate/2

nwsims <- epiNet.simNet(est, nsteps=100, nsims=10)

x <- epiNet.simTrans(
  nwsims, 
  type = 'SI',
  vital = FALSE,
  trans.rate = trans.rate, 
  trans.rate.m2 = trans.rate.m2, 
  i.ids = c(1:10, 251),
  sims.per.nw = 1,
  nsteps = 50)

sim.nw1 <- x$network$sim1
sim.nw10 <- x$network$sim10

get.stat(sim.nw1, stat=1, at=1, out='all')
get.stat(sim.nw10, stat=1, at=1, out='all')


x
head(as.data.frame(x), 25)
head(as.data.frame(x, out='vals'), 25)
summary(x, time=20)
comp.plot(x, time=20, digits=0)
head(x$trans$sim1, 25)
test.epiNet(x)

plot(x)
plot(x, y='i.num', leg=T)

plot(x, type='network', at=10, col.inf=T, shp.bip='triangle')
plot(x, type='network', at=50, col.inf=T, shp.bip='square')
plot(x, type='network', at=50, sim=2, col.inf=T)
