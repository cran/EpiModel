
## ----setup, include=FALSE------------------------------------------------
require(knitr)
require(EpiModel)
opts_chunk$set(dev="pdf", fig.align="center", cache=TRUE, comment=NA, 
               message=FALSE, fig.width=8, out.width="0.9\\linewidth", 
               fig.height=5, fig.pos="ht!", tidy=FALSE, highlight=TRUE)
options(replace.assign=TRUE,width=70)
par(mar=c(3.5,3,1,1), mgp=c(2,1,0))
version <- packageDescription('EpiModel')$Version


## ----help, eval=FALSE----------------------------------------------------
## help(package = "EpiModel")
## ?EpiModel
## browseVignettes(package = "EpiModel")


## ----dcmSi---------------------------------------------------------------
mod <- epiDCM(type = "SI", s.num = 500, i.num = 1, 
              trans.rate = 0.2, act.rate = 0.25, nsteps = 500)


## ----dcmSiPrint----------------------------------------------------------
mod


## ----dcmSiHead-----------------------------------------------------------
head(mod$i.num)
head(mod$si.flow)


## ----dcmSiPlot-----------------------------------------------------------
plot(mod)


## ----dcmSiSumm-----------------------------------------------------------
summary(mod, time = 150)


## ----dcmSir--------------------------------------------------------------
mod <- epiDCM(type = "SIR", s.num = 1000, i.num = 1, r.num = 0, 
              trans.rate = 0.2, act.rate = 1, rec.rate = 1/20,
              b.rate = 1/95, ds.rate = 1/100, di.rate = 1/80, 
              dr.rate = 1/100, nsteps = 500, dt = 0.5)


## ----dcmSirPlot, fig.height=4--------------------------------------------
par(mar = c(3.2,3,2,1), mgp = c(2,1,0), mfrow = c(1,2))
plot(mod, popfrac = FALSE, alpha = 0.5, 
     lwd = 4, main = "Compartment Sizes")
plot(mod, y = "si.flow", lwd = 4, col = "firebrick", 
     main = "Disease Incidence", leg = "n")


## ----dcmSirCPlot---------------------------------------------------------
par(mfrow = c(1,1))
comp.plot(mod, time = 50, digits = 1)


## ----dcmSis--------------------------------------------------------------
mod <- epiDCM(type = "SIS", s.num = 500, i.num = 1,
              trans.rate = 0.2, act.rate = seq(0.25, 0.5, 0.05),
              rec.rate = 1/50, nsteps = 350)


## ----dcmSisPrint---------------------------------------------------------
mod


## ----dcmSisHead----------------------------------------------------------
head(as.data.frame(mod, run = 5))


## ----dcmSisPlot, fig.height=4--------------------------------------------
par(mfrow = c(1,2), mar = c(3.2,3,2.5,1))
plot(mod, alpha = 1, main = "Disease Prevalence")
plot(mod, y = "si.flow", col = "Greens", alpha = 0.8,
     main = "Disease Incidence")


## ----dcmSisPlotOpts, eval=FALSE------------------------------------------
## plot(mod, col = "black")
## plot(mod, col = 1:6)
## plot(mod, col = c("black", "red", "blue", "green", "purple", "pink"))
## plot(mod, col = rainbow(6))


## ----dcmSisSensOpts, eval=FALSE------------------------------------------
## act.rates <- c(0.2, 0.2, 0.4, 0.4, 0.6, 0.6)
## trans.rates <- c(0.1, 0.2, 0.1, 0.2, 0.1, 0.2)
## mod <- epiDCM(type = "SIS", s.num = 500, i.num = 1,
##               trans.rate = trans.rates, act.rate = act.rates,
##               rec.rate=1/50, nsteps = 350)


## ----dcmSi2g, results="hide"---------------------------------------------
mod <- epiDCM(type = "SI", groups = 2, 
              s.num = 500, i.num = 1, 
              s.num.g2 = 500, i.num.g2 = 0, 
              trans.rate = 0.4,  trans.rate.g2 = 0.1, 
              act.rate = 0.25, balance = "g1",
              b.rate = 1/100, b.rate.g2 = NA,
              ds.rate = 1/100, ds.rate.g2 = 1/100,
              di.rate = 1/50, di.rate.g2 = 1/50,
              nsteps = 500)


## ----dcmSi2gPrint--------------------------------------------------------
mod


## ----dcmSi2gPlot---------------------------------------------------------
plot(mod)


## ----dcmGui, eval=FALSE--------------------------------------------------
## gui.epiDCM()


## ----icmSi, results="hide"-----------------------------------------------
mod <- epiICM(type = "SI", s.num = 500, i.num = 1, 
              trans.rate = 0.2, act.rate = 0.25, 
              nsims = 10, nsteps = 300)


## ----icmSiprint----------------------------------------------------------
mod


## ----icmSiSumm-----------------------------------------------------------
summary(mod, time = 125)


## ----icmSIAsDf-----------------------------------------------------------
head(as.data.frame(mod, out='mean'))


## ----icmSiPlot-----------------------------------------------------------
plot(mod, sim.lines = TRUE, sim.col = c("steelblue", "firebrick"))


## ----icmSirDet-----------------------------------------------------------
det <- epiDCM(type = "SIR", s.num = 1000, i.num = 100,
              trans.rate = 0.2, act.rate = 0.8, rec.rate = 1/50,
              b.rate = 1/100, ds.rate = 1/100, di.rate = 1/90,
              dr.rate = 1/100, nsteps = 300)


## ----icmSirSto, results="hide"-------------------------------------------
sim <- epiICM(type = "SIR", s.num = 1000, i.num = 100,
              trans.rate = 0.2, act.rate = 0.8, rec.rate = 1/50,
              b.rate = 1/100, ds.rate = 1/100, di.rate = 1/90,
              dr.rate = 1/100, nsteps = 300, nsims = 10)


## ----icmSirSto2, results="hide"------------------------------------------
sim2 <- epiICM(type = "SIR", s.num = 1000, i.num = 100,
               trans.rate = 0.2, act.rate = 0.8, rec.rate = 1/50,
               b.rate = 1/100, ds.rate = 1/100, di.rate = 1/90,
               dr.rate = 1/100, nsteps = 300, nsims = 10,
               b.rand = FALSE, d.rand = FALSE)


## ----icmSirPlot----------------------------------------------------------
plot(det, alpha = 0.75, lwd = 4, main = "DCM and ICM Comparison")
plot(sim, qnts = FALSE, add = TRUE, mean.lty = c(2,2,2), leg = FALSE)
plot(sim2, qnts = FALSE, add = TRUE, mean.lty = c(3,3,3), leg = FALSE)


## ----icmSirPlot2, fig.height=4-------------------------------------------
par(mfrow = c(1,2), mar = c(3,3,2,1), mgp = c(2,1,0))
plot(sim, y = "di.flow", mean.line = FALSE, 
     sim.lines = TRUE, sim.alpha = 0.5, 
     ylim = c(0, max(sim$di.flow)), 
     main = "di.flow: Full Stochastic Model")
plot(sim2, y = "di.flow", mean.line = FALSE, 
     sim.lines = TRUE, sim.alpha = 0.5, 
     ylim = c(0, max(sim$di.flow)), 
     main = "di.flow: Limited Stochastic Model")


## ----icmSirAdf-----------------------------------------------------------
icm.compare <- rbind(round(as.data.frame(sim, out = "sd")[50,], 2),
                     round(as.data.frame(sim2, out = "sd")[50,], 2))
row.names(icm.compare) <- c("full", "lim")
icm.compare


## ----icmSisSim, results="hide"-------------------------------------------
set.seed(12345)
sim <- epiICM(type = "SIS", groups = 2,
              s.num = 500, i.num = 1,
              s.num.g2 = 500, i.num.g2 = 1,
              trans.rate = 0.2, trans.rate.g2 = 0.1,
              act.rate = 0.5, balance = "g1",
              rec.rate = 1/25, rec.rate.g2 = 1/50,
              b.rate = 1/100, b.rate.g2 = NA,
              ds.rate = 1/100, ds.rate.g2 = 1/100,
              di.rate = 1/90, di.rate.g2 = 1/90,
              nsteps = 500, nsims = 10)


## ----icmSISPlot1---------------------------------------------------------
par(mfrow = c(1,1))
plot(sim)


## ----icmSisPlot2---------------------------------------------------------
plot(sim, y = c("i.num", "i.num.g2"), mean.lwd=3,
     sim.lines = TRUE, sim.col = c('steelblue', 'firebrick'), 
     main = "Disease Prevalence: Means and Individual Simulations",
     leg = TRUE)


## ----netEstParams1-------------------------------------------------------
nw <- network.initialize(n = 100, directed = FALSE)  
nw %v% "race" <- rep(0:1, each = 50)


## ----netEstParams2-------------------------------------------------------
formation <- ~ edges + nodematch("race") + degree(0) + concurrent
target.stats <- c(45, 37.4, 36, 18)
dissolution <- ~ offset(edges)


## ----netEstParams3-------------------------------------------------------
duration <- 20
coef.diss <- dissolution.coefs(dissolution, duration)
coef.diss


## ----netEst, results="hide", warning=FALSE-------------------------------
est1 <- epiNet.est(nw, 
                   formation, 
                   dissolution, 
                   target.stats, 
                   coef.diss)


## ----netEstPrint---------------------------------------------------------
est1


## ----netEstPlot1---------------------------------------------------------
plot(est1)


## ----netEstPlot2---------------------------------------------------------
plot(est1, plots.joined = FALSE, dx.start = 100, dx.end = 200)


## ----netSim, results="hide"----------------------------------------------
nwsims <- epiNet.simNet(est1, nsteps = 50, nsims = 5)


## ----netSimPrint---------------------------------------------------------
nwsims


## ----netSimPlot----------------------------------------------------------
plot(nwsims, plots.joined = FALSE)


## ----netSimPlot2---------------------------------------------------------
plot(nwsims, type = "duration")


## ----netIndSimTrans, results="hide"--------------------------------------
sim1 <- epiNet.simTrans(nwsims, 
                        type = "SIS",
                        sims.per.nw = 1,
                        trans.rate = 0.5,
                        act.rate = 3,
                        rec.rate = 0.1,
                        i.num = 10)


## ----netIndSimTransPrint-------------------------------------------------
sim1


## ----netIndSimTransSumm--------------------------------------------------
summary(sim1, time = 25, comp.plot = TRUE)


## ----netIndSimTransHead--------------------------------------------------
head(as.data.frame(sim1, out = 'vals', sim = 2))


## ----netIndSimTransNw----------------------------------------------------
( nw1 <- sim1$network$sim1 )


## ----netIndSimTransTrans-------------------------------------------------
trans1 <- sim1$trans$sim1
head(trans1, 10)


## ----netIndSimTransPlot1-------------------------------------------------
plot(sim1)


## ----netIndSimTransPlot2-------------------------------------------------
plot(sim1, y = c("si.flow", "is.flow"), leg = TRUE)


## ----netIndSimTransPlot3-------------------------------------------------
par(mfrow = c(1,2), mar = c(0,0,2,0))
plot(sim1, type = "network", at = 1, col.inf = TRUE, 
     zeromarg = FALSE, main = "Prevalence at t1")
plot(sim1, type = "network", at = 50, col.inf = TRUE, 
     zeromarg = FALSE, main = "Prevalence at t50")


## ----netIndSimTransPlot4-------------------------------------------------
time <- 25
lsim <- which.min(sim1$i.num[time, ])
hsim <- which.max(sim1$i.num[time, ])
lprev <- sim1$i.num[time, lsim]
hprev <- sim1$i.num[time, hsim]
par(mfrow = c(1,2), mar = c(0,0,2,0))
plot(sim1, type = "network", at = 25, sim = lsim, 
     col.inf = TRUE, zeromarg = FALSE, 
     main = paste("Sim", lsim, ", prev=", lprev, sep=""))
plot(sim1, type = "network", at = 25, sim = hsim, 
     col.inf = TRUE, zeromarg = FALSE, 
     main = paste("Sim", hsim, ", prev=", hprev, sep=""))


## ----netDepEstParams1----------------------------------------------------
num.m1 <- 50 
num.m2 <- 50
nw <- network.initialize(num.m1 + num.m2, 
                         bipartite = num.m1, directed = FALSE)


## ----netDepEstParams2----------------------------------------------------
deg.dist.m1 <- c(0.40, 0.55, 0.04, 0.01)
deg.dist.m2 <- c(0.48, 0.41, 0.08, 0.03)
par(mar=c(3,3,2,1))
barplot(cbind(deg.dist.m1, deg.dist.m2), beside = TRUE, 
        legend.text = paste("deg", 0:3, sep=""), ylim = c(0,0.6))


## ----netDepEstParams2b---------------------------------------------------
bip.degdist.check(num.m1, num.m2, 
                  deg.dist.m1, deg.dist.m2)


## ----netDepParams3-------------------------------------------------------
formation <- ~ edges + b1degree(0:1) + b2degree(0:1)
target.stats <- c(33, 20, 27.5, 24, 20.5)
dissolution <- ~ offset(edges)
coef.diss <- dissolution.coefs(dissolution, duration=25, d.rate=0.01)
coef.diss


## ----netDepEst, results="hide", warning=FALSE----------------------------
dx.stats <- ~ edges + b1degree(0:5) + b2degree(0:5)
set.seed(12345)
est2 <- epiNet.est(nw, 
                   formation, 
                   dissolution, 
                   target.stats, 
                   coef.diss, 
                   stats.formula = dx.stats)


## ----netDepEstPrint------------------------------------------------------
est2


## ----netDepEstPlot-------------------------------------------------------
plot(est2)


## ----netDepSimTransParams------------------------------------------------
i.num <- 10
trans.rate <- 0.5
trans.rate.m2 <- 0.1
b.rate <- 2/100
ds.rate <- 1/100
di.rate <- 1/90


## ----netDepSimTrans, results="hide"--------------------------------------
sim2 <- epiNet.simTrans(est2,
                        type = "SI",
                        vital = TRUE,
                        i.num = i.num,
                        trans.rate = trans.rate,
                        trans.rate.m2 = trans.rate.m2,
                        b.rate = b.rate,
                        ds.rate = ds.rate,
                        di.rate = di.rate,
                        sims.per.nw = 3,
                        nsteps = 50)


## ----netDepSimTransPrint-------------------------------------------------
sim2


## ----netDeptSimTransTrans------------------------------------------------
trans1 <- sim2$trans$sim1
head(trans1)


## ----netDepSimTransPlot1-------------------------------------------------
plot(sim2, popfrac=FALSE)


## ----netDepSimTransPlot3-------------------------------------------------
par(mfrow=c(1,2))
plot(sim2, type = "network", at = 1, col.inf = TRUE, shp.bip = "triangle")
plot(sim2, type = "network", at = 50, col.inf = TRUE, shp.bip = "square")


## ----ndtvLoad, eval=FALSE------------------------------------------------
## library(ndtv)


## ----ndtvRecode, eval=FALSE----------------------------------------------
## nw <- sim2$network$sim1
## nw <- colorTEA(nw)


## ----ndtvCoords, eval=FALSE----------------------------------------------
## slice.par <- list(start = 1, end = 50, interval = 1,
##                   aggregate.dur = 1, rule = "any")
## compute.animation(nw, slice.par = slice.par,
##                   animation.mode = "MDSJ")


## ----ndtvRender, eval=FALSE----------------------------------------------
## render.par=list(tween.frames = 10,
##                 show.time = FALSE)
## plot.par=list(mar = c(0,0,0,0))
## render.animation(nw,
##                  render.par = render.par,
##                  plot.par = plot.par,
##                  vertex.cex = 0.9,
##                  vertex.col = "ndtvcol",
##                  edge.col = "darkgrey",
##                  vertex.border = "lightgrey",
##                  displaylabels = FALSE)


## ----ndtvSaveGif, eval=FALSE---------------------------------------------
## saveGIF(ani.replay(),
##         ani.width = 600,
##         ani.height = 600,
##         outdir = getwd())


## ----ndtvSaveVideo, eval=FALSE-------------------------------------------
## saveVideo(ani.replay(),
##           video.name = "EpiModelndtv.mp4",
##           other.opts = "-b 5000k",
##           clean = TRUE,
##           ani.width = 1200,
##           ani.height = 1200)


## ----eval=FALSE----------------------------------------------------------
## help('epiNetModules')


