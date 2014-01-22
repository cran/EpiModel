
## ----setup, include=FALSE, cache=FALSE-----------------------------------
library(EpiModel)
library(knitr)
opts_chunk$set(cache=FALSE, comment=NA, tidy=FALSE)
par(mfrow=c(1,1),mar=c(3,3,1,1), mgp=c(2,1,0))


## ------------------------------------------------------------------------
dissolution <- ~offset(edges)
duration <- 25
coefs <- dissolution.coefs(dissolution, duration)
coefs


## ------------------------------------------------------------------------
dissolution.coefs(dissolution, duration, d.rate = 0.001)


## ------------------------------------------------------------------------
# An imbalanced distribution          
bip.degdist.check(num.m1 = 500, 
                  num.m2 = 500, 
                  deg.dist.m1 = c(0.40, 0.55, 0.03, 0.02),
                  deg.dist.m2 = c(0.48, 0.41, 0.08, 0.03))


## ------------------------------------------------------------------------
# A balanced distribution
targets <- bip.degdist.check(
                  num.m1 = 500, 
                  num.m2 = 500, 
                  deg.dist.m1 = c(0.40, 0.55, 0.04, 0.01),
                  deg.dist.m2 = c(0.48, 0.41, 0.08, 0.03))
targets


