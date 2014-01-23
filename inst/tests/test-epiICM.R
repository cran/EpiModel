
##
## Testing Series for epiICM
##


run.tests.epiICM <- function(plots=FALSE, seed=FALSE) {

library(EpiModel)
library(testthat)

if (plots == TRUE) pdf('EpiModel/epiICM.pdf', h=5, w=10)
par(mar=c(5,5,3,1), mgp=c(3,2,1))

cat('\n===============================')
cat('\nTesting Series for epiICM')
cat('\n===============================')

# SI Models ---------------------------------------------------------------

test <- 'SI model, 1 group, no vital, 1 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = 'SI', 
    s.num = 500, 
    i.num = 1,
    trans.rate = 0.2, 
    act.rate = 0.25, 
    nsteps = 500, 
    nsims = 1, 
    verbose = F
  )
  expect_equal(max(x$i.num), 501)
  expect_equal(max(x$s.num+x$i.num), 501)
  expect_is(as.data.frame(x), 'data.frame')
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SI model, 1 group, no vital, 5 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = 'SI', 
    s.num = 500, 
    i.num = 1,
    trans.rate = 0.2, 
    act.rate = 0.25, 
    nsteps = 500, 
    nsims = 5, 
    verbose = F
  )
  expect_equal(max(x$i.num), 501)
  expect_equal(max(x$s.num+x$i.num), 501)
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SI model, 2 groups, no vital, 1 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = "SI", 
    groups = 2, 
    s.num = 500, 
    i.num = 1, 
    trans.rate = 0.2, 
    act.rate = 0.25, 
    s.num.g2 = 500, 
    i.num.g2 = 0, 
    trans.rate.g2 = 0.1, 
    balance = "g1",
    nsteps = 500, 
    nsims = 1, 
    verbose = F
  )
  expect_equal(max(x$i.num), 501)
  expect_equal(max(x$s.num+x$i.num), 501)
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SI model, 2 groups, no vital, 5 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = "SI", 
    groups = 2, 
    s.num = 500, 
    i.num = 1, 
    trans.rate = 0.2, 
    act.rate = 0.25, 
    s.num.g2 = 500, 
    i.num.g2 = 0, 
    trans.rate.g2 = 0.1, 
    balance = "g1",
    nsteps = 500, 
    nsims = 5, 
    verbose = F
  )
  expect_equal(max(x$i.num), 501)
  expect_equal(max(x$s.num+x$i.num), 501)
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SI model, 1 group, vital, 1 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = 'SI', 
    s.num = 500, 
    i.num = 1,
    trans.rate = 0.2, 
    act.rate = 0.25,
    b.rate = 1/100,
    ds.rate = 1/100,
    di.rate = 1/90,
    nsteps = 500, 
    nsims = 1, 
    verbose = F
  )
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SI model, 1 group, vital, 5 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = 'SI', 
    s.num = 500, 
    i.num = 1,
    trans.rate = 0.2, 
    act.rate = 0.25,
    b.rate = 1/100,
    ds.rate = 1/100,
    di.rate = 1/90,
    nsteps = 500, 
    nsims = 5, 
    verbose = F
  )
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})


test <- 'SI model, 2 groups, vital, 1 sim' 
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = "SI", 
    groups = 2, 
    s.num = 500, 
    i.num = 1, 
    trans.rate = 0.2, 
    act.rate = 0.25, 
    s.num.g2 = 500, 
    i.num.g2 = 0, 
    trans.rate.g2 = 0.1, 
    balance = "g1",
    b.rate = 1/100,
    b.rate.g2 = NA,
    ds.rate = 1/100,
    ds.rate.g2 = 1/100,
    di.rate = 1/90,
    di.rate.g2 = 1/90,
    nsteps = 500, 
    nsims = 1, 
    verbose = F
  )
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SI model, 2 groups, vital, 5 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = "SI", 
    groups = 2, 
    s.num = 500, 
    i.num = 1, 
    trans.rate = 0.2, 
    act.rate = 0.25, 
    s.num.g2 = 500, 
    i.num.g2 = 0, 
    trans.rate.g2 = 0.1, 
    balance = "g1",
    b.rate = 1/100,
    b.rate.g2 = NA,
    ds.rate = 1/100,
    ds.rate.g2 = 1/100,
    di.rate = 1/90,
    di.rate.g2 = 1/90,
    nsteps = 500, 
    nsims = 5, 
    verbose = F
  )
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})


# SIR Models --------------------------------------------------------------

test <- 'SIR model, 1 group, no vital, 1 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = 'SIR', 
    s.num = 500, 
    i.num = 1,
    trans.rate = 0.2, 
    act.rate = 0.25, 
    rec.rate = 1/50,
    nsteps = 500, 
    nsims = 1, 
    verbose = F
  )
  expect_equal(max(x$s.num+x$i.num+x$r.num), 501)
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SIR model, 1 group, no vital, 5 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = 'SIR', 
    s.num = 500, 
    i.num = 1,
    trans.rate = 0.2, 
    act.rate = 0.25, 
    rec.rate = 1/50,
    nsteps = 500, 
    nsims = 5, 
    verbose = F
  )
  expect_equal(max(x$s.num+x$i.num+x$r.num), 501)
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)    
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SIR model, 2 groups, no vital, 1 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = "SIR", 
    groups = 2, 
    s.num = 500, 
    i.num = 1, 
    s.num.g2 = 500, 
    i.num.g2 = 0, 
    trans.rate = 0.2, 
    act.rate = 0.25, 
    trans.rate.g2 = 0.1, 
    balance = "g1",
    rec.rate = 1/100,
    rec.rate.g2 = 1/100,
    nsims = 1,
    nsteps = 1000, 
    verbose = F
  )
  expect_equal(max(x$s.num+x$i.num+x$r.num+
                   x$s.num.g2+x$i.num.g2+x$r.num.g2), 1001)
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SIR model, 2 groups, no vital, 5 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = "SIR", 
    groups = 2, 
    s.num = 500, 
    i.num = 1, 
    s.num.g2 = 500, 
    i.num.g2 = 0, 
    trans.rate = 0.2, 
    act.rate = 0.25, 
    trans.rate.g2 = 0.1, 
    balance = "g1",
    rec.rate = 1/100,
    rec.rate.g2 = 1/100,
    nsims = 5,
    nsteps = 1000, 
    verbose = F
  )
  expect_equal(max(x$s.num+x$i.num+x$r.num+
                   x$s.num.g2+x$i.num.g2+x$r.num.g2), 1001)
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SIR model, 1 group, vital, 1 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = 'SIR', 
    groups = 1,
    s.num = 500, 
    i.num = 1,
    trans.rate = 0.2, 
    act.rate = 3,
    rec.rate = 1/50,
    b.rate = 1/100,
    ds.rate = 1/100,
    di.rate = 1/90,
    dr.rate = 1/100,
    nsteps = 500, 
    nsims = 1, 
    verbose = F
  )
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SIR model, 1 group, vital, 5 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = 'SIR', 
    groups = 1,
    s.num = 500, 
    i.num = 1,
    trans.rate = 0.2, 
    act.rate = 3,
    rec.rate = 1/50,
    b.rate = 1/100,
    ds.rate = 1/100,
    di.rate = 1/90,
    dr.rate = 1/100,
    nsteps = 500, 
    nsims = 5, 
    verbose = F
  )
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SIR model, 2 group, vital, 1 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = 'SIR', 
    groups = 2,
    s.num = 500, 
    i.num = 1,
    s.num.g2 = 500,
    i.num.g2 = 1,
    trans.rate = 0.2, 
    trans.rate.g2 = 0.1,
    act.rate = 3,
    balance = 'g1',
    rec.rate = 1/50,
    rec.rate.g2 = 1/50,
    b.rate = 1/100,
    b.rate.g2 = NA,
    ds.rate = 1/100,
    ds.rate.g2 = 1/100,
    di.rate = 1/90,
    di.rate.g2 = 1/90,
    dr.rate = 1/100,
    dr.rate.g2 = 1/100,
    nsteps = 500, 
    nsims = 1, 
    verbose = F
  )
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SIR model, 2 group, vital, 5 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = 'SIR', 
    groups = 2,
    s.num = 500, 
    i.num = 1,
    s.num.g2 = 500,
    i.num.g2 = 1,
    trans.rate = 0.2, 
    trans.rate.g2 = 0.1,
    act.rate = 3,
    balance = 'g1',
    rec.rate = 1/50,
    rec.rate.g2 = 1/50,
    b.rate = 1/100,
    b.rate.g2 = NA,
    ds.rate = 1/100,
    ds.rate.g2 = 1/100,
    di.rate = 1/90,
    di.rate.g2 = 1/90,
    dr.rate = 1/100,
    dr.rate.g2 = 1/100,
    nsteps = 500, 
    nsims = 5, 
    verbose = F
  )
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})


# SIS Models --------------------------------------------------------------

test <- 'SIS model, 1 group, no vital, 1 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = 'SIS', 
    s.num = 500, 
    i.num = 1,
    trans.rate = 0.2, 
    act.rate = 0.25, 
    rec.rate = 1/50,
    nsteps = 500, 
    nsims = 1, 
    verbose = F
  )
  expect_equal(max(x$s.num+x$i.num), 501)
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SIS model, 1 group, no vital, 5 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = 'SIS', 
    s.num = 500, 
    i.num = 1,
    trans.rate = 0.2, 
    act.rate = 0.25, 
    rec.rate = 1/50,
    nsteps = 500, 
    nsims = 5, 
    verbose = F
  )
  expect_equal(max(x$s.num+x$i.num), 501)
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SIS model, 2 groups, no vital, 1 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = "SIS", 
    groups = 2, 
    s.num = 500, 
    i.num = 1, 
    s.num.g2 = 500, 
    i.num.g2 = 0, 
    trans.rate = 0.2, 
    act.rate = 0.25, 
    trans.rate.g2 = 0.1, 
    balance = "g1",
    rec.rate = 1/100,
    rec.rate.g2 = 1/100,
    nsteps = 1000, 
    nsims = 1,
    verbose = F
  )
  expect_equal(max(x$s.num+x$i.num), 501)
  expect_equal(max(x$s.num+x$i.num+x$s.num.g2+x$i.num.g2), 1001)
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SIS model, 2 groups, no vital, 5 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = "SIS", 
    groups = 2, 
    s.num = 500, 
    i.num = 1, 
    s.num.g2 = 500, 
    i.num.g2 = 0, 
    trans.rate = 0.2, 
    act.rate = 0.25, 
    trans.rate.g2 = 0.1, 
    balance = "g1",
    rec.rate = 1/100,
    rec.rate.g2 = 1/100,
    nsteps = 1000, 
    nsims = 5,
    verbose = F
  )
  expect_equal(max(x$s.num+x$i.num), 501)
  expect_equal(max(x$s.num+x$i.num+x$s.num.g2+x$i.num.g2), 1001)
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SIS model, 1 group, vital, 1 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = 'SIS', 
    groups = 1,
    s.num = 500, 
    i.num = 1,
    trans.rate = 0.2, 
    act.rate = 0.5,
    rec.rate = 1/50,
    b.rate = 1/100,
    ds.rate = 1/100,
    di.rate = 1/90,
    nsteps = 500, 
    nsims = 1, 
    verbose = F
  )
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SIS model, 1 group, vital, 5 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = 'SIS', 
    groups = 1,
    s.num = 500, 
    i.num = 1,
    trans.rate = 0.2, 
    act.rate = 0.5,
    rec.rate = 1/50,
    b.rate = 1/100,
    ds.rate = 1/100,
    di.rate = 1/90,
    nsteps = 500, 
    nsims = 5, 
    verbose = F
  )
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SIS model, 2 group, vital, 1 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = 'SIS', 
    groups = 2,
    s.num = 500, 
    i.num = 1,
    s.num.g2 = 500,
    i.num.g2 = 1,
    trans.rate = 0.2, 
    trans.rate.g2 = 0.1,
    act.rate = 0.5,
    balance = 'g1',
    rec.rate = 1/50,
    rec.rate.g2 = 1/50,
    b.rate = 1/100,
    b.rate.g2 = NA,
    ds.rate = 1/100,
    ds.rate.g2 = 1/100,
    di.rate = 1/90,
    di.rate.g2 = 1/90,
    nsteps = 500, 
    nsims = 1, 
    verbose = F
  )
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

test <- 'SIS model, 2 group, vital, 5 sim'
EpiModel:::mcat(test)
test_that(test, {
  if (seed == TRUE) set.seed(98102)
  x <- epiICM(
    type = 'SIS', 
    groups = 2,
    s.num = 500, 
    i.num = 1,
    s.num.g2 = 500,
    i.num.g2 = 1,
    trans.rate = 0.2, 
    trans.rate.g2 = 0.1,
    act.rate = 0.5,
    balance = 'g1',
    rec.rate = 1/50,
    rec.rate.g2 = 1/50,
    b.rate = 1/100,
    b.rate.g2 = NA,
    ds.rate = 1/100,
    ds.rate.g2 = 1/100,
    di.rate = 1/90,
    di.rate.g2 = 1/90,
    nsteps = 500, 
    nsims = 5, 
    verbose = F
  )
  capture.output(print(x),
                 summary(x, time=150),
                 file='NULL')
  if (plots == TRUE) {
    plot(x, sub=test)
    plot(x, sub=test, popfrac=F)
  }
  cat('passed')
  EpiModel:::test.epiICM(x)
})

invisible(file.remove('NULL'))
if (plots == TRUE) invisible(dev.off())
}
