\name{transco}
\alias{transco}
\title{Obtain Transparent Colors}
\usage{
  transco(hexcol, alpha = 1, invis = TRUE)
}
\arguments{
  \item{hexcol}{hexidecimal color}

  \item{alpha}{transparency level (where 0 is transparent
  and 1 is opaque)}

  \item{invis}{supresses printing of the RGB color}
}
\description{
  This function returns a RGB transparent color from a
  defined hexidecimal color.
}
\examples{
n <- 10
x <- sort(sample(1:200, n))
y <- 10 + 2*x + rnorm(n, 0, 10)
z <- rpois(n, 10)

tcol <- transco('steelblue', 0.5)

plot(x,y, cex=z/4, pch=21, col='black', bg=tcol, lwd=1.2, axes=FALSE,
     ylim=c(0,500), xlim=c(0,250), yaxs='r', xaxs='r')
axis(2, seq(0,500,100), col='white', las=2, cex.axis=0.9, mgp=c(2,0.5,0))
axis(1, seq(0,250,50), cex.axis=0.9, mgp=c(2,0.5,0))
abline(h=seq(100,500,100), col=transco('black', 0.35))
}
\author{
  Samuel M. Jenness <sjenness@uw.edu>
}
\keyword{color}

