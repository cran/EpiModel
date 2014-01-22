\name{transco}
\alias{transco}
\title{Obtain Transparent Colors}
\usage{
transco(col, alpha = 1, invisible = TRUE)
}
\arguments{
  \item{col}{vector of any of the three kinds of \code{R}
  color specifications (named, hexadecimal, or positive
  integer; see \code{\link{col2rgb}})}

  \item{alpha}{transparency level, where 0 is transparent
  and 1 is opaque}

  \item{invisible}{supresses printing of the RGB color}
}
\value{
A vector of length equal to the input \code{col} vector,
containing the transformed color values in hexidemical
format.
}
\description{
This function returns an RGB transparent color from an
input R color
}
\details{


The purpose of this function is to facilitate color
transparency, which is used widely in \code{EpiModel}
plots. This is an internal function that is not ordinarily
called by the end-user.
}
\examples{
# Three-variable bubble plot to show transparency
n <- 25
x <- sort(sample(1:200, n))
y <- 10 + 2*x + rnorm(n, 0, 10)
z <- rpois(n, 10)

# Input named R color and transparency level
tcol <- transco("steelblue", 0.5)

plot(x,y, cex=z/4, pch=21, col="black", bg=tcol, lwd=1.2, axes=FALSE,
     ylim=c(0,500), xlim=c(0,250), yaxs="r", xaxs="r")
axis(2, seq(0,500,100), col="white", las=2, cex.axis=0.9, mgp=c(2,0.5,0))
axis(1, seq(0,250,50), cex.axis=0.9, mgp=c(2,0.5,0))
abline(h=seq(100,500,100), col=transco("black", 0.35))
}
\seealso{
\code{\link{rgb}}, \code{\link{col2rgb}}
}
\keyword{colorUtils}
\keyword{internal}

