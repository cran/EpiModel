\name{dropdf}
\alias{dropdf}
\title{Drop Variables from a Data Frame}
\usage{
  dropdf(df, vars)
}
\arguments{
  \item{df}{data frame with multiple variables}

  \item{vars}{vector of variables in character format}
}
\description{
  This function drops named variables from a data frame
}
\examples{
w <- x <- y <- z <- 1:10
df <- data.frame(w,x,y,z)
( newdf <- dropdf(df, c('x', 'y')) )
}
\author{
  Samuel M. Jenness <sjenness@uw.edu>
}

