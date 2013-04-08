\name{expanddf}
\alias{expanddf}
\title{Expand a Frequency-Weighted Data Frame}
\usage{
  expanddf(df, count)
}
\arguments{
  \item{df}{data frame with one row per combination of
  variables}

  \item{count}{variable on the data frame with the count
  variable}
}
\description{
  This function expands a frequency-weighted data frame by
  the specified weights.
}
\examples{
df <- expand.grid(a=0:1, b=0:1, c=0:1)
df$count <- rpois(nrow(df), 5)
( df2 <- expanddf(df, df$count) )
}
\author{
  Samuel M. Jenness <sjenness@uw.edu>
}

