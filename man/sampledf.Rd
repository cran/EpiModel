\name{sampledf}
\alias{sampledf}
\title{Sample Rows of a Data Frame}
\usage{
  sampledf(df, size, replace = FALSE, sortby = "")
}
\arguments{
  \item{df}{data frame to be sampled}

  \item{size}{number of rows to sample}

  \item{replace}{sample with replacement}

  \item{sortby}{sort the data frame by a variable in the
  data frame}
}
\description{
  This function samples rows from a data frame.
}
\examples{
df <- expand.grid(a=0:1, b=0:1, c=0:1)
df$count <- rpois(nrow(df), 5)
df2 <- expanddf(df, df$count)
( df3 <- sampledf(df2, size=5, replace=FALSE, sortby='c') )
}
\author{
  Samuel M. Jenness <sjenness@uw.edu>
}

