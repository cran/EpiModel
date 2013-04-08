\name{mysummary}
\alias{mysummary}
\title{Summarizing Regression Model Fits}
\usage{
  mysummary(fit, exp = FALSE, robust = FALSE, digits = 3)
}
\arguments{
  \item{fit}{fitted regression model}

  \item{exp}{exponentiate the summary coefficients}

  \item{robust}{robust standard errors}

  \item{digits}{digits to provide in output}
}
\description{
  Provides summary statistics from a regression model fit,
  including exponentiated coefficients, robust standard
  errors, and confidence intervals.
}
\examples{
y <- rbinom(100, 1, 0.5)
x <- ifelse(y == 1, rnorm(length(y[y==1]), 100, 20),
            rnorm(length(y[y==0]), 50, 20))
mod <- glm(y ~ x, family=binomial)
summary(mod)
mysummary(mod, exp=TRUE, robust=TRUE)
}
\author{
  Samuel M. Jenness <sjenness@uw.edu>
}

