% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/net.accessor.R
\name{unique_id-tools}
\alias{unique_id-tools}
\alias{get_unique_ids}
\alias{get_posit_ids}
\title{Convert Unique identifiers from/to Positional identifiers}
\usage{
get_unique_ids(dat, posit_ids = NULL)

get_posit_ids(dat, unique_ids = NULL)
}
\arguments{
\item{dat}{a Master list object of network models}

\item{posit_ids}{a vector of node positional identifiers (default = NULL)}

\item{unique_ids}{a vector of node unique identifiers (default = NULL)}
}
\value{
a vector of unique or positional identifiers
}
\description{
EpiModel refers to its nodes either by positional identifiers
             (\code{posit_ids}), the position of the node in the \code{attr}
             vectors or by unique identifiers (\code{unique_ids}), allowing
             to refer to nodes even after they are deactivated
}
\section{All elements}{

  When  \code{unique_ids} or \code{get_posit_ids} is NULL (default)
  the full list of positional IDs or unique IDs is returned
}

\section{Deactivated nodes}{

  When providing \code{unique_ids} of deactivated nodes to
  \code{get_posit_ids}, \code{NA}s are returned instead and a warning is
  produced.
}

