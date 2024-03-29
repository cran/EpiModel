
#' @title Definition for absdiffnodemix ERGM Term
#'
#' @description This function defines and initializes the absdiffnodemix ERGM
#'              term that allows for targeting homophily based on a non-binary
#'              attribute (e.g., age) by combinations of a binary attribute
#'              (e.g., race).
#'
#' @param nw An object of class \code{network}.
#' @param arglist A list of arguments as specified in the \code{ergm.userterms}
#'        package framework.
#' @param ... Additional data passed into the function as specified in the
#'        \code{ergm.userterms} package framework.
#'
#' @details
#' This ERGM user term was written to allow for age-based homophily in
#' partnership formation that is heterogeneous by race. The \code{absdiff}
#' component targets the distribution of age mixing on that continuous
#' variable, and the \code{nodemix} component differentiates this for
#' black-black, black-white, and white-white couples.
#'
#' @aliases absdiffnodemix
#'
InitErgmTerm.absdiffnodemix <- function(nw, arglist, ...) {

  a <- check.ErgmTerm(nw,
                      arglist,
                      directed = FALSE,
                      bipartite = FALSE,
                      varnames = c("attr", "by"),
                      vartypes = c(ERGM_VATTR_SPEC, ERGM_VATTR_SPEC),
                      defaultvalues = list(NULL, NULL),
                      required = c(TRUE, TRUE))

  nodecov <- ergm_get_vattr(a$attr, nw, accept = "numeric")
  nodecovby <- ergm_get_vattr(a$by, nw)
  nodecovbyname <- attr(nodecovby, "name")
  u <- sort(unique(nodecovby))
  if (any(is.na(nodecovby))) {
    u <- c(u, NA)
  }

  nodecovby <- match(nodecovby, u, nomatch = length(u) + 1)
  ui <- seq(along = u)

  uui <- matrix(seq_along(ui) ^ 2, length(ui), length(ui))
  urm <- t(sapply(ui, rep, length(ui)))
  ucm <- sapply(ui, rep, length(ui))
  uun <- outer(u, u, paste, sep = ".")
  uui <- uui[upper.tri(uui, diag = TRUE)]
  urm <- urm[upper.tri(urm, diag = TRUE)]
  ucm <- ucm[upper.tri(ucm, diag = TRUE)]
  uun <- uun[upper.tri(uun, diag = TRUE)]

  inputs <-  c(length(nodecov), length(urm), nodecov, nodecovby, urm, ucm)

  list(name = "absdiffnodemix",
       coef.names = paste("absdiffnodemix", attr(nodecov, "name"),
                          nodecovbyname, uun, sep = "."),
       pkgname = "EpiModel",
       inputs = inputs,
       dependence = FALSE)
}


#' @title Definition for absdiffby ERGM Term
#'
#' @description This function defines and initializes the absdiffby ERGM term
#'              that allows for representing homophily with respect to a
#'              non-binary attribute (e.g., age) differentially by a binary
#'              attribute (e.g., sex).
#'
#' @param nw An object of class \code{network}.
#' @param arglist A list of arguments as specified in the \code{ergm.userterms}
#'        package framework.
#' @param ... Additional data passed into the function as specified in the
#'        \code{ergm.userterms} package framework.
#'
#' @details
#' This ERGM user term was written to allow for age-based homophily in
#' partnership formation that is asymmetric by sex. The \code{absdiff} component
#' targets age-based homophily while the \code{by} component allows that to be
#' structured by a binary attribute such as "male", in order to enforce an
#' offset in the average difference. This allows, for example, a average age
#' difference in partnerships, but with males (on average) older than females.
#'
#' @aliases absdiffby
#'
InitErgmTerm.absdiffby <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw,
                      arglist,
                      directed = FALSE,
                      bipartite = FALSE,
                      varnames = c("attr", "by", "assym"),
                      vartypes = c(ERGM_VATTR_SPEC, ERGM_VATTR_SPEC, "numeric"),
                      required = c(TRUE, TRUE, TRUE),
                      defaultvalues = list(NULL, NULL, NULL))

  nodecov <- ergm_get_vattr(a$attr, nw, accept = "numeric")
  nodeby <- ergm_get_vattr(a$by, nw)
  coef.names <- paste("absdiffby", attr(nodecov, "name"),
                      attr(nodeby, "name"), sep = ".")

  list(name = "absdiffby",
       coef.names = coef.names,
       pkgname = "EpiModel",
       inputs = c(a$assym, nodecov, nodeby),
       dependence = FALSE,
       emptynwstats = 0
  )
}


#' @title Definition for fuzzynodematch ERGM Term
#'
#' @description This function defines and initializes the fuzzynodematch ERGM
#'              term that allows for generalized homophily.
#'
#' @param nw An object of class \code{network}.
#' @param arglist A list of arguments as specified in the \code{ergm.userterms}
#'        package framework.
#' @param ... Additional data passed into the function as specified in the
#'        \code{ergm.userterms} package framework.
#'
#' @details
#' This ERGM user term was written to allow for generalized homophily.The
#' \code{attr} term argument should specify a character vertex attribute
#' encoding the "venues" associated to each node.  The \code{split} argument
#' should specify a string that separates different "venues" in the attribute
#' value for each node, as handled by \code{strsplit} with \code{fixed = TRUE}.
#' For example, if \code{split} is \code{"|"} (the default), and the attribute
#' value for a given node is \code{"a12|b476"}, then the associated venues for
#' this node are \code{"a12"} and \code{"b476"}.  The empty string \code{""} is
#' interpreted as "no venues".
#'
#' If the \code{binary} term argument is \code{FALSE} (the default), the change
#' statistic for an on-toggle is the number of unique venues associated to both
#' nodes (informally speaking, this could be described as the number of venues
#' on which the two nodes "match"); if \code{binary} is \code{TRUE}, the change
#' statistic for an on-toggle is \code{1} if any venue is associated to both
#' nodes, and \code{0} otherwise.
#'
#' @aliases fuzzynodematch
#'
InitErgmTerm.fuzzynodematch <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("attr", "split", "binary"),
                      vartypes = c(ERGM_VATTR_SPEC, "character", "logical"),
                      defaultvalues = list(NULL, "|", FALSE),
                      required = c(TRUE, FALSE, FALSE))

  nodecov <- ergm_get_vattr(a$attr, nw, accept = "character")
  venues <- strsplit(nodecov, split = a$split, fixed = TRUE)

  ## drop "" from venues and enforce uniqueness of venues for each node
  venues <- lapply(venues, function(x) unique(x[nchar(x) > 0L]))

  ## record number of venues and offset in position for each node
  lengths <- unlist(lapply(venues, length))
  positions <- cumsum(lengths) - lengths

  ## convert venues to vector
  venues <- unlist(venues)

  ## convert venues from strings to integers
  levels <- sort(unique(venues))
  venues <- match(venues, levels)

  ## sort venues for each node
  venues <- unlist(lapply(seq_len(network.size(nw)), function(i) sort(venues[positions[i] + seq_len(lengths[i])])))

  binary <- a$binary

  list(name = "fuzzynodematch",
       coef.names = paste("fuzzynodematch", attr(nodecov, "name"),
                          binary, sep = "."),
       binary = as.integer(binary),
       venues = as.integer(venues),
       lengths = as.integer(lengths),
       positions = as.integer(positions),
       dependence = FALSE,
       minval = 0)
}
