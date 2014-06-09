#' 
#' Mathematical Modeling of Infectious Disease
#' 
#' \tabular{ll}{
#'    Package: \tab EpiModel\cr
#'    Type: \tab Package\cr
#'    Version: \tab 1.0\cr
#'    Date: \tab 2014-06-08\cr
#'    License: \tab GPL (>= 2)\cr
#'    LazyLoad: \tab yes\cr
#' }
#' 
#' @details 
#' The \code{EpiModel} package provides functions for building, solving, and 
#' plotting mathematical models of infectious disease. The goals of the 
#' package are to provide basic tools for modeling in multiple
#' frameworks for pedagogical purposes, and to support users to 
#' develop and expand these tools using the package's utility functions 
#' for their own research.
#' 
#' \code{EpiModel} currently provides functionality for three classes of epidemic
#' models:
#' \itemize{
#'  \item \strong{Deterministic Compartmental Models:} these continuous-time
#'        models are solved using ordinary differential equations. \code{EpiModel} 
#'        allows for easy specification of sensitivity models to compare multiple 
#'        runs of the same model with different parameter values.
#'  \item \strong{Stochastic Individual Contact Models:} a novel class of microsimulation
#'        models were developed to mirror the deterministic models but add random
#'        variation in all components of the transmission dynamics system,
#'        from infection to recovery to vital dynamics (births and deaths).
#'  \item \strong{Stochastic Network Models:} using the underlying statistical framework
#'        of dynamic exponential random graph models (ERGMs) recently developed in
#'        the \code{Statnet} suite of software in R, our network models  
#'        simulate partnership formation and dissolution stochastically according to 
#'        a user-defined statistical model, as well as disease spread on that network.
#'  }
#'  Future additions to the package will expand the varieties of models within
#'  each of these classes that may be run "out-of-the-box," as well as the
#'  helper utility functions that support users' own expansion of network 
#'  models specifically.
#'  
#' \code{EpiModel} supports three infectious disease types to be run across all of the
#'  three classes:
#'  \itemize{
#'   \item \strong{Susceptible-Infectious (SI):} a two-state disease in which there
#'         is life-long infection without recovery. HIV/AIDS is one example, 
#'         although for this case it is more common to model infection stages as 
#'         separate compartments (forthcoming in a future release).
#'   \item \strong{Susceptible-Infectious-Recovered (SIR):} a three-stage disease in 
#'         which one has life-long recovery with immunity after infection. Measles 
#'         is one example, but modern models for the disease also require 
#'         consideration of vaccination patterns in the population (forthcoming 
#'         in a future release).
#'   \item \strong{Susceptible-Infectious-Susceptible (SIS):} a two-stage disease 
#'         in which one may transition back and forth from the susceptible to 
#'         infected states throughout life. Examples include bacterial sexually 
#'         transmitted diseases like gonorrhea.
#'  }
#'
#' The core functions for parameterizing and solving the three model classes 
#' are:
#' \itemize{
#'  \item \code{\link{dcm}} for deterministic compartmental epidemic models.
#'  \item \code{\link{icm}} for individual contact epidemic models.
#'  \item \code{\link{netest}} for estimating the statistical models underlying 
#'        partnership formation and dissolution used in stochastic network epidemic 
#'        models. This function is a wrapper around the \code{ergm} and \code{stergm} 
#'        functions in the \code{ergm} and \code{tergm} packages, respectively, with 
#'        additional diagnostic tables and plots useful for epidemic modeling.
#'  \item \code{\link{netsim}} for stochastic network epidemic models, 
#'        with a given network model fit or set of network simulations from 
#'        \code{\link{netest}}. For models involving dependence of disease 
#'        trajectories on the network structure (e.g., disease causes death, 
#'        which dissolves partnerships), it is not necessary to pre-simulate the 
#'        networks since each disease simulation re-simulates the network at each 
#'        time step. A help page providing an overview of the internals of 
#'        \code{\link{netsim}} that may be useful for adapting and expanding the 
#'        software for novel research is available at \code{\link{modules.net}}.
#' }
#'
#' @references \url{http://www.statnet.org/trac/wiki/EpiModel}
#'
#' @name EpiModel-package
#' @aliases EpiModel
#' @import deSolve ergm network networkDynamic RColorBrewer statnet.common 
#'         tergm 
#' @docType package
#' @keywords package
#' 
NULL
