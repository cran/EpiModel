
#' @title Modules for Stochastic Network Models
#' 
#' @description Stochastic network models of infectious disease epidemics 
#' simulated through the \code{epiNet} class models require a complex series of
#' statistical modeling, random simulations, and bookkeeping. After the network
#' model is estimated, disease is simulated upon networks using the 
#' \code{\link{epiNet.simTrans}} function. Within this function are a series of
#' steps to initialize the simulation, and then randomly set new infections, 
#' recoveries, and vital dynamics on the network. 
#' 
#' Addressing novel research questions will likely require modifying or expanding
#' the existing functionality of \code{epiNet.simTrans} code. Our software design
#' goal has been to "modularize" the internal processes within 
#' \code{epiNet.simTrans} in order to present the procedures of these network
#' models clearly for easy editing and supplementation. 
#' 
#' A set of related procedures are known as modules. In this help page, we present
#' a brief overview of the module functions in the order in which they are used
#' within \code{epiNet.simTrans} to help guide users on the "big picture" for 
#' their own code expansion. Note that these functions are not shown on the main
#' help page index since they are not called directly by the end-user. Review the
#' help pages for these functions for further details.
#' 
#' @section Initialization Functions:
#' These functions set up the nodal attributes like disease status on the network
#' at the starting time step of disease simualtion, \eqn{t_1}. For multiple-simulation
#' function calls, these are reset at the beginning of each individual simulation.
#' \itemize{
#'  \item \code{\link{init.status}}: sets which nodes are initially infected,
#'    either through the \code{i.num} or \code{i.ids} parameters.
#'  \item \code{\link{init.inf.time}}: sets the time of infection for those
#'    nodes which were infected in \code{init.status}. 
#'  \item \code{\link{init.pids}}: establishes persistent ID numbers for bipartite
#'    network simulations with vital dynamics.
#' }
#'  
#' @section Infection and Recovery Functions:
#' The main disease simulation occurs at each time step given the current state
#' of the network at that step. Infection of nodes is simulated as a function of 
#' attributes of the nodes and the dyad (e.g., the mode number or infection time
#' of the infected node in a partnership). Recovery of nodes is likewise simulated
#' as a function of nodal attributes of those infected nodes. These functions 
#' also conduct analysis on the network to save summary disease measures such as 
#' disease incidence.
#' \itemize{
#'  \item \code{\link{infection}}: simulates disease transmission given an edgelist
#'    of serodisdant partnerships by calculating the relevant transmission and act
#'    rates for each nodal pair, and then updated the nodal attributes and summary
#'    statistics.
#'  \item \code{\link{discord.edgelist}}: determines from the network object which
#'    edges are active, and then the subset of active edges which are comprised of
#'    a serodiscordant (one node susceptible and one node infected) nodal pairs.
#'  \item \code{\link{recovery}}: simulates recovery from infection either to a
#'    lifelong immune state (for SIR models) or back to the susceptible state 
#'    (for SIS models), as a function of the recovery rate specified in the
#'    \code{rec.rate} parameter.
#' }
#'  
#'
#' @section Vital Dynamics Functions:
#' Vital dynamics such as birth and death processes are simulated at each time
#' step to update entries into and exits from the network. These are used in 
#' dependent network models.
#' \itemize{
#'  \item \code{\link{deaths.sus}}: randomly simulates death for susceptible status
#'    nodes given the death rate specified in the \code{ds.rate} parameter. This
#'    involves deactivating susceptible nodes.
#'  \item \code{\link{deaths.inf}}: randomly simulates death for infected status
#'    nodes given the death rate specified in the \code{di.rate} parameter. This
#'    involves deactivating infected nodes.
#'  \item \code{\link{deaths.rec}}: randomly simulates death for recovered status
#'    nodes given the death rate specified in the \code{dr.rate} parameter. This
#'    involves deactivating recovered nodes.
#'  \item \code{\link{births}}: randomly simulates new births into the network
#'    given the current population size and the birth rate specified in the 
#'    \code{b.rate} parameter. This involves adding new nodes into the network.
#' }
#'
#' @name epiNetModules
#' @aliases epiNetModules
#' 
NULL