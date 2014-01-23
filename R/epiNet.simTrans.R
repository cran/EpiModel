#' 
#' @title Stochastic Network Epidemic Models 
#'
#' @description This function simulates stochastic network-based epidemic models 
#'  for infectious disease.
#' 
#' @param x an \code{EpiModel} object of class \code{epiNet.est} for dependent
#'   disease-network simulations or class \code{epiNet.simNet} for independent
#'   disease-network simulations (see details).
#' @param type disease type to be modeled, with the choice of \code{"SI"} for
#'   Susceptible-Infected diseases, \code{"SIR"} for Susceptible-Infected-Recovered
#'   diseases, and \code{"SIS"} for Susceptible-Infected-Susceptible diseases.
#' @param vital if \code{TRUE}, incorporate vital dynamics (births and deaths) into
#'   the model. This requires passing a birth rate and death rate for each 
#'   compartment in the model (see details).
#' @param i.num number in the population initially infected. In bipartite models,
#'   infections are randomly allocated across modes with equal probability.
#' @param i.ids a vector of node ID numbers to set those specific nodes as 
#'   initially infected at baseline. One must use either \code{i.num} or 
#'   \code{i.ids} to set initial infected. Setting infected IDs here overrides 
#'   \code{i.rand=TRUE}.
#' @param i.rand if \code{TRUE} and using \code{i.num}, sets infection based on 
#'   random binomial draws from distribution implied by i.num. 
#' @param trans.rate probability of transmission given an act or contact between
#'   a susceptible and an infected person in the population. In bipartite models
#'   this is the rate governing the infection to the mode 1 members.
#' @param trans.rate.m2 probability of transmission given an act or contact 
#'   between a susceptible mode 2 person and an infected mode 1 person in the 
#'   population. It is the rate governing the infection to the mode 2 members.
#' @param act.rate average number of acts governing transmission \emph{per 
#'   partnership} per unit time, regardless of disease status (see details).
#' @param rec.rate average rate of recovery with immunity (in \code{SIR} models) 
#'   or re-susceptibility (in \code{SIS} models). The recovery rate is a function 
#'   of the disease duration as 1/duration of disease. The recovery rate is a function of 
#'   the disease duration as 1/duration of disease.
#' @param b.rate birth rate into the network. For one-mode models, the birth rate is
#'   the probability of a new birth per person per unit time. For bipartite models, 
#'   the birth rate is the rate per "mode 1 persons" per unit time.
#' @param ds.rate mortality rate for susceptibles, i.e. the probability of death 
#'   per susceptible person per unit time.
#' @param di.rate mortality rate for infecteds, i.e., the probability of death per 
#'   infected person per unit time.
#' @param dr.rate mortality rate for recovered, i.e., the probability of death per 
#'   recovered person per unit time.
#' @param nsteps number of time steps to simulate the model over. This must be a
#'   positive integer. In independent models (\code{vital=FALSE}), \code{nsteps} 
#'   is set by default through the \code{\link{epiNet.simNet}} network 
#'   simulation, although \code{nsteps} may be set here to a smaller number of 
#'   steps. For dependent models (\code{vital=FALSE}), \code{nsteps} must be set 
#'   since there is no default.
#' @param sims.per.nw In independent models (\code{vital=FALSE}), the number of 
#'   disease simulations per network simulated in \code{\link{epiNet.simNet}}. 
#'   In dependent models, (\code{vital=TRUE}), the total number of disease 
#'   simulations, since the network is resimulated at each time step.
#' @param tea if \code{TRUE}, use temporally extended attributes to store disease 
#'   status information. This tends to slow the simulations but allows for more
#'   flexible network plotting and animations.
#' @param save.statmat if \code{TRUE} and \code{tea=FALSE}, save a disease 
#'   status matrix to the main output object. This is a time by n sized matrix
#'   containing the disease status of all nodes at each time; it may be used 
#'   instead of temporally extended attributes to store all disease history for
#'   later plotting and analysis. 
#' @param save.stats if \code{TRUE} and \code{vital=TRUE}, save network statistics 
#'   in a \code{data.frame} out to the main output object. The specific statistics
#'   to be requested are set in the \code{stats.formula} parameter. Note that in
#'   independent models, \code{save.stats=TRUE} is overridden because the networks
#'   are their statistics of interest have already been generated during the
#'   \code{\link{epiNet.simNet}} simulation.
#' @param stats.formula formula to specify diagnostic statistics, where the 
#'   default is formation formula (specified during the estimation phase).
#' @param save.trans if \code{TRUE}, save a transmission \code{data.frame} for 
#'   each simulation out to the main object. This object contains one row for
#'   each transmission event that has occurred (see \code{\link{discord.edgelist}}).
#' @param save.network if \code{TRUE}, save a \code{networkDynamic} object 
#'   containing full edge and nodal history for each simulation out to the 
#'   main object.
#' @param verbose if \code{TRUE}, print simulation progress to the console.
#' @param plot.prog if \code{TRUE}, dynamically plot simulation progress to the 
#'   plot window.
#' @param ... additional arguments to pass to \code{simulate.stergm}.
#' 
#' @details
#' This is the main disease simulation function for stochastic network-based 
#' epidemic models in \code{EpiModel}. This class of models differs from both 
#' deterministic compartmental models in \code{epiDCM} and their stochastic 
#' counterparts in \code{epiICM} by explicitly tracking partnership dyads and 
#' changes within those dyads over time. All three disease types (SI, SIR, and 
#' SIS) may be modeled in \code{epiNet.simTrans}.
#' 
#' If one is interested in extended these network models for novel research, it
#' will be necessary to understand what goes on "under-the-hood" of this simulation
#' function. To start, it is recommended that you review the 
#' \code{\link{epiNetModules}} help page, which gives an overview of all the 
#' functions that update the disease simulation at each time step. Since these
#' functions are not called directly by the end-user, they are not shown on the 
#' main help index.
#' 
#' @section Dependent versus Independent Models:
#' One important modelling choice is between a dependent and independent disease 
#' model. In an independent model, the disease simulation does not change the 
#' network structure, but in a dependent model demographics and disease influence 
#' the network structure. See Section 4 of \href{../doc/Tutorial.pdf}{EpiModel Tutorial} 
#' for further clarification. 
#' 
#' In this current version of \code{EpiModel}, dependent models is synonymous 
#' with incorporating vital dynamics (births and deaths), however that is only 
#' one type of dependence, and other versions will be incorporated in future 
#' \code{EpiModel} releases. One incorporates vital dynamics into an 
#' \code{epiNet.simTrans} simulation by setting \code{vital=TRUE}. In these dependent 
#' simulations, the \code{x} parameter is the network estimation object from 
#' \code{\link{epiNet.est}}, whereas in independent simulations the \code{x} 
#' parameter is the network simulation object from \code{\link{epiNet.simNet}}. With 
#' the latter, the dynamic network is simulated in advance and outside of this 
#' disease simulation; with the former, the network is resimulated at each 
#' timestep, and only the model fit is needed to get started. 
#' 
#' @section The act.rate Parameter:
#' A key difference between these network models and the DCM/ICM models is their 
#' treatment of transmission events. In DCM and ICM, contacts or partnerships are 
#' mathematically instantaneous events: they have no duration in time, and thus
#' no changes may occur within them over time. In contrast, network models allow
#' for partnership durations defined by the dynamic network model, summarized in
#' the model dissolution coefficients calculated in \code{\link{dissolution.coefs}}.
#' Therefore, the \code{act.rate} parameter has a different interpretation here, 
#' where it is the number of acts defining transmission per partnership per unit 
#' time. 
#' 
#' @return
#' This function returns a list object of class \code{epiNet.simTrans} with the 
#' following elements:
#' \itemize{
#'  \item \strong{type:} disease type as specified in model parameter.
#'  \item \strong{modes:} network modes as set in the starting network for 
#'    estimation in the \code{\link{epiNet.est}} call.
#'  \item \strong{time:} a vector of time steps over which the model was solved.
#'  \item \strong{nruns:} number of independent model runs.
#'  \item \strong{*.num:} a vector (if \code{nruns=1}) or data frame (if 
#'    \code{nruns>1}) of compartment or state sizes over time, for each model run, 
#'    where * may equal \code{s}, \code{i}, or \code{r} for susceptible, infected, 
#'    or recovered, respectively. Bipartite models have equivalent elements with a
#'    \code{.m2} suffix.
#'  \item \strong{*.flow:} a vector (if \code{nruns=1}) or data frame (if 
#'    \code{nruns>1}) of transition sizes between compartments, for each model run, 
#'    where * may equal \code{si} for susceptible to infected flows, \code{is} for 
#'    infected to susceptible flows, \code{ir} for infected to recovered flows, 
#'    \code{b} for birth in-flows, \code{ds} for susceptible death out-flows, 
#'    \code{di} for infected death out-flows, \code{dr} for recovered death out-flows.
#'    Bipartite models have equivalent elements with a \code{.m2} suffix.
#'  \item \strong{stat.mat:} a list of matrices, one for each simulation, containing
#'    the disease status history, specified to be saved with the \code{save.statmat}
#'    parameter.
#'  \item \strong{trans:} list of data frames, one for each simulation, containing
#'    the transmission history, specified to be saved with the \code{save.trans}
#'    parameter.
#'  \item \strong{network:} a list of \code{networkDynamic} objects, one for each 
#'    simulation, containing the complete edge and nodal histories, specified to 
#'    be saved with the \code{save.network} parameter.
#'  \item \strong{vital:} logical, whether vital dynamics were specified in the 
#'      parameterization.
#'  \item \strong{nsims:} number of disease simulations conducted.
#'  \item \strong{call:} exact model call.
#' }
#' 
#' @references
#' Goodreau SM, Carnegie NB, et al. What drives the US and Peruvian HIV epidemics 
#' in men who have sex with men (MSM)? PloS One. 2012; 7(11): e50522.
#' 
#' @seealso \code{\link{epiNet.est}} for network estimation, 
#'  \code{\link{epiNet.simNet}} for network simulation, 
#'  \code{\link{plot.epiNet.simTrans}} for plots, and 
#'  \code{\link{summary.epiNet.simTrans}} for data summaries.
#' 
#' @keywords model
#' @export
#' 
#' @examples
#' \dontrun{
#' ## See the EpiModel Tutorial vignette for more details and examples ##
#' 
#' ## Independent SI Model ##
#' # Initialize network and set network model parameters
#' nw <- network.initialize(n=500, bipartite=250, directed=FALSE)
#' formation <- ~ edges + b1degree(0:1) + b2degree(0:1)
#' target.stats <- c(165, 100, 137.5, 120, 102.5)
#' dissolution <- ~ offset(edges)
#' coef.diss <- dissolution.coefs(dissolution, duration=25)
#' 
#' # Estimate the ERGM models (see help for epiNet.est)
#' # Skipping model diagnostics for this, but should always run
#' est1 <- epiNet.est(
#'   nw, 
#'   formation, 
#'   dissolution, 
#'   target.stats, 
#'   coef.diss)
#' 
#' # Simulate 10 networks from model fit (see help for epiNet.simNet)
#' # Skipping model diagnostics for this, but should always run
#' nwsims <- epiNet.simNet(est1, nsteps = 250, nsims = 10)
#' 
#' # Independent epidemic model simulation for SI disease
#' # Transmissibility is higher to mode 1 than mode 2
#' sim.SI <- epiNet.simTrans(
#'   nwsims, 
#'   type = "SI",
#'   vital = FALSE,
#'   trans.rate = 0.3, 
#'   trans.rate.m2 = 0.15, 
#'   i.num = 50,
#'   sims.per.nw = 1)
#' 
#' # Print, plot, and summarize the results
#' sim.SI
#' plot(sim.SI, type="sim")
#' plot(sim.SI, y=c("i.num", "i.num.m2"), sim.lines=TRUE, 
#'      sim.col=c("steelblue", "firebrick"), qnts=0.5)
#' plot(sim.SI, type="network", at=50, col.inf=TRUE, shp.bip="square")
#' summary(sim.SI, time=50)
#' 
#' ## Dependent SIR Model ##
#' coef.diss <- dissolution.coefs(dissolution, duration=25, d.rate=0.0021)
#' est2 <- epiNet.est(
#'   nw, 
#'   formation, 
#'   dissolution, 
#'   target.stats, 
#'   coef.diss)
#' 
#' sim.SIR <- epiNet.simTrans(
#'   est2,
#'   type = "SIR",
#'   vital = TRUE,
#'   i.num = 50,
#'   trans.rate = 0.3,
#'   trans.rate.m2 = 0.15,
#'   rec.rate = 0.02,
#'   b.rate = 0.0021*2,
#'   ds.rate = 0.0021,
#'   di.rate = 0.0021*1.1,
#'   dr.rate = 0.0021,
#'   nsteps = 250,
#'   sims.per.nw = 10)
#' 
#' # Print, plot, and summarize the results
#' sim.SIR
#' plot(sim.SIR, type="sim")
#' plot(sim.SIR, type="network")
#' summary(sim.SIR, time=100)
#' }
#' 
epiNet.simTrans <- function(x, 
                            type,
                            vital = FALSE,
                            i.num,
                            i.ids,
                            i.rand = TRUE,
                            trans.rate,
                            trans.rate.m2,
                            act.rate = 1,
                            rec.rate,
                            b.rate,
                            ds.rate,
                            di.rate,
                            dr.rate,
                            nsteps, 
                            sims.per.nw = 1,
                            tea = FALSE,
                            save.statmat = TRUE,
                            save.stats = FALSE,
                            stats.formula,
                            save.trans = TRUE,
                            save.network = TRUE,
                            verbose = TRUE,
                            plot.prog = FALSE,
                            ...) {
  
  ## Warnings and Global Settings ##
  {
    if (missing(type)) stop("Supply model type")
    if (vital == FALSE & class(x) == "epiNet.est")
      stop("For vital=FALSE models, x must be a epiNet.simNet object")
    if (vital == TRUE & class(x) == "epiNet.simNet")
      stop("For vital=TRUE models, x must be a epiNet.est object")
    if (tea == TRUE) save.statmat <- FALSE
  }
  
  
  ## Set Simulation Numbers ##
  {
    if (vital == FALSE) {
      nws <- 1:x$nsims
      nw.to.sim <- rep(nws, each=sims.per.nw)
      tot.sims <- length(nw.to.sim)
      
      # Default nsteps to nD object, error if longer than nD
      if (missing(nsteps)) nsteps <- x$nsteps
      if (nsteps > x$nsteps) {
        stop("nsteps here must be <= nsteps from network simulation: ", x$nsteps)
      }
    }
    if (vital == TRUE) {
      tot.sims <- sims.per.nw
    }   
  }
  
  
  ## Graphical/Console Settings
  {
    ops <- list(mar=par()$mar, mfrow=par()$mfrow, mgp=par()$mgp)
    if (verbose == TRUE) {
      if (tot.sims == 1) {
        cat("===============================")
        cat("\nStarting 1 Disease Simulation")
        cat("\n===============================")
      } else {
        cat("===============================")
        cat("\nStarting", tot.sims, "Disease Simulations")
        cat("\n===============================")
      }
    }
  }
  
  
    ### SIMULATION LOOP ###
    for (s in 1:tot.sims) {
    
      ## Master List for Data ##
      all <- list()
      
      
      ## Load NW Model Data ##
      {
        if (vital == FALSE) {
          nw <- x[[nw.to.sim[s]]]
          nw <- activate.vertices(nw, onset=1, terminus=Inf)
        }
        if (vital == TRUE) {
          nw <- simulate(x$fit)
          if (class(x$fit) == "stergm") {
            nw <- network.collapse(nw, at = 1)
          }
          nw <- activate.vertices(nw, onset=1, terminus=Inf)
        } 
        all$nw <- nw
        formation <- all$formation <- x$formation
        dissolution <- all$dissolution <- x$dissolution
        coef.form <- all$coef.form <- x$coef.form
        coef.diss <- all$coef.diss <- x$coef.diss
        constraints <-  x$constraints
          if (is.null(constraints)) constraints <- ~ .
          all$constraints <- constraints
        modes <- all$modes <- ifelse(nw %n% "bipartite", 2, 1)
      }
      
      
      ### Warnings ###
      {
        if (vital == TRUE) {
          if (dissolution != ~ offset(edges))
            stop("Currently only ~offset(edges) dissolution models supported")
          if (missing(nsteps)) stop("Must supply nsteps if vital = TRUE")
          if (missing(b.rate)) 
            stop("Must supply b.rate if vital = TRUE")
          all$b.rate <- b.rate
          if (missing(ds.rate))
            stop("Must supply ds.rate if vital = TRUE")
          all$ds.rate <- ds.rate
          if (missing(di.rate))
            stop("Must supply di.rate if vital = TRUE")
          all$di.rate <- di.rate
          if (type == "SIR") {
            if (missing(dr.rate)) {
              stop("Must supply dr.rate if vital = TRUE and type = SIR")
            } else {
              all$dr.rate <- dr.rate
            }
          }
        } 
        if (modes == 1 && !(missing(trans.rate.m2)) && !is.null(trans.rate.m2)) {
          warning("Currently multiple trans.rates only with bipartite networks.
                   Using only trans.rate for this simulation.")
        }
        if (missing(trans.rate)) 
          stop("Specify a trans.rate")
        if (missing(trans.rate.m2)) 
          trans.rate.m2 <- NULL
        if (!is.null(trans.rate.m2) && (length(trans.rate.m2) != length(trans.rate))) 
          stop("trans.rate and trans.rate.m2 parameters should be the same length")
        ## Set/warnings for rec.rate
        if (type %in% c("SIR", "SIS")) {
          if (missing(rec.rate)) {
            stop("Must supply a rec.rate for SIR/SIS type models")
          } else {
            all$rec.rate <- rec.rate
          }
        }
      }
      
      
      ## Set Simulation Objects on all ##
      {
        all$type <- type
        all$nsteps <- nsteps
        all$tot.sims <- tot.sims
        all$vital <- vital
        all$act.rate <- act.rate
        all$trans.rate <- trans.rate
        all$trans.rate.m2 <- trans.rate.m2
        if (!missing(i.num)) all$i.num.t0 <- i.num
        if (!missing(i.ids)) all$i.ids.t0 <- i.ids
        all$i.rand <- i.rand
        all$save.statmat <- save.statmat
      }
    
      
      ## Infection Status and Time Modules ##
      all <- init.status(all, tea)
      all <- init.inf.time(all)
    
      
      ## Set Vital Dynamics ##
      if (vital == TRUE) {
        
        # If bipartite, activate persistant IDs
        if (modes == 2) all$nw <- init.pids(all$nw)
        
        # Output model stats
        if (missing(stats.formula)) 
          stats.formula <- formation
        
        # Activate edges
        suppressWarnings(
          all$nw <- simulate(all$nw,
                             formation = formation, 
                             dissolution = dissolution,
                             coef.form = coef.form, 
                             coef.diss = coef.diss$coef.crude,
                             constraints = constraints,
                             time.start = 1,
                             time.slices = 1,
                             time.offset = 0,
                             monitor = stats.formula, ...))
      }

      
      ### TIME LOOP ###  
      for (ts in 2:nsteps) {
      
        
        ### Infection Module ###
        all <- infection(all, at=ts)
        
        
        ### Recovery Module ###
        if (type %in% c("SIR", "SIS")) 
          all <- recovery(all, at=ts)
        
        
        ### Vital Dynamics Modules ###
        if (vital == TRUE) {
          all <- deaths.sus(all, ts)        # Susceptible deaths
          all <- deaths.inf(all, ts)        # Infected deaths
          all <- deaths.rec(all, ts)        # Recovered deaths
          all <- births(all, ts)            # Births
        }
        
        
        ### Resimulate network in dependent network
        if (vital == TRUE) {
          
          ## Population size check
          if ((modes == 1 && node.active(all$nw, at=ts, out="prev") == 0) ||
                (modes == 2 && node.active(all$nw, at=ts, out="prev", mode=1) == 0)) {
            zeropop <- TRUE 
          } else {
            zeropop <- FALSE
          }
          
          if (zeropop == FALSE) {
            suppressWarnings(
            all$nw <- simulate(all$nw,
                               formation = formation, 
                               dissolution = dissolution,
                               coef.form = coef.form, 
                               coef.diss = coef.diss$coef.adj,
                               constraints = constraints,
                               time.start = ts,
                               time.slices = 1,
                               time.offset = 0,
                               monitor = stats.formula, ...))
          }
        }
        
        
        ### Save Prevalence Vectors ###
        all <- get.prev(all, at = ts, set.all = TRUE)
        
        
        ### Popsize Edges Correction ###
        if (vital == TRUE)
          coef.form <- edges.correct(all, ts, coef.form)
        
        
        ### Progress Console and Plot ###
        if (verbose == TRUE) cons.prog(all, s, ts)
        if (plot.prog == TRUE) prog.plot(all, ts)
        
      } # end timestep loop
      
      # For sim1 create lists
      if (s == 1) {
        out <- list(type = type,
                    time = 1:nsteps,
                    modes = modes,
                    s.num = data.frame(all$s.num.m1),
                    i.num = data.frame(all$i.num.m1),
                    si.flow = data.frame(all$si.flow.m1))
        if (type == "SIR") {
          out$r.num <- data.frame(all$r.num.m1)
          out$ir.flow <- data.frame(all$recovs.m1)
        }
        if (type == "SIS") {
          out$is.flow <- data.frame(all$recovs.m1)
        }
        if (vital == TRUE) {
          out$b.flow <- data.frame(all$b.flow.m1)
          out$ds.flow <- data.frame(all$ds.flow.m1)
          out$di.flow <- data.frame(all$di.flow.m1)
          if (type == "SIR") {
            out$dr.flow <- data.frame(all$dr.flow.m1)
          }
        }
        if (modes == 2) {
          out$s.num.m2 = data.frame(all$s.num.m2)
          out$i.num.m2 = data.frame(all$i.num.m2)
          out$si.flow.m2 = data.frame(all$si.flow.m2)
          if (type == "SIR") {
            out$r.num.m2 <- data.frame(all$r.num.m2)
            out$ir.flow.m2 <- data.frame(all$recovs.m2)
          }
          if (type == "SIS") {
            out$is.flow.m2 <- data.frame(all$recovs.m2)
          }
          if (vital == TRUE) {
            out$b.flow.m2 <- data.frame(all$b.flow.m2)
            out$ds.flow.m2 <- data.frame(all$ds.flow.m2)
            out$di.flow.m2 <- data.frame(all$di.flow.m2)
            if (type == "SIR") {
              out$dr.flow.m2 <- data.frame(all$dr.flow.m2)
            }
          }
        }
        if (tea == FALSE && save.statmat == TRUE) 
          out$stat.mat <- list(all$stat.mat)
        if (vital == TRUE && save.stats == TRUE) {
          out$stats <- list(as.data.frame(attributes(all$nw)$stats))
        }
        if (save.trans == TRUE) {
          if (!is.null(all$transdf)) {
            row.names(all$transdf) <- 1:nrow(all$transdf)
            out$trans <- list(all$transdf)
          } else {
            out$trans <- list(data.frame())
          }
        }
        if (save.network == TRUE) {
          out$network <- list(all$nw)
        }
      } # end sim == 1 out
      # For sim2+, bind to lists
      if (s > 1) {
        out$s.num[,s] <- all$s.num.m1
        out$i.num[,s] <- all$i.num.m1
        out$si.flow[,s] <- all$si.flow.m1
        if (type == "SIR") {
          out$r.num[,s] <- all$r.num.m1
          out$ir.flow[,s] <- all$recovs.m1
        }
        if (type == "SIS") {
          out$is.flow[,s] <- all$recovs.m1
        }
        if (vital == TRUE) {
          out$b.flow[,s] <- all$b.flow.m1
          out$ds.flow[,s] <- all$ds.flow.m1
          out$di.flow[,s] <- all$di.flow.m1
          if (type == "SIR") {
            out$dr.flow[,s] <- all$dr.flow.m1
          }
        }
        if (modes == 2) {
          out$s.num.m2[,s] <- all$s.num.m2
          out$i.num.m2[,s] <- all$i.num.m2
          out$si.flow.m2[,s] <- all$si.flow.m2
          if (type == "SIR") {
            out$r.num.m2[,s] <- all$r.num.m2
            out$ir.flow.m2[,s] <- all$recovs.m2
          }
          if (type == "SIS") {
            out$is.flow.m2[,s] <- all$recovs.m2
          }
          if (vital == TRUE) {
            out$b.flow.m2[,s] <- all$b.flow.m2
            out$ds.flow.m2[,s] <- all$ds.flow.m2
            out$di.flow.m2[,s] <- all$di.flow.m2
            if (type == "SIR") {
              out$dr.flow.m2[,s] <- all$dr.flow.m2
            }
          }
        }
        if (tea == FALSE && save.statmat == TRUE) 
          out$stat.mat[[s]] <- all$stat.mat
        if (vital == TRUE && save.stats == TRUE) {
          out$stats[[s]] <- as.data.frame(attributes(all$nw)$stats)
        }
        if (save.trans == TRUE) {
          if (!is.null(all$transdf)) {
            row.names(all$transdf) <- 1:nrow(all$transdf)
            out$trans[[s]] <- all$transdf
          } else {
            out$trans[[s]] <- data.frame()
          }
        }
        if (save.network == TRUE) {
          out$network[[s]] <- all$nw
        }
      } # end sim > 1 out
      
  } # end sim loop  

  ## Final bookkeeping
  {
    # Set names for out
    simnames <- paste("sim", 1:tot.sims, sep="")
    for (i in as.vector(which(lapply(out, class) == "data.frame"))) {
      colnames(out[[i]]) <- simnames
    }
    if (tea == FALSE && save.statmat == TRUE) names(out$stat.mat) <- simnames
    if (vital == TRUE && save.stats == TRUE) names(out$stats) <- simnames
    if (save.trans == TRUE) names(out$trans) <- simnames[1:length(out$trans)]
    if (save.network == TRUE) names(out$network) <- simnames
    
    # If only 1 sim, then change output from df to vectors
    if (tot.sims == 1) {
      for (i in as.vector(which(lapply(out, class) == "data.frame"))) {
        out[[i]] <- out[[i]][,1]
      }
    }
    
    out$vital <- vital
    out$nsims <- tot.sims
    out$call <- match.call()
    
    class(out) <- "epiNet.simTrans"
    on.exit(par(ops))
    invisible(out)
  }
  
}

