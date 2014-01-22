#' 
#' @title Concurrency Microsimulation Model for HIV-1 Transmission Dynamics
#'
#' @description This function simulates an HIV epidemic in a population of men 
#'  and women with purely heterosexual mixing under varying scenarios of sexual
#'  partnership concurrency.
#'
#' @param s.num.f number of initial susceptible females in the population.
#' @param i.num.f number of initial infected females in the population.
#' @param s.num.m number of initial susceptible males in the population.
#' @param i.num.m number of initial infected females in the population.
#' @param monog.f if \code{TRUE}, enforce a momentary degree constraint of monogamy 
#'  for females (females not allowed concurrent partnerships).
#' @param monog.m if \code{TRUE}, enforce a momentary degree constraint of monogamy
#'  for males (males not allowed concurrent partnerships).
#' @param meandeg average momentary mean degree (number of current partnerships)
#'  in the population.
#' @param part.duration average length of partnerships in months.
#' @param nsteps number of time steps to simulate the model over. This must be a 
#'  positive integer. 
#' @param nsims number of simulations to run.
#' @param verbose print model simulation progress to the console.
#'  
#' @details 
#' This is a microsimulation model of a dynamic network to see how the presence 
#' or absence of relational concurrency affects the prevalence of HIV prevalence
#' at the population level.
#' 
#' HIV infection is simulated based on a four-stage disease progression model in
#' which persons transition from acute to latent to pre-AIDS to AIDS stages. These
#' transitions occur at deterministic intervals based on estimates of the average
#' time per stage. Also, the transmission probability to uninfected partners varies
#' by stage of the infected partner: it is highest in the acute stage and lowest
#' in the AIDS stage when no sexual acts are simulated to occur. See the 
#' Hollingsworth reference below for further details.
#' 
#' The main parameters of the model include the initial number of susceptible 
#' females and males, the initial number of infected females and males, whether 
#' males and females are allowed concurrency, the mean degree of all persons in 
#' the population, and the average duration of partnerships (in momths). As the 
#' four examples below show, we can hold constant all parameters but toggle whether
#' men, women, or both sexes exhibit concurrency.
#' 
#' This model makes several simplifying assumptions about partnership formation 
#' and dissolution, such as the phenomenon of all dissolving partnerships being 
#' immediately replaced by a new partnership in the network. Additionally, the
#' user may specify whether concurrency is allowed, but not the level of concurrency 
#' (it is automatically calculated here based on a binomial distribution with the 
#' probability set by the mean degree parameter; see the full tutorial link below). 
#' Therefore, this model as an introduction to network modeling featured in the 
#' \code{epiNet} class of functions in \code{EpiModel}: there the user has much 
#' more control over the network parameterization and evolution.
#' 
#' @references 
#' A web-based implementation of this model is available at 
#' \url{http://statnet.org/apps/Conc}. The background and details of this model 
#' are explained in a full tutorial on concurrency at  
#' \url{https://statnet.csde.washington.edu/trac/wiki/ConcurrencyIndex}.
#' 
#' Hollingsworth TD, Anderson RM, Fraser C. HIV-1 transmission, by stage of 
#' infection. Journal of Infectious Diseases. 2008; 198(5): 687-693.
#'  
#' @keywords model
#' @export
#' 
#' @examples
#' \dontrun{
#' # No concurrency model
#' no.conc <- conc.microsim(
#'    s.num.f = 1000, 
#'    i.num.f = 50, 
#'    s.num.m = 1000, 
#'    i.num.m = 50, 
#'    monog.f = TRUE, 
#'    monog.m = TRUE,
#'    meandeg = 0.8, 
#'    part.duration = 10, 
#'    nsteps = 2000, 
#'    nsims = 10, 
#'    verbose = TRUE)
#' 
#' # Male concurrency only model
#' male.conc <- conc.microsim(
#'    s.num.f = 1000, 
#'    i.num.f = 50, 
#'    s.num.m = 1000, 
#'    i.num.m = 50, 
#'    monog.f = TRUE, 
#'    monog.m = FALSE,
#'    meandeg = 0.8, 
#'    part.duration = 10, 
#'    nsteps = 2000, 
#'    nsims = 10, 
#'    verbose = TRUE)
#'    
#' # Female concurrency only model
#' feml.conc <- conc.microsim(
#'    s.num.f = 1000, 
#'    i.num.f = 50, 
#'    s.num.m = 1000, 
#'    i.num.m = 50, 
#'    monog.f = FALSE, 
#'    monog.m = TRUE,
#'    meandeg = 0.8, 
#'    part.duration = 10, 
#'    nsteps = 2000, 
#'    nsims = 10, 
#'    verbose = TRUE)  
#'   
#' # Both sexes concurrency model 
#' both.conc <- conc.microsim(
#'    s.num.f = 1000, 
#'    i.num.f = 50, 
#'    s.num.m = 1000, 
#'    i.num.m = 50, 
#'    monog.f = FALSE, 
#'    monog.m = FALSE,
#'    meandeg = 0.8, 
#'    part.duration = 10, 
#'    nsteps = 2000, 
#'    nsims = 10, 
#'    verbose = TRUE)
#'    
#' # Plot the results
#' par(mfrow=c(2,2), mar=c(3,3,3,1), mgp=c(2,1,0))
#' plot(no.conc, alpha=0.5, ylim=c(0, 0.5), main="No Concurrency")
#' plot(male.conc, alpha=0.5, ylim=c(0, 0.5), main="Male Concurrency")
#' plot(feml.conc, alpha=0.5, ylim=c(0, 0.5), main="Female Concurrency")
#' plot(both.conc, alpha=0.5, ylim=c(0, 0.5), main="Both Concurrency")
#'}
#'
conc.microsim <- function(s.num.f, 
                          i.num.f, 
                          s.num.m, 
                          i.num.m, 
                          monog.f = TRUE, 
                          monog.m = TRUE,
                          meandeg, 
                          part.duration, 
                          nsteps, 
                          nsims = 1, 
                          verbose = TRUE) {
  
  
  # Probability of transmission per month for active relationship
  # From Hollingsworth et al. 2008
  beta.by.time.since.inf <- c(rep(0.2055, 3), rep(0.0088, 100),
                              rep(0.0614, 9), rep(0, 10))

  
  out <- list()
  for (sim in 1:nsims){
    ptm <- proc.time()
    
    # Basic calculations ------------------------------------------------------
    
    # Total pop size
    n.femls <- s.num.f + i.num.f
    n.males <- s.num.m + i.num.m
    n.pop <- n.femls + n.males					
    
    # Expected # of relationships ("edges") in the population at any time
    expected.edges <- round(meandeg * n.pop / 2)				
    
    # Daily prob of dissolution for an existing relationship
    prob.dissolution <- 1 / part.duration		
    
    # Length of time from HIV infection until death (in months)
    time.of.aids.death <- length(beta.by.time.since.inf)					
    
    # Cumulative # of female AIDS deaths
    cum.num.feml.aids.deaths <- 0				
    
    # Cumulative # of female AIDS deaths
    cum.num.male.aids.deaths <- 0											
    
    # Sets up vectors to store prevalence at each time step
    feml.prev <- vector()		
    male.prev <- vector()														
    
    
    # Female data frame -------------------------------------------------------
    
    # Creates a data frame to store info about the females.  One row per female.
    femls <- data.frame(row.names = 1:n.femls)					
    
    # Create a variable in the female data frame called "hiv.status"
    # Give all females the value 0 (for the moment).
    femls$hiv.status <- rep(0, n.femls)										
    
    # Create a variable in the female data frame called "inf.time" 
    #   (infection time). Give all females an NA (for the moment).
    femls$inf.time <- NA							
    
    # If there are any initially infected females,
    # give them a value of 1 for hiv.status,
    # and an infection time in the past
    if (i.num.f > 0) {		
      # Sample from vector of female ids with size of initially infected females
      init.inf.f <- sample(1:n.femls, i.num.f)
      femls$hiv.status[init.inf.f] <- 1				
      
      # Sample times backwards from present to length of infection
      init.inf.time.f <- sample(0:(-time.of.aids.death+2), i.num.f, replace=TRUE)
      femls$inf.time[init.inf.f] <- init.inf.time.f	
    }
    
    
    # Male data frame ---------------------------------------------------------
    
    # All parallel to female data frame above
    males <- data.frame(row.names = 1:n.males)
    males$hiv.status <- 0
    males$inf.time <- NA
    if (i.num.m > 0) {
      init.inf.m <- sample(1:n.males, i.num.m)
      males$hiv.status[init.inf.m] <- 1 
      init.inf.time.m <- sample(0:(-time.of.aids.death+2), i.num.m, replace=TRUE)
      males$inf.time[init.inf.m] <- init.inf.time.m
    }
    
    
    # Initial contact network -------------------------------------------------
    
    # Create a data frame that will be used to hold the IDs of the relationship pairs
    edgelist <- data.frame(row.names = 1:expected.edges)							
    # If we are *not* enforcing female monogamy,
    # sample from the female IDs with replacement, and assign them to a column 
    # in the data frame called "f"
    if (monog.f == FALSE) {
      edgelist$f <- sample(1:n.femls, expected.edges, replace = TRUE)				
      #  otherwise,	sample from the female IDs without replacement, and assign 
      #   them to a column in the data frame called "f"
    } else {
      edgelist$f <- sample(1:n.femls, expected.edges, replace = FALSE)																			
    }
    
    # Same for males
    if (monog.m == FALSE) {
      edgelist$m <- sample(1:n.males, expected.edges, replace = TRUE)
    } else {
      edgelist$m <- sample(1:n.males, expected.edges, replace = FALSE)
    }
    
    # These rather odd lines provide a check to see a table of relationships 
    #   per person in the population. When monogamy is enforced for
    #   either sex, the corresponding table should be limited to 0's and 1's; 
    #   when it os not enforced, the table is not limited.
    #table(tabulate(edgelist$f))													
    #table(tabulate(edgelist$m))
    
    
    # Time loop ---------------------------------------------------------------
    
    for (time in 1:nsteps) {
      
      ## Transmissions ##
      # Create a vector containing the HIV status for the female partner in each 
      #  relationship
      hiv.status.feml.partners <- femls$hiv.status[edgelist$f]
      # Create a vector containing the HIV status for the male partner in each 
      #  relationship
      hiv.status.male.partners <- males$hiv.status[edgelist$m]
      
      # Make a vector of IDs for the relationships that are serodiscordant with 
      #  positive female (SDPF)
      sdpf <- which(hiv.status.male.partners == 0 &
                    hiv.status.feml.partners == 1)
      
      # Make a vector of IDs for the females in the SDPF relationships
      f.in.sdpf <- edgelist$f[sdpf]
      
      # Make a vector that identifies how long these females were infected
      inf.time.sdpf <- time - femls$inf.time[f.in.sdpf]		
      
      # Determine probability of transmission in the SDPM relationships, based 
      #  on time since the male was infected
      prob.trans.sdpf <- beta.by.time.since.inf[inf.time.sdpf]
      
      # Flip a weighted coin for each SDPM relationship to determine if 
      #  transmission occurs
      trans.sdpf <- rbinom(length(sdpf), 1, prob.trans.sdpf)
      
      # Double-indexing (an R specialty!)  This lets you get the identities of 
      #  the newly infected males in one step.
      newly.inf.males <- edgelist$m[sdpf[trans.sdpf == 1]]
      
      # Assign the newly infected males an hiv.status of 1
      males$hiv.status[newly.inf.males] <- 1
      
      # Assign the newly infected males an infection time of time
      males$inf.time[newly.inf.males] <- time
      
      # All parallel for serodiscordant with positive male (SDPM)
      sdpm <- which(hiv.status.feml.partners == 0 &
                    hiv.status.male.partners == 1)
      m.in.sdpm <- edgelist$m[sdpm]
      inf.time.sdpm <- time - males$inf.time[m.in.sdpm]
      prob.trans.sdpm <- beta.by.time.since.inf[inf.time.sdpm]
      trans.sdpm <- rbinom(length(sdpm), 1, prob.trans.sdpm)
      newly.inf.femls <- edgelist$f[sdpm[trans.sdpm == 1]]
      femls$hiv.status[newly.inf.femls] <- 1
      femls$inf.time[newly.inf.femls] <- time
      
      
      ## Deaths to AIDS ##
      # Which females have been HIV+ long enough to die of AIDS?
      femls.dying.of.AIDS <- which((time - femls$inf.time) == time.of.aids.death)
      # Which females have been HIV+ long enough to die of AIDS?
      males.dying.of.AIDS <- which((time - males$inf.time) == time.of.aids.death) 
      
      # Increase cumulative # of female AIDS deaths by the newly dying females
      cum.num.feml.aids.deaths <- cum.num.feml.aids.deaths + 
                                  length(femls.dying.of.AIDS)
      # Increase cumulative # of male AIDS deaths by the newly dying males
      cum.num.male.aids.deaths <- cum.num.male.aids.deaths + 
                                  length(males.dying.of.AIDS)
      
      ## End of ties because of death ##
      # Determine the IDs of those relationships involving dying women
      edges.with.feml.dying.of.aids <- which(edgelist$f %in% femls.dying.of.AIDS)
      # Determine the IDs of those relationships involving dying women
      edges.with.male.dying.of.aids <- which(edgelist$m %in% males.dying.of.AIDS)					
      # Combine the two
      edges.with.either.dying.of.aids <- c(edges.with.feml.dying.of.aids, 
                                           edges.with.male.dying.of.aids)	
      # Remove those edges from the edgelist
      if (length(edges.with.either.dying.of.aids > 0)){						
        edgelist <- edgelist[-edges.with.either.dying.of.aids,]
      }
      
      
      ## End of other ties randomly ##
      # Flip a weighted coin for each remaning relationship
      edges.coinflip <- rbinom(nrow(edgelist), 1, prob.dissolution)	
      # Make vector of those edges to break
      edges.to.break <- which(edges.coinflip == 1)						
      # If there are any edges to break, then break them.
      if (length(edges.to.break > 0)) {
        edgelist <- edgelist[-edges.to.break, ]
      }
      
      
      ## Add new edges ##
      # Determine how many edges to add. 
      # We assume same as # broken in this simple model.
      num.ties.to.add <- expected.edges - nrow(edgelist)				
      
      # If there are edges to add,
      if (num.ties.to.add > 0) {
        #  and if we are *not* enforcing female monogamy,
        if (monog.f == FALSE) {	
          # sample from the female IDs with replacement,
          f <- sample(1:n.femls, num.ties.to.add, replace = TRUE)		
          #	else if we are enforcing female monogamy,
        } else {											
          # determine which women do not currently have a relationship, and
          femls.with.no.ties <- setdiff(1:n.femls, edgelist$f)		
          # sample from them without replacement
          f <- sample(femls.with.no.ties, num.ties.to.add, replace = FALSE)	
        }
        # and if we are *not* enforcing male monogamy,
        if (monog.m == FALSE) {
          # sample from the male IDs with replacement,
          m <- sample(1:n.males, num.ties.to.add, replace = TRUE)
          # else if we are enforcing male monogamy,  
        } else {
          # determine which men do not currently have a relationship, and
          males.with.no.ties <- setdiff(1:n.males, edgelist$m)		
          # sample from them without replacement
          m <- sample(males.with.no.ties, num.ties.to.add, replace = FALSE)	
        }
        # combine the new males and females into pairs in a data frame
        new.edges <- data.frame(f, m)
        # "Bind" that data frame onto the end of our existing esgelist
        edgelist <- data.frame(rbind(edgelist, new.edges))		
        # Rename the rows of the updated edgelist so they are consecutive 
        #  numbers again.
        row.names(edgelist) <- 1:nrow(edgelist)					
      }
      
      
      ## Insert new births ##
      # Make a vector of the IDs of the women who just died. We shall "replace" 
      #  them with a new arrival, by
      new.feml.ids <- femls.dying.of.AIDS
      #   setting the HIV status of that node back to 0, and
      femls$hiv.status[new.feml.ids] <- 0		
      #   setting the infection time for that node back to NA.
      femls$inf.time[new.feml.ids] <- NA								
      
      # All parallel with males
      new.male.ids <- males.dying.of.AIDS								
      males$hiv.status[new.male.ids] <- 0
      males$inf.time[new.male.ids] <- NA
          
      
      ## Track prevalence ##
      # Calculate current female prevalence and append it onto the female 
      #  prevalence vector.
      feml.prev[time] <- mean(femls$hiv.status)
      # Calculate current male prevalence and append it onto the male prevalence 
      #  vector.
      male.prev[time] <- mean(males$hiv.status)
    }    
    # Progress tracker
    time <- (proc.time()[3] - ptm[3])
    done <- (proc.time()[3] - ptm[3]) * (nsims-sim)
    if (verbose == TRUE) {
      if (done/60 < 60) {
        cat("SIM = ", sim, "/", nsims,
            " | SIM Time = ", round(time/60, 1), " Min",
            " | SIM Done = ", round(done/60), " Min", "\n", sep="")
      } else {
        cat("SIM = ", sim, "/", nsims,
            " | SIM Time = ", round(time/60, 1), " Min",
            " | SIM Done = ", round(done/60/60, 1), " Hrs", "\n", sep="")
      }
    }
    
    df <- data.frame(feml.prev, male.prev)
    out[[sim]] <- df
  
    } # end sim loop
  
  class(out) <- "conc.microsim"
  return(out)
}



#' @title Plot Values from a Concurrency Microsimulation Model
#' 
#' @description This function plots values from an the concurrency microsimulation 
#' epidemic model simulated with \code{conc.microsim}.
#'
#' @param x an \code{EpiModel} object of class \code{conc.microsim}.
#' @param xlim x-axis scale limits for plot, with default calculated based on 
#'   model time steps.
#' @param ylim y-axis scale limits for plot, with default calculated based on 
#'   range of data.
#' @param alpha transparency level for simulation lines, where 0 = transparent 
#'   and 1 = opaque (see \code{\link{transco}}).
#' @param ... additional arguments to pass to main plot 
#'   (see \code{\link{plot.default}}).
#' 
#' @details
#' This function will extract the simulated values output values from an 
#' \code{conc.microsim} model and plot relevant time series data of disease prevalence 
#' and other results. Currently, the plot function is limited to individual 
#' simulation lines of disease prevalence; future releases will standardize the 
#' plotting options to those available with \code{\link{epiICM}} class models.
#' 
#' @method plot conc.microsim
#' @export 
#' 
plot.conc.microsim <- function(x, 
                               xlim, 
                               ylim, 
                               alpha = 1, 
                               ...) {
  
  nsims <- length(x)
  nsteps <- nrow(x[[1]])
  pal <- transco(brewer.pal(3, "Set1"), alpha)
  
  if (missing(ylim)) ylim <- 0:1
  if (missing(xlim)) xlim <- c(0, nsteps)
  plot(1, 1, type="n", xlim = xlim, ylim = ylim, bty="n", 
       xlab="Time (months)", ylab="Infected", ...)
  for (i in 1:nsims) {
    lines(1:nsteps, x[[i]][,1], col=pal[1])
    lines(1:nsteps, x[[i]][,2], col=pal[2])
  }
  legend("topleft", c("Female", "Male"), lwd=3, col=pal[1:2], bty="n", cex=0.9)
  
}