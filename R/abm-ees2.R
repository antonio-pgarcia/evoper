##================================================================================
## This file is part of the evoper package - EvoPER
##
## (C)2016, 2017 Antonio Prestes Garcia <@>
## For license terms see DESCRIPTION and/or LICENSE
##
## @file: abm-ees2.R
##
## This file contains the EvoPER Evolutionary Strategy 1 (EES2) metaheuristic
##================================================================================


#' @title EvoPER Evolutionary Strategy 2
#'
#' @description This function tries to provide a rough approximation to
#' best solution when no information is available for the correct range
#' of input parameters for the objective function. It can useful for
#' studying the behavior of individual-based models with high
#' variability in the output variables showing nonlinear behaviors.
#'
#' @param objective An instance of ObjectiveFunction (or subclass) class \link{ObjectiveFunction}
#' @param options An apropiate instance from a sublclass of \link{Options} class
#'
#' @export
abm.ees2<- function(objective, options= NULL) {
  ## Handling the heuristic specific options
  options<- OptionsFactory("ees2", options)
  elog.info("Options(%s): %s", options$getType(), options$toString())

  ## --- Creating the estimation object for returning results
  estimates<- Estimates$new()

  ## --- Metaheuristic options
  N<- options$getValue("N")
  rho<- options$getValue("rho")
  iterations<- options$getValue("iterations")

  n<- trunc(N*rho)

  elog.info("Initializing solution")
  parameters<- objective$getParameters()
  s<- initSolution(parameters, N, "lhs")

  ## -- Evaluate
  s0<- es.evaluate(objective, s)

  elog.info("Starting metaheuristic")
  for(iteration in 1:(iterations-1)) {
    s<- data.frame(getSolution(s0)[1:n,], stringsAsFactors=FALSE)
    mmin<- apply(s,2,min)
    mmax<- apply(s,2,max)
    mmean<- apply(s,2,gm.mean)
    ssd<- apply(s,2,sd)
    interval<-  abs(mmax-mmin)/2

    for(k in parameters$name) {
      if(runif(1) < 1/5) {
        parameters[which(parameters$name == k),"min"]<- as.numeric((mmean-interval-runif(1,0,1))[k])
        parameters[which(parameters$name == k),"max"]<- as.numeric((mmean+interval+runif(1,0,1))[k])
      } else {
        parameters[which(parameters$name == k),"min"]<- as.numeric(( mmean - interval * runif(1,1,3))[k])
        parameters[which(parameters$name == k),"max"]<- as.numeric(( mmean + interval * runif(1,1,3))[k])
      }
    }

    s<- initSolution(parameters, N, "lhs")

    ## -- Evaluate
    s1<- es.evaluate(objective, s)

    ## -- Save the best iteration value
    estimates$addIterationBest(iteration, s1[1,])

    s0<- rbind(s0,s1)
    s0<- s0[with(s0,order(fitness)),]

    elog.info("Iteration=%d/%d, fitness=%g, iteration fitness=%g", iteration, iterations, s0[1,"fitness"], s1[1,"fitness"])

  }

  ## -- Saving the whole solution space
  estimates$addVisitedSpace(s0)
  #for(i in 1:length(s0[,1])) {
  #estimates$addVisitedSpace(s0[i,])
  #}

  estimates$setBest(s0[1,])
  estimates
}
