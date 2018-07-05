##================================================================================
## This file is part of the evoper package - EvoPER
##
## (C)2016, 2017, 2018 Antonio Prestes Garcia <@>
## For license terms see DESCRIPTION and/or LICENSE
##
## @file: abm-ts.R
##
## This file contains the Tabu Search optimization metaheuristic and
## its associated auxiliary functions.
##================================================================================


#' @title Tabu Search metaheuristic
#'
#' @description An implementation of Tabu Search algorithm for parameter estimation
#'
#' @param objective An instance of ObjectiveFunction (or subclass) class \link{ObjectiveFunction}
#' @param options An apropiate instance from a sublclass of \link{Options} class
#'
#' @examples \dontrun{
#'  f<- PlainFunction$new(f0.rosenbrock2)
#'
#'  f$Parameter(name="x1",min=-100,max=100)
#'  f$Parameter(name="x2",min=-100,max=100)
#'
#'  or
#'
#'  f$Parameter0(name="x1",levels=c(0:4))
#'  f$Parameter0(name="x2",levels=c(-2,-1,0,1,2))
#'
#'  extremize("tabu", f)
#' }
#'
#' @references
#'
#' [1] Fred Glover (1989). "Tabu Search – Part 1". ORSA Journal on Computing,
#' 190–206. doi:10.1287/ijoc.1.3.190.
#' [2] Fred Glover (1990). "Tabu Search – Part 2". ORSA Journal on Computing,
#' 4–32. doi:10.1287/ijoc.2.1.4.
#'
#' @export
abm.tabu<- function(objective, options= NULL) {
  ## Handling the heuristic specific options
  options<- OptionsFactory("tabu", options)
  elog.info("Options(%s): %s", options$getType(), options$toString())

  ## --- Adjusting parameter types
  parameterz<- paramconverter(objective$parameters, options$isDiscrete(), options$getLevels())
  objective$parameters<- parameterz

  ## --- Creating the estimation object for returning results
  estimates<- Estimates$new()

  ## --- Configure algorithm parameters
  iterations<- options$getValue("iterations")
  N<- options$getValue("N")
  tabusize<- options$getValue("tabu_size")


  ## Tabu list initialization
  tabu<- c()

  ## Generates an initial solution
  elog.info("Initializing solution")
  s<- generateSolution(parameterz, 1)

  ## -- Evaluate
  s0<- es.evaluate(objective, s)
  push(tabu, getSolution(s0))

  elog.info("Starting metaheuristic")
  for(iteration in 1:iterations) {
    candidates<- tabu.getNeighbors(tabu, parameterz, getSolution(s0), N)

    s1<- es.evaluate(objective, candidates)

    if(s1[1, "fitness"] < s0[1, "fitness"]) {
      s0<- s1[1,]
      push(tabu, getSolution(s0))
      if(nrow(tabu) > tabusize) {
        v<- pop.first(tabu)
      }
    }

    ## --- Storing the best of this iteration
    estimates$addIterationBest(iteration, s1[1,])

    ## --- Save the complete visited space
    estimates$addVisitedSpace(s1)

    elog.info("Iteration=%d/%d, fitness=%g, iteration fitness=%g", iteration, iterations, s0[1,"fitness"], s1[1,"fitness"])

    ## Check for algorithm convergence
    if(objective$isConverged(s0[1, "fitness"])) break
  }

  estimates$setBest(s0[1,])
  estimates
}

#' @title tabu.getNeighbors
#'
#' @description create neighbor solutions
#'
#' @param parameters The parameter set
#' @param solution The current solution
#' @param size The neigborhood size
#'
#' @return
#'
#' @export
tabu.getNeighbors<- function(tabu, parameters, solution, size) {
  nnames<- unlist(parameters[,"name"])
  mylevels<- c()
  for(nname in nnames) {
    mylevels<- rbind(mylevels, GetFactorLevels(parameters, nname))
  }

  current<- solution
  neighbors<- c()
  for(i in 1:size) {

    repeat {
      neighbor<- c()
      for(p in 1:ncol(solution)) {
        if(runif(1) < 0.25) {
          #index<- sample(c(1:length(mylevels[p,])),1)
          delta<- trunc(length(mylevels[p,]) * 0.5)
          index<- which(mylevels[p,] == as.numeric(current[p])) + sample(c(-delta,0,delta),1)
          index<- ifelse(index < 1, 1, ifelse(index > length(mylevels[p,]), length(mylevels[p,]), index))
        } else {
          index<- which(mylevels[p,] == as.numeric(current[p])) + sample(c(-1,0,1),1)
          index<- ifelse(index < 1, 1, ifelse(index > length(mylevels[p,]), length(mylevels[p,]), index))
        }
        neighbor<- cbind(neighbor, mylevels[p, index])
      }
      if(!tabu.istabu(tabu, neighbor) ) break
    }

    current<- neighbor
    neighbors<- rbind(neighbors, neighbor)
  }
  neighbors<- as.data.frame(neighbors)
  names(neighbors)<- names(solution)
  neighbors
}


#' @title tabu.istabu
#'
#' @description Check whether a solution is present on tabulist
#'
#' @param tabulist The matrix of tabu solutions
#' @param solution The solution value to be checked
#'
#' @return Boolean TRUE tabulist contains the solution
#'
#' @export
tabu.istabu<- function(tabulist, solution) {
    any(searchrow(tabulist, solution))
}

