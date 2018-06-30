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

  ## --- Creating the estimation object for returning results
  estimates<- Estimates$new()

  ## --- Configure algorithm parameters
  iterations<- options$getValue("iterations")
  N<- options$getValue("N")
  tabusize<- options$getValue("tabu_size")
  candsize<- options$getValue("cand_size")

  ## Tabu list initialization
  tabu<- c()

  ## Generates an initial solution
  elog.info("Initializing solution")
  s<- generateSolution(objective$parameters, N)

  ## -- Evaluate
  s0<- es.evaluate(objective, s)
  push(tabu, s0)

  elog.info("Starting metaheuristic")
  for(iteration in 1:iterations) {

    s1<- generateSolution(objective$parameters, candsize)

    #sortSolution


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
