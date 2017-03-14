##================================================================================
## This file is part of the evoper package - EvoPER
##
## (C)2016, 2017 Antonio Prestes Garcia <@>
## For license terms see DESCRIPTION and/or LICENSE
##
## @file: abm-saa.R
##
## This file contains the Simulated Annealing Algorithm (SAA) metaheuristic
## and its associated auxiliary functions.
##================================================================================


#' @title abm.saa
#'
#' @description An implementation of Simulated Annealing Algorithm
#' optimization method for parameter estimation of Individual-based
#' models.
#'
#' @param objective An instance of ObjectiveFunction (or subclass) class \link{ObjectiveFunction}
#' @param options An apropiate instance from a sublclass of \link{Options} class
#'
#' @return The best solution.
#'
#' @references
#'
#' [1] Kirkpatrick, S., Gelatt, C. D., & Vecchi, M. P. (1983).
#' Optimization by Simulated Annealing. Science, 220(4598).
#'
#' @examples \dontrun{
#'  ## A Repast defined function
#'  f<- RepastFunction$new("/usr/models/BactoSim(HaldaneEngine-1.0)","ds::Output",300)
#'
#'  ## or a plain function
#'
#'  f1<- function(x1,x2,x3,x4) {
#'    10 * (x1 - 1)^2 + 20 * (x2 - 2)^2 + 30 * (x3 - 3)^2 + 40 * (x4 - 4)^2
#'  }
#'
#'  f<- PlainFunction$new(f1)
#'
#'  f$addFactor(name="cyclePoint",min=0,max=90)
#'  f$addFactor(name="conjugationCost",min=0,max=100)
#'  f$addFactor(name="pilusExpressionCost",min=0,max=100)
#'  f$addFactor(name="gamma0",min=1,max=10)
#'
#'  abm.saa(f, 100, 1,  100, 0.75)
#' }
#'
#' @importFrom stats runif
#' @export
abm.saa<- function(objective, options= NULL) {
  ## Handling the heuristic specific options
  options<- OptionsFactory("saa", options)

  ## --- Creating the estimation object for returning results
  estimates<- Estimates$new()

  ## --- Configure algorithm parameters
  iterations<- options$getValue("iterations")
  t<- options$getValue("t0")
  TMIN<- options$getValue("t.min")
  L<- options$getValue("L")
  d<- options$getValue("d")
  f.neighborhood<- options$neighborhoodFunction()
  f.temp<- options$getTemperatureF()

  ## Generates an initial solution
  elog.info("Initializing solution")
  s<- initSolution(objective$parameters, 1)
  f<- objective$EvalFitness(s)
  S<- sortSolution(s, f)
  bestS<- S

  elog.info("Starting metaheuristic")
  for(iteration in 1:iterations) {

    ## --- Calculate the current temperature
    t<- f.temp(t, iteration)
    S<- bestS
    for(l in 1:L) {
      ## Evaluate a neighbor of s
      s1<- f.neighborhood(objective,getSolution(S), d)
      f1<- objective$EvalFitness(s1)
      S1<- sortSolution(s1, f1)

      ## --- Saving visited solution points
      estimates$addVisitedSpace(S1)

      ## --- New solution is beter than the previous one
      if(bestFitness(S1) < bestFitness(S)) {
        S<- S1
        if(bestFitness(S1) < bestFitness(bestS)) {
          bestS<- S1
        }
      } else {
        ## --- Calculate the cost delta
        delta<- bestFitness(S1) - bestFitness(S)
        if(runif(1) < exp(-delta/t)) {
          S<- S1
        }
      }
    }

    ## Show current bests
    elog.info("Iteration=%d/%d, best fitness=%g, iteration best fitness=%g", iteration, iterations, bestFitness(bestS), bestFitness(S))

    ## --- Storing the best of this iteration
    estimates$addIterationBest(iteration, S1[1,])

    if(t <= TMIN) break;
    if(objective$isConverged(bestFitness(bestS))) break
  }

  estimates$setBest(bestS)
  estimates
}


## ##################################################################
##
## --- Simulated Annealing Algorithm helper functions -
##
## ##################################################################


#' @title saa.neighborhood
#'
#' @description Generates neighbor solutions for simulated annealing
#'
#' @param f An instance of ObjectiveFunction (or subclass) class \link{ObjectiveFunction}
#' @param S The current solution to find a neighbor
#' @param d The distance from current solution S \code{distance = (max - min) * d}
#' @param n The number of parameters to be perturbed
#'
#' @return The neighbor of solution S
#'
#' @export
saa.neighborhood<- function(f, S, d, n) {
  assert(n > 0 && n <= ncol(S),"Invalid number of parameters to be perturbed!")
  newS<- S
  f.range<- f$getParameterRange

  for(i in sample(1:ncol(S),n)) {
    k<- colnames(S)[i]
    distance<- f.range(k) * d
    m.mean<- mean( c(as.numeric(f$getParameterValue(k,"max")), as.numeric(f$getParameterValue(k,"min"))) )

    #newS[,i]<- newS[,i] + runif(1,as.numeric(f$getParameterValue(k,"min")),as.numeric(f$getParameterValue(k,"max"))) * distance
    #88# newS[,i]<- newS[,i] + newS[,i] * runif(1,-1,1) * distance

    if(runif(1) < 1/5) {
      #newS[,i]<- newS[,i] + stats::rnorm(1)
      newS[,i]<- stats::rnorm(1,mean=m.mean,sd=distance)
    } else {
      #newS[,i]<- stats::rnorm(1,mean=m.mean,sd=distance)
      newS[,i]<- newS[,i] + newS[,i] * runif(1,-1,1)
    }


    #99# newS[,i]<- newS[,i] + .01 * f.range(k) * runif(1,0,1)
  }
  enforceBounds(as.data.frame(newS), f$parameters)
}


#' @title saa.neighborhood1
#'
#' @description Generates neighbor solutions perturbing one parameter from current
#' solution \code{S} picked randonly.
#'
#' @param f An instance of ObjectiveFunction (or subclass) class \link{ObjectiveFunction}
#' @param S The current solution to find a neighbor
#' @param d The distance from current solution S \code{distance = (max - min) * d}
#'
#' @return The neighbor of solution of S
#'
#' @export
saa.neighborhood1<- function(f, S, d) {
  saa.neighborhood(f, S, d, 1)
}

#' @title saa.neighborhoodH
#'
#' @description Generates neighbor solutions perturbing half parameters from current
#' solution \code{S}.
#'
#' @param f An instance of ObjectiveFunction (or subclass) class \link{ObjectiveFunction}
#' @param S The current solution to find a neighbor
#' @param d The distance from current solution S \code{distance = (max - min) * d}
#'
#' @return The neighbor of solution of S
#'
#' @export
saa.neighborhoodH<- function(f, S, d) {
  saa.neighborhood(f, S, d, floor(ncol(S)/2))
}

#' @title saa.neighborhoodN
#'
#' @description Generates neighbor solutions perturbing all parameters from current
#' solution \code{S}.
#'
#' @param f An instance of ObjectiveFunction (or subclass) class \link{ObjectiveFunction}
#' @param S The current solution to find a neighbor
#' @param d The distance from current solution S \code{distance = (max - min) * d}
#'
#' @return The neighbor of solution of S
#'
#' @export
saa.neighborhoodN<- function(f, S, d) {
  saa.neighborhood(f, S, d, ncol(S))
}

#' @title saa.tbyk
#'
#' @description Temperature function t/k
#'
#' @param t0 The current temperature
#' @param k The annealing value
#'
#' @return The new temperature
#'
#' @export
saa.tbyk<- function(t0, k) {
  t0 / k
}

#' @title saa.texp
#'
#' @description Temperature function exponential
#'
#' @param t0 The current temperature
#' @param k The annealing value
#'
#' @return The new temperature
#'
#' @export
saa.texp<- function(t0, k) {
  t0 * 0.95^k
}

#' @title saa.bolt
#'
#' @description Temperature function boltzmann
#'
#' @param t0 The current temperature
#' @param k The annealing value
#'
#' @return The new temperature
#'
#' @export
saa.bolt<- function(t0, k) {
  t0 / log(k)
}

#' @title saa.tcte
#'
#' @description Temperature function cte * t0
#'
#' @param t0 The current temperature
#' @param k The annealing value
#'
#' @return The new temperature
#'
#' @export
saa.tcte<- function(t0, k) {
  .95 * t0
}
