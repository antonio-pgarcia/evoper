##================================================================================
## This file is part of the evoper package - EvoPER
##
## (C)2020 Antonio Prestes Garcia <@>
## For license terms see DESCRIPTION and/or LICENSE
##
## @file: abm-ga.R
##
## This file contains the Genetic Algorithm optimization metaheuristic and
## its associated auxiliary functions.
##================================================================================


#' @title Genetic Algorithm metaheuristic
#'
#' @description An implementation of Genetic Algorithm metaheuristic for parameter estimation
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
#' [1] John Henry Holland (1992). "Adaptation in Natural and Artificial Systems; An Introductory Analysis
#' with Applications to Biology, Control and Artificial Intelligence". MIT Press, Cambridge, MA, USA.
#' ISBN 0262082136.
#' [2] Zbigniew Michalewiczx (1994). "Genetic Algorithms + Data Structures = Evolution Programs (2nd Ed.)".
#' Springer-Verlag, Berlin, Heidelberg. ISBN 3540580905.
#'
#' @export
abm.ga<- function(objective, options= NULL) {
  ## Handling the heuristic specific options
  options<- OptionsFactory("ga", options)
  elog.info("Options(%s): %s", options$getType(), options$toString())

  ## --- Adjusting parameter types
  parameterz<- paramconverter(objective$parameters, options$isDiscrete(), options$getLevels())
  objective$parameters<- parameterz

  ## --- Creating the estimation object for returning results
  estimates<- Estimates$new()

}

