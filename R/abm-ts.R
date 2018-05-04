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
#' @references
#'
#' [1] Fred Glover (1989). "Tabu Search – Part 1". ORSA Journal on Computing,
#' 190–206. doi:10.1287/ijoc.1.3.190.
#' [2] Fred Glover (1990). "Tabu Search – Part 2". ORSA Journal on Computing,
#' 4–32. doi:10.1287/ijoc.2.1.4.
#'
#' @export
abm.acor<- function(objective, options= NULL) {
}
