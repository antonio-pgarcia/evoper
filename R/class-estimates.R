##================================================================================
## This file is part of the evoper package - EvoPER
##
## (C)2016, 2017 Antonio Prestes Garcia <@>
## For license terms see DESCRIPTION and/or LICENSE
##
## @file: class-estimates.R
##
## This file contains the classes abstracting the return type for all
## optimization metaheuristics
##================================================================================

#' @title Estimates
#'
#' @description A simple class for encapsulating the return of metaheuristic
#' methods
#'
#' @importFrom data.table rbindlist
#' @importFrom methods new
#' @importFrom boot boot boot.ci
#' @export Estimates
#' @exportClass Estimates
Estimates<- setRefClass("Estimates",

  fields = list(
    ci = 'ANY',
    overall.best = 'ANY',
    iteration.best = 'ANY',
    visited.space = 'ANY'
  ),


  methods = list(

    initialize = function() {
      overall.best<<- Inf
      iteration.best<<- c()
      visited.space<<- c()
    },

    setBest = function(v) {
      overall.best<<- v
    },

    getBest = function() {
      overall.best
    },

    addIterationBest = function(iteration, solution) {
      #iteration.best<<- rbind(iteration.best, c(iteration,solution))
      if(length(iteration.best) == 0) {
        iteration.best<<- c(setNames(iteration,"iteration"),solution)
      } else {
        iteration.best<<- rbindlist(list(iteration.best, c(setNames(iteration,"iteration"),solution)), use.names=TRUE)
      }

    },

    getIterationBest = function() {
      data.frame(iteration.best, stringsAsFactors=FALSE)
    },

    addVisitedSpace = function(solution) {
      m<- nrow(solution)

      if(m > 1) {
        for(i in 1:m) {
          visited.space<<- rbindlist(list(visited.space,solution[i,]), use.names=TRUE)
          #visited.space<<- rbind(visited.space, solution[i,])
        }
      } else {
        visited.space<<- rbindlist(list(visited.space, solution), use.names=TRUE)
        #visited.space<<- rbind(visited.space, solution)
      }

    },

    getVisitedSpace = function() {
      data.frame(visited.space, stringsAsFactors=FALSE)
    }


  )
)

