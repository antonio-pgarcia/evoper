##================================================================================
## This file is part of the evoper package - EvoPER
##
## (C)2016, 2017 Antonio Prestes Garcia <@>
## For license terms see DESCRIPTION and/or LICENSE
##
## @file: class-objective.R
##
## This file contains the classes abstracting the objective function
##================================================================================


#' @title ObjectiveFunction class
#'
#' @description The base class for optimization functions.
#'
#' @field object The raw output of objective function
#' @field objective The objective function
#' @field parameters The parameter list for objective function
#' @field value The results from objective function
#'
#' @importFrom methods new
#' @export ObjectiveFunction
#' @exportClass ObjectiveFunction
ObjectiveFunction<- setRefClass("ObjectiveFunction",

  fields = list(
    object = 'ANY',
    objective = 'function',
    replicates = 'numeric',
    objective.defaults = 'ANY',
    parameters = 'ANY',
    value = 'ANY',
    rawdata = 'ANY',
    tolerance = 'ANY',
    converged = 'ANY',
    maximize = 'ANY',
    counter = 'ANY'),

  methods = list(

    initialize = function(funct) {
      object<<- NULL
      objective<<- funct
      replicates<<- 1
      objective.defaults<<- NULL
      parameters<<- NULL
      value<<- NULL
      rawdata<<- NULL
      tolerance<<- .Machine$double.eps^0.30
      converged<<- FALSE
      maximize<<- FALSE
      counter<<- 0
    },

    stats = function() {
      cbind(total_evals=counter,converged=converged,tolerance=tolerance)
    },

    isConverged = function(v) {
      if(!converged) {
        converged<<- (v <= tolerance)
      }
      converged
    },

    setTolerance = function(v) {
      tolerance<<- v
    },

    setMaximize = function(v=FALSE) {
      maximize<<- v
    },

    Parameter0 = function(...) {
      assert(hasArg(name) && ((hasArg(min) && hasArg(min)) || hasArg(levels)), "Please provide the required parameters")
      if(is.null(parameters)) {
        parameters<<- rrepast::AddFactor0(factors=c(), ...)
      } else {
        parameters<<- rrepast::AddFactor0(factors=parameters, ...)
      }
    },

    Parameter = function(name, min, max, forceint=FALSE) {
      assert(hasArg(name) && hasArg(min) && hasArg(min), "Please provide the required parameters")
      if(is.null(parameters)) {
        parameters<<- rrepast::AddFactor(c(), name= name, min= min, max= max, int= forceint)
      } else {
        parameters<<- rrepast::AddFactor(parameters, name= name, min= min, max= max, int= forceint)
      }
    },

    getParameterRange = function(key) {
      as.numeric(parameters[which(parameters[,"name"] == key),"max"]) - as.numeric(parameters[which(parameters[,"name"] == key),"min"])
    },

    getParameter = function(key) {
      p<- parameters[which(parameters[,"name"] == key),]
      p<- as.data.frame(as.list(p))

      for(k in names(p)) {
        if(is.factor(p[,k])) {
          v<- levels(p[,k])
        } else {
          v<- p[,k]
        }

        if(k %in% c("min","max")) {
          p[,k]<- as.numeric(v)
        } else {
          p[,k]<- v
        }
      }
        p
    },

    getParameterNames = function() {
      parameters[,"name"]
    },

    getParameterValue = function(key, name) {
      parameters[which(parameters[,"name"] == key),name]
    },

    getParameters = function() {
      p<- data.frame(parameters,  stringsAsFactors=FALSE)

      for(i in 1:length(p[,1])) {
        for(k in names(p)) {
          v<- p[,k]
          if(k %in% c("min","max")) {
            p[,k]<- as.numeric(v)
          } else {
            p[,k]<- v
          }
        }
      }
      p
    },

    ParametersSize = function() {
      length(parameters[,1])
    },


    defaults = function(v = NULL) {
      if(!is.null(v)) {
        objective.defaults<<- v
      } else {
        objective.defaults
      }
    },

    RawData = function(v = NULL) {
      if(!is.null(v)) {
        rawdata<<- v
      }
      rawdata
    },

    Replicates = function(v= NULL) {
      if(!is.null(v)) {
          replicates<<- v
      }
      replicates
    },

    Value = function(v = NULL) {
      if(!is.null(v)) {
        if(maximize) {
          v<- 1/v
        }
        value<<- v
      }
      value
    },

    Evaluate = function(swarm) {
      assert(!is.null(objective),"The objective has not been provided!")
      assert(!is.null(swarm),"Solution must not be null!")

      ## Stats: Update the call counter
      counter<<- counter + nrow(swarm)
    },

    ## Evaluate population and return results
    EvalFitness = function(swarm) {
      Evaluate(swarm)
      Value()
    }

  )
)

#' @title PlainFunction
#'
#' @description PlainFunction Class
#'
#' @importFrom methods new
#' @export PlainFunction
#' @exportClass PlainFunction
PlainFunction<- setRefClass("PlainFunction", contains = "ObjectiveFunction",

  methods = list(

    initialize = function(o= NULL) {
      assert(!is.null(o),"NUll function!")
      assert(is.function(o),"Invalid function!")

      if(is.null(o)) {
        o<- objective
      }
      callSuper(o)
    },

    Evaluate = function(swarm) {
      callSuper(swarm)
      v<- Value(as.data.frame(t(apply(expand.grid(j=1,i=1:nrow(swarm)),1,function(k) { c(pset=as.integer(k[2]),fitness=do.call(objective,swarm[k[2],])) }))))
    }

  )
)

#' @title RepastFunction
#'
#' @description RepastFunction class
#'
#' @importFrom methods new
#' @export RepastFunction
#' @exportClass RepastFunction
RepastFunction<- setRefClass("RepastFunction", contains = "ObjectiveFunction",

  fields = list(
    model = 'ANY',
    directory = 'character',
    datasource = 'character',
    endAt = 'numeric'
  ),

  methods = list(

    initialize = function(d= NULL, ds= NULL, t= NULL, o= NULL) {
      if(is.null(o)) {
        o<- objective
      } else {
        directory<<- d
        datasource<<- ds
        endAt<<- t
      }
      callSuper(o)
    },

    Evaluate = function(swarm) {
      callSuper(swarm)

      ## Initialization of Repast Model instance
      rrepast::Easy.Setup(directory)
      ##model<<- rrepast::Model(directory,endAt,datasource,TRUE)
      my.model<- rrepast::Model(directory,endAt,datasource,TRUE)

      ## --- Update if needed the default parameters
      if(!is.null(defaults())) {
        rrepast::UpdateDefaultParameters(my.model, defaults())
      }


      ## Building up parameter set
      p<- rrepast::GetSimulationParameters(my.model)

      if(!is.null(swarm)){
        tmp<- p
        p<- rrepast::BuildParameterSet(swarm, tmp)
      }

      ## Evaluate model
      object<<- RunExperiment(my.model,r=replicates,p, objective)

      ## Clean UP
      Engine.Finish(my.model)

      ## Sum the objective output and change the column name
      object$output<<- col.sum(object$output)
      n<- names(object$output)
      names(object$output)<<- replace(n, which(n == "total"),c("fitness"))

      ## Store the available model's raw data
      RawData(object)
      Value(object$output)
    },

    show = function() {
      print(paste("Model directory is .... [",directory,"]"))
      print(paste("Model datasource is ... [",datasource,"]"))
      print(paste("Simulation time is .... [",endAt,"]"))
    })
)
