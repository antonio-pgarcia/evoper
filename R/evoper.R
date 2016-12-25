##================================================================================
## This file is part of the evoper package - EvoPER
##
## (C)2016 Antonio Prestes Garcia <@>
## For license terms see DESCRIPTION and/or LICENSE
##
## $Id$
##================================================================================


# ------------------------------------------------------------
# .onLoad, Hook for loading package namespace
#
# ------------------------------------------------------------
.onLoad<- function(libname, pkgname) {
  assign("pkg.globals", new.env(), envir=parent.env(environment()))

  # The Random Seed. You may want to change this.
  set.seed(exp(1)*10^6)

  # Internal variables
  assign("pkg.basedir", NA, pkg.globals)
}

#' @title assert
#'
#' @description The assert function stop the execution if the logical
#' expression given by the parameter \code{expresion} is false.
#'
#' @param expresion Some logical expression
#' @param string The text message to show if expression does not hold
#'
#' @export
assert<- function(expresion, string) {
  if(!expresion) stop(string, call. = FALSE)
}


##
## ----- Wrapping the logging system
##

#' @title elog.level
#'
#' @description Configure the current log level
#'
#' @param level The log level (ERROR|WARN|INFO|DEBUG)
#'
#' @return The log level
#'
#' @importFrom futile.logger flog.threshold
#' @export
elog.level<- function(level=NULL) {
  if(!is.null(level)) {
    flog.threshold(level)
  }
  flog.threshold()
}

#' @title elog.error
#'
#' @description Wrapper for logging error messages.
#'
#' @param ... Variable number of arguments including a format string.
#'
#' @importFrom futile.logger flog.error
#' @export
elog.error<- function(...) {
  flog.error(...)
}

#' @title elog.info
#'
#' @description Wrapper for logging info messages.
#'
#' @param ... Variable number of arguments including a format string.
#'
#' @importFrom futile.logger flog.info
#' @export
elog.info<- function(...) {
  flog.info(...)
}

#' @title elog.debug
#'
#' @description Wrapper for logging debug messages.
#'
#' @param ... Variable number of arguments including a format string.
#'
#' @importFrom futile.logger flog.debug
#' @export
elog.debug<- function(...) {
  flog.debug(...)
}


##
## ----- Classes for the objective function
##


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

    Parameter = function(name, min, max) {
      assert(hasArg(name) && hasArg(min) && hasArg(min), "Please provide the required parameters")
      if(is.null(parameters)) {
        parameters<<- rrepast::AddFactor(c(), name= name, min= min, max= max)
      } else {
        parameters<<- rrepast::AddFactor(parameters, name= name, min= min, max= max)
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
      assert(!is.null(swarm),"Swarm must not be null!")

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
      model<<- rrepast::Model(directory,endAt,datasource,TRUE)

      ## Building up parameter set
      p<- rrepast::GetSimulationParameters(model)
      if(!is.null(swarm)){
        tmp<- p
        p<- rrepast::BuildParameterSet(swarm, tmp)
      }

      ## Evaluate model
      object<<- RunExperiment(model,r=replicates,p, objective)

      # Sum the objective output and change the column name
      object$output<<- col.sum(object$output)
      n<- names(object$output)
      names(object$output)<<- replace(n, which(n == "total"),c("fitness"))

      RawData(object)
      Value(object$output)
    },

    show = function() {
      print(paste("Model directory is .... [",directory,"]"))
      print(paste("Model datasource is ... [",datasource,"]"))
      print(paste("Simulation time is .... [",endAt,"]"))
    })
)


#' @title Estimates
#'
#' @description Estimates class
#'
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
      iteration.best<<- rbind(iteration.best, c(iteration,solution))
    },

    getIterationBest = function() {
      data.frame(iteration.best, stringsAsFactors=FALSE)
    },

    addVisitedSpace = function(solution) {
      m<- nrow(solution)

      if(m > 1) {
        for(i in 1:m) {
          visited.space<<- rbind(visited.space, solution[i,])
        }
      } else {
        visited.space<<- rbind(visited.space, solution)
      }

    },

    getVisitedSpace = function() {
      data.frame(visited.space, stringsAsFactors=FALSE)
    }

  )


)


##
## ----- Options classes
##


#' @title Options
#'
#' @description The base class for the options for the optimization methods.
#'
#' @field type The configuration type
#' @field container The object holding the configuration otions
#'
#' @importFrom methods new
#' @export Options
#' @exportClass Options
Options<- setRefClass("Options",
  fields = list(
    type = 'ANY',
    neighborhood = 'ANY',
    container = 'list'
  ),

  methods = list(
    initialize = function() {
      type<<- 'none'
      container<<- list(iterations=100, trace=FALSE)
    },

    ## Set/Get the neighborhood function
    neighborhoodFunction = function(f=NULL) {
      if(!is.null(f)) {
        neighborhood<<- f
      }
      neighborhood
    },

    ##setNeighborhoodF = function(v) {
    ##  neighborhood<<- v
    ##},

    ##getNeighborhoodF = function() {
    ##  neighborhood
    ##},

    setType = function(v) {
      type<<- v
    },

    getType = function() {
      type
    },

    setValue = function(k , v) {
      container[k]<<- v
    },

    getValue = function(k) {
      (container[k])[[1]]
    },

    Keys = function() {
      names(container)
    },

    toString = function() {
      sstring<- c()
      values<- unlist(container)
      for(k in names(values)) { sstring<- paste0(sstring,k,"=",values[k],sep='\t') }
      sstring
    }
  )
)

#' @title OptionsPSO
#'
#' @description Options for PSO method
#'
#' @importFrom methods new
#' @export OptionsPSO
#' @exportClass OptionsPSO
OptionsPSO<- setRefClass("OptionsPSO", contains = "Options",
  methods = list(
    initialize = function() {
      callSuper()
      setType("pso")
      #setValue("iterations", 1000)
      setValue("N",16)
      setValue("phi1",1.193)
      setValue("phi2",1.193)
      setValue("W",0.721)
      neighborhoodFunction(pso.neighborhood.KN)
    }
  )
)

#' @title OptionsSAA
#'
#' @description Options for SAA method
#'
#' @field temperature The temperature function
#'
#' @importFrom methods new
#' @export OptionsSAA
#' @exportClass OptionsSAA
OptionsSAA<- setRefClass("OptionsSAA", contains = "Options",
  fields = list(
    temperature = 'ANY'
  ),

  methods = list(
    initialize = function() {
      callSuper()
      setType("saa")
      setValue("t0", 1*10^4)
      setValue("t.min", 10^-5)
      setValue("L", 50)
      setValue("d", 0.05)
      neighborhoodFunction(saa.neighborhoodN)
      setTemperatureF(saa.tcte)
    },

    setTemperatureF = function(v) {
      temperature<<- v
    },

    getTemperatureF = function() {
      temperature
    }
  )
)

#' @title OptionsACOR
#'
#' @description Options for ACOR method
#'
#' @importFrom methods new
#' @export OptionsACOR
#' @exportClass OptionsACOR
OptionsACOR<- setRefClass("OptionsACOR", contains = "Options",
  methods = list(
    initialize = function() {
      callSuper()
      setType("acor")
      setValue("n.ants", 64)  ## The number of simulated ants
      setValue("k", 32)       ## The archive size
      setValue("q", 0.2)      ## Locality of the search process
      setValue("Xi", 0.5)     ## Equivalent to evaporation rate, higher Xi reduce convergence speed
    }
  )
)


#' @title OptionsEES1
#'
#' @description Options for EvoPER Evolutionary Stratety 1
#'
#'
#' @importFrom methods new
#' @export OptionsEES1
#' @exportClass OptionsEES1
OptionsEES1<- setRefClass("OptionsEES1", contains = "Options",
  fields = list(
  ),

  methods = list(
    initialize = function() {
      callSuper()
      setType("ees1")
      setValue("N", 10)             ## Solution size
      setValue("mu", 0.3)           ## Fitness preference strenght
      setValue("rho", 0.01)         ## Mutation probability
      setValue("kkappa", 0.2)       ## Selective pressure
      setValue("iterations", 50)    ## Total number of iterations
    }
  )
)


#' @title OptionsEES2
#'
#' @description Options for Serial Dilutions method
#'
#' @field dilutions The desired dilutions
#'
#' @importFrom methods new
#' @export OptionsEES2
#' @exportClass OptionsEES2
OptionsEES2<- setRefClass("OptionsEES2", contains = "Options",
  fields = list(
  ),

  methods = list(
    initialize = function() {
      callSuper()
      setType("ees2")
      setValue("N", 100)           ## Solution size
      setValue("rho", 0.05)        ## Solution size
      setValue("iterations", 10)   ## Total number of iterations
    }
  )
)


## ##################################################################
##
## --- Entry point function -
##
## ##################################################################


#' @title extremize
#'
#' @description Entry point for optimization functions
#'
#' @param type The optimization method (aco,pso,saa,sda)
#' @param objective An instance of ObjectiveFunction (or subclass) class \link{ObjectiveFunction}
#' @param options An apropiate instance from a sublclass of \link{Options} class
#'
#' @examples \dontrun{
#'  f<- PlainFunction$new(f0.rosenbrock2)
#'
#'  f$Parameter(name="x1",min=-100,max=100)
#'  f$Parameter(name="x2",min=-100,max=100)
#'
#'  extremize("pso", f)
#' }
#'
#' @export
extremize<- function(type, objective, options = NULL) {

  switch(type,
    pso={
      optimization.fun<- abm.pso
    },

    saa={
      optimization.fun<- abm.saa
    },

    acor={
      optimization.fun<- abm.acor
    },

    ees1={
      optimization.fun<- abm.ees1
    },

    ees2={
      optimization.fun<- abm.ees2
    },

    {
      stop("Invalid optimization function!")
    }
  )

  optimization.fun(objective, options)
}


## ##################################################################
##
## --- Particle Swarm Optimization functions -
##
## ##################################################################



#' @title abm.pso
#'
#' @description An implementaion of Particle Swarm Optimization
#' method for parameter estimation of Individual-based models.
#'
#' @param objective An instance of ObjectiveFunction (or subclass) class \link{ObjectiveFunction}
#' @param options An apropiate instance from a sublclass of \link{Options} class
#'
#' @examples \dontrun{
#'  f<- RepastFunction$new("c:/usr/models/BactoSim(HaldaneEngine-1.0)","ds::Output",300)
#'
#'  f$Parameter(name="cyclePoint",min=0,max=90)
#'  f$Parameter(name="conjugationCost",min=0,max=100)
#'  f$Parameter(name="pilusExpressionCost",min=0,max=100)
#'  f$Parameter(name="gamma0",min=1,max=10)
#'
#'  abm.pso(f)
#' }
#'
#' @references
#'
#' [1] Kennedy, J., & Eberhart, R. (1995). Particle swarm optimization.
#' In Proceedings of ICNN 95 - International Conference on Neural
#' Networks (Vol. 4, pp. 1942-1948). IEEE.
#'
#' [2] Poli, R., Kennedy, J., & Blackwell, T. (2007). Particle swarm optimization.
#' Swarm Intelligence, 1(1), 33-57.
#'
#'
#' @importFrom rrepast col.sum
#' @export
#abm.pso<- function(objective, iterations=100, N=16, phi1=1.193, phi2=1.193, W=0.721, f.neighborhood=pso.neighborhood.KN) {
abm.pso<- function(objective, options = NULL) {
  ## Handling the heuristic specific options
  options<- OptionsFactory("pso", options)

  ## --- Creating the estimation object for returning results
  estimates<- Estimates$new()

  iterations<- options$getValue("iterations")
  N<- options$getValue("N")
  phi1<- options$getValue("phi1")
  phi2<- options$getValue("phi2")
  W<- options$getValue("W")
  f.neighborhood<- options$neighborhoodFunction()


  elog.info("Initializing solution")
  Pg<- Pi<- x<- initSolution(objective$parameters,N)
  lbest<- pbest<- objective$EvalFitness(x)
  Pi<- x

  for(i in 1:nrow(x)){
    lbest[i,]<- pso.lbest(i,pbest,f.neighborhood)
    Pg[i,]<- Pi[lbest[i,"pset"],]
  }

  vi<- Pi * 0.1

  elog.info("Starting metaheuristic")
  for(index in 1:iterations) {
    #elog.info("Iteration=%d/%d, fitness=%g", index, iterations, (pso.best(lbest, Pg))["fitness"])

    vi<- pso.Velocity(W,vi,phi1,phi2,Pi,Pg,x)
    x<- enforceBounds((x + vi), objective$parameters)

    f1<- objective$EvalFitness(x)

    for(i in 1:nrow(x)){
      if(f1[i,"fitness"] < pbest[i,"fitness"]) {
        pbest[i,"fitness"]<- f1[i,"fitness"]
        Pi[i,]<- x[i,]
      }
    }

    for(i in 1:nrow(x)){
      gbest<- pso.lbest(i,pbest,f.neighborhood)
      if(gbest[1,"fitness"] < lbest[i,"fitness"]) {
        lbest[i,]<- gbest
        Pg[i,]<- Pi[lbest[i,"pset"],]
      }
    }

    ## --- Exit if convergence criteria is met
    if(objective$isConverged(gbest[1,"fitness"])) break

    s1<- sortSolution(x,f1)
    ## Show current bests
    elog.info("Iteration=%d/%d, best fitness=%g, iteration best fitness=%g", index, iterations, (pso.best(lbest, Pg))["fitness"], s1[1,"fitness"])

    estimates$addVisitedSpace(s1)

    ## --- Storing the best of this iteration
    estimates$addIterationBest(index, s1[1,]) # pso.best(lbest, Pg))
  }
  estimates$setBest(pso.best(lbest, Pg))
  estimates
}

#' @title pso.best
#'
#' @description Search for the best particle solution which minimize
#' the objective function.
#'
#' @param objective The results of evaluating the objective function
#' @param particles The particles tested
#'
#' @return The best particle
#'
#' @export
pso.best<- function(objective, particles) {
  best<- which(objective$fitness == min(objective$fitness),arr.ind = TRUE)
  best<- ifelse(length(best) > 1, best[1],best)
  particle<- merge(particles[best,],objective[best,c("pset","fitness")])
  return(particle)
}

#' @title pso.printbest
#'
#' @description Shows the best particle of each of simulated generations
#'
#' @param objective An instance of ObjectiveFunction (or subclass) class \link{ObjectiveFunction}
#' @param particles The current particle population
#' @param generation The current generation
#' @param title Some informational text to be shown
#'
pso.printbest<- function(objective, particles, generation, title) {
  particle<- pso.best(objective, particles)
  print(paste0("[", generation, "] ", title, " --------------------"))
  print(particle)
}

#' @title pso.velocity
#'
#' @description Calculates the PSO Velocity
#'
#' @param W Weight (Inertia weight or constriction coefficient)
#' @param Vi Current Velocity vector
#' @param phi1 Acceleration coefficient toward the previous best
#' @param phi2 Acceleration coefficient toward the global best
#' @param Pi Personal best
#' @param Pg Neighborhood best
#' @param Xi Particle vector
#'
#' @return Updated velocity
#'
#' @export
pso.Velocity<- function(W=1, Vi, phi1, phi2, Pi, Pg, Xi) {
  W * Vi + runif(length(Xi),0,phi1) * (Pi-Xi) + runif(length(Xi),0,phi2) * (Pg-Xi)
}

#' @title pso.chi
#'
#' @description Implementation of constriction coefficient
#'
#' @param phi1 Acceleration coefficient toward the previous best
#' @param phi2 Acceleration coefficient toward the global best
#'
#' @return The calculated constriction coefficient
#'
#' @export
pso.chi<- function(phi1, phi2) {
  phi<- phi1 + phi2
  assert(phi>4,"Constriction coefficient must be greater than zero")
  return(2/(phi - 2 + sqrt(phi^2 + 4 * phi)))
}

#' @title pso.neighborhood.K2
#'
#' @description The neighborhood function for a simple linear topology
#' where every particle has k = 2 neighbors
#'
#' @param i The particle position
#' @param n the size of particle population
#'
#' @export
pso.neighborhood.K2<- function(i,n) {
  assert((i >= 1 && i <= n),"Invalid particle position")
  i0<- ifelse((i-1) > 0, i-1, n)
  i1<- ifelse((i+1) > n, 1, i+1)
  return(c(i0,i,i1))
}

#' @title pso.neighborhood.K4
#'
#' @description The von neumann neighborhood function for
#' a lattice-based topology where every particle
#' has k = 4 neighbors
#'
#' @param i The particle position
#' @param n the size of particle population
#'
#' @export
pso.neighborhood.K4<- function(i,n) {
  assert((i >= 1 && i <= n),"Invalid particle position")
  m<- matrix(seq(1,n),nrow=(floor(sqrt(n))))
  p<- which(m == i, arr.ind = TRUE)
  p.row<- as.integer(p[1,"row"])
  p.col<- as.integer(p[1,"col"])

  x0<- ifelse((p.col-1) > 0, p.col-1, ncol(m))
  x1<- ifelse((p.col+1) > ncol(m), 1, p.col+1)
  y0<- ifelse((p.row-1) > 0, p.row-1, nrow(m))
  y1<- ifelse((p.row+1) > nrow(m), 1, p.row+1)

  return(c(m[p.row,x0],m[y0,p.col],i,m[p.row,x1],m[y1,p.col]))
}

#' @title pso.neighborhood.KN
#'
#' @description Simple helper method for 'gbest' neighborhood
#'
#' @param i The particle position
#' @param n the size of particle population
#'
#' @export
pso.neighborhood.KN<- function(i,n) {
  return(seq(1,n))
}

pso.neighborhood.KR<- function(i,n) {
  sample(setdiff(seq(1,n),c(i)),4)
}

#' @title pso.lbest
#'
#' @description Finds the lbest for the particle 'i' using the
#' topology function given by the topology parameter.
#'
#' @param i The particle position
#' @param pbest The pbest particle collection
#' @param topology The desired topology function
#'
#' @return The lbes for i th particle
#'
#' @export
pso.lbest<- function(i,pbest, topology) {
  neighbors<- topology(i,nrow(pbest))
  p<- which(pbest == min(pbest[neighbors,"fitness"]),arr.ind = TRUE)
  p<- as.integer(p[1,"row"])
  return(pbest[p,])
}


## ##################################################################
##
## --- Simulated Annealing Algorithm functions -
##
## ##################################################################


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
  s<- initSolution(objective$parameters,1)
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
      s1<- f.neighborhood(objective,getSolution(S),0.01)
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
    #newS[,i]<- newS[,i] + runif(1,as.numeric(f$getParameterValue(k,"min")),as.numeric(f$getParameterValue(k,"max"))) * distance
    newS[,i]<- newS[,i] + runif(1,-1,1) * distance
    #>>newS[,i]<- newS[,i] + 0.01 * f.range(k) * stats::rnorm(1)
    #newS[,i]<- newS[,i] + .01 * f.range(k) * runif(1,0,1)
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

## ##################################################################
##
## --- Ant colony optimization functions -
##
## ##################################################################

#' @title Ant colony optimization for continuous domains
#'
#' @description An implementation of Ant Colony Optimization algorithm
#' for continuous variables.
#'
#' @param objective An instance of ObjectiveFunction (or subclass) class \link{ObjectiveFunction}
#' @param options An apropiate instance from a sublclass of \link{Options} class
#'
#' @references
#'
#' [1] Socha, K., & Dorigo, M. (2008). Ant colony optimization for continuous domains.
#' European Journal of Operational Research, 185(3), 1155-1173.
#' http://doi.org/10.1016/j.ejor.2006.06.046
#'
#' @export
abm.acor<- function(objective, options= NULL) {
  ## Handling the heuristic specific options
  options<- OptionsFactory("acor", options)

  ## --- Creating the estimation object for returning results
  estimates<- Estimates$new()


  ## --- Configure algorithm parameters
  iterations<- options$getValue("iterations")
  n.ants<- options$getValue("n.ants")
  k<- options$getValue("k")
  q<- options$getValue("q")
  Xi<- options$getValue("Xi")

  ## --- Create the weigth vector
  W<- acor.weigth(q,k,1:k)

  ## --- Initialize the solution
  elog.info("Initializing solution")
  S<- initSolution(objective$parameters,n.ants)
  C<- objective$EvalFitness(S)
  T<- acor.archive(S, C, W, k)

  elog.info("Starting metaheuristic")
  for(index in 1:iterations) {
    s.s<- acor.S(T)
    s.sd<- acor.sigma(Xi, k, T)

    ## --- AntBasedSolutionConstruction
    S<- acor.updateants(S, n.ants, W, s.s, s.sd)
    C<- objective$EvalFitness(S)

    ## --- PheromoneUpdate
    T<- acor.archive(S, C, W, k, T)

    ## Show current bests
    elog.info("Iteration=%d/%d, best fitness=%g, iteration best fitness=%g", index, iterations, T[1, "fitness"], (sortSolution(S, C))[1,"fitness"])

    ## --- Storing the best of this iteration
    estimates$addIterationBest(index, sortSolution(S,C)[1,])

    estimates$addVisitedSpace(sortSolution(S,C))

    ## Check for algorithm convergence
    if(objective$isConverged(T[1, "fitness"])) break
  }

  estimates$setBest(T[1,])
  estimates
}

#' @title acor.updateants
#'
#' @description Update the solution using the gaussian kernel
#'
#' @param S The current solution ants
#' @param N The numnber of required ants in solution
#' @param W The weight vector
#' @param t.mu The 'mean' from solution archive
#' @param t.sigma The value of sigma from solution archive
#'
#' @return The new solution ants
#'
#' @references
#'
#' [1] Socha, K., & Dorigo, M. (2008). Ant colony optimization for continuous domains.
#' European Journal of Operational Research, 185(3), 1155-1173.
#' http://doi.org/10.1016/j.ejor.2006.06.046
#'
#' @export
acor.updateants<- function(S, N, W, t.mu, t.sigma) {
  t.mu<- as.matrix(t.mu)
  t.sigma<- as.matrix(t.sigma)
  n<- length(S[1,])
  for(i in 1:N) {
    l<- acor.lthgaussian(W)
    S[i,]<- stats::rnorm(n, t.mu[l,], t.sigma[l,])
  }
  S
}

#' @title Weight calculation for ant colony optimization
#'
#' @description Calculates the weight element of ACOr algorithm for the
#' solution archive.
#'
#' @param q The Algorithm parameter. When small best-ranked solution is preferred
#' @param k The Archive size
#' @param l The lth element of algorithm solution archive T
#'
#' @return A scalar or a vector with calculated weigth.
#'
#' @references
#'
#' [1] Socha, K., & Dorigo, M. (2008). Ant colony optimization for continuous domains.
#' European Journal of Operational Research, 185(3), 1155-1173.
#' http://doi.org/10.1016/j.ejor.2006.06.046
#'
#' @export
acor.weigth<- function(q, k, l) {
  ( 1 / (q * k * sqrt(2 * pi)) ) * exp(1) ^ - ( ( (l - 1) ^ 2 ) / ( 2 * q ^ 2 + k ^ 2) )
}

#' @title Gaussian kernel choosing probability
#'
#' @description Calculate the probability of choosing the lth Gaussian function
#'
#' @param W The vector of weights
#' @param l The lth element of algorithm solution archive T
#'
#' @return The vector of probabilities 'p'
#'
#' @references
#'
#' [1] Socha, K., & Dorigo, M. (2008). Ant colony optimization for continuous domains.
#' European Journal of Operational Research, 185(3), 1155-1173.
#' http://doi.org/10.1016/j.ejor.2006.06.046
#'
#' @export
acor.probabilities<- function(W, l= NULL) {
  if(is.null(l)) {
    l<- 1:length(W)
  }
  W[l]/sum(W)
}

#' @title Select the lth gaussian function
#'
#' @description Given a weight vector calculate the probabilities of selecting
#' the lth gaussian function and return the index of lht gaussian selected with
#' probability p
#'
#' @param W The vector of weights
#'
#' @return The index of lht gaussian function
#'
#' @references
#'
#' [1] Socha, K., & Dorigo, M. (2008). Ant colony optimization for continuous domains.
#' European Journal of Operational Research, 185(3), 1155-1173.
#' http://doi.org/10.1016/j.ejor.2006.06.046
#'
#' @export
acor.lthgaussian<- function(W) {
  which(runif(1) < cumsum(acor.probabilities(W)))[1]
}

#' @title Sigma calculation for ACOr
#'
#' @description Calculate the value of sigma
#'
#' @param Xi The algorithm parameter
#' @param k The solution archive size
#' @param T The solution archive
#'
#' @return The sigma value
#'
#' @references
#'
#' [1] Socha, K., & Dorigo, M. (2008). Ant colony optimization for continuous domains.
#' European Journal of Operational Research, 185(3), 1155-1173.
#' http://doi.org/10.1016/j.ejor.2006.06.046
#'
#' @export
acor.sigma<- function(Xi, k, T) {
  #assert((length(i) <= acor.N(T) && length(l) <= k),"Invalid i or l dimensions")
  s<- as.matrix(acor.S(T))
  ssigma<- matrix(nrow = k, ncol = acor.N(T))
  for(l in 1:k) {
    summatory<- 0
    for(e in 1:k) {
      summatory<- summatory + abs(s[e,] - s[l,]) / (k - 1)
    }
    ssigma[l,]<- Xi * summatory
  }
  ssigma
}

#' @title acor.archive
#'
#' @description This function is used for creating and maintaining the ACOr archive 'T'.
#' The function keeps the track of 'k' solotion in the archive.
#'
#' @param s The solution 'ants'
#' @param f The evaluation of solution
#' @param w The weight vector
#' @param k The archive size
#' @param T The current archive
#'
#' @return The solution archive
#'
#' @references
#'
#' [1] Socha, K., & Dorigo, M. (2008). Ant colony optimization for continuous domains.
#' European Journal of Operational Research, 185(3), 1155-1173.
#' http://doi.org/10.1016/j.ejor.2006.06.046
#'
#' @export
acor.archive<- function(s, f, w, k, T= NULL) {
  ## Drop "pset" column from f and leave only the fitness value
  f[,c("pset")]<- list(NULL)

  if(is.null(T)) {
    archive<- cbind(s,f)
  } else {
    t.s<- acor.S(T)
    t.f<- acor.F(T)
    archive<- cbind(t.s,t.f)
    archive<- rbind(archive,cbind(s,f))
  }
  T<- (archive[with(archive,order(fitness)),])[1:k,]
  cbind(T,w)
}

#' @title acor.S
#'
#' @description Helper function for extracting solution 'S' from archive 'T'
#'
#' @param T The solution archive
#'
#' @return The solution matrix
#'
#' @export
acor.S<- function(T) {
  T[,c("fitness","w")]<- list(NULL)
  T
}

#' @title acor.F
#'
#' @description Helper function for extracting the 'F' function evaluations
#' from archive ACOr 'T'
#'
#' @param T The solution archive
#'
#' @return The F matrix
#'
#' @export
acor.F<- function(T) {
  T[,c(names(acor.S(T)), "w")]<- list(NULL)
  T
}

#' @title acor.W
#'
#' @description Helper function for extracting the 'W' function evaluations
#' from archive ACOr 'T'
#'
#' @param T The solution archive
#'
#' @return The weight vector
#'
#' @export
acor.W<- function(T) {
  T[,c(names(acor.S(T)), "fitness")]<- list(NULL)
  T
}

#' @title acor.N
#'
#' @description Helper function for getting the size of solution
#'
#' @param T The solution archive
#'
#' @return The size 'n' of a solution 's'
#'
#' @export
acor.N<- function(T) {
  length(T[1,])-2
}


## ##################################################################
##
## --- EvoPER Evolutionary Strategy 1
##
## ##################################################################

#' @title EvoPER Evolutionary Strategy 1
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
abm.ees1<- function(objective, options= NULL) {
  ## Handling the heuristic specific options
  options<- OptionsFactory("ees1", options)
  elog.info("Options(%s): %s", options$getType(), options$toString())

  ## --- Creating the estimation object for returning results
  estimates<- Estimates$new()

  ## --- Metaheuristic options
  k<- objective$ParametersSize()
  iterations<- options$getValue("iterations")
  N<- options$getValue("N")
  mu<- options$getValue("mu")                 ##
  rho<- options$getValue("rho")               ##
  kkappa<- options$getValue("kkappa")         ##

  #iterations<- 100
  #N<- 6

  elog.info("Initializing solution")
  s<- initSolution(objective$parameters, N, "lhs")
  #s<- ees1.challenge(s,2)

  ## -- Evaluate
  s0<- es.evaluate(objective, s)

  elog.info("Starting metaheuristic")
  for(iteration in 1:iterations) {
    mates<- ees1.mating1(getSolution(s0), mu)
    s<- ees1.recombination(s0, mates)
    s<- ees1.mutation(s, mates, rho)

    ## --- Evaluate solution
    s1<- es.evaluate(objective, getSolution(s))

    ## --- Selection
    s0<- ees1.selection(s0, s1, kkappa)


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

#' @title es.evaluate
#'
#' @description For each element in solution 's' evaluate the respective
#' fitness.
#'
#' @param f A reference to an instance of objective function
#' @param s The set of solutions
#' @param enforce If true the values are enforced to fall within provided range
#'
#' @return The solution ordered by its fitness.
#' @export
es.evaluate<- function(f, s, enforce=TRUE) {
  if(enforce) {
    s<- enforceBounds(s, f$parameters)
  }
  c<- f$EvalFitness(s)
  sortSolution(s, c)
}

#' @title ees1.mating
#'
#' @description This function 'mix' the elements present in the solution. The
#' parameter 'mu' controls the intensity of mixing. Low values give preference
#' to best solution components and high values make the values being select randomly.
#'
#' @param solution The Problem solution
#' @param mu The mixing intensity ratio, from 0 to 1. The mix intensity
#' controls de the probability of chosing a worst solutions
#'
#' @export
ees1.mating<- function(solution, mu) {
  P<- function(p, x) { (p^x) }
  k<- length(solution[,1])
  p<- c( P(mu,1:k ) )

  sampling<- c()
  indexes<- 1:k

  ## Repeat until sampling have the right size
  while(length(sampling) < (k/2)) {
    for(i in 1:length(indexes)) {
      if(runif(1) < (1/k * p[i])) {
        sampling<- c(sampling, indexes[i])
        indexes<- indexes[-c(i)]
        p<- p[-c(i)]
        break
      }
    }
  }
  solution[sampling,]
}

#' @title ees1.mating1
#'
#' @description This function 'mix' the elements present in the solution. The
#' parameter 'mu' controls the intensity of mixing. Low values give preference
#' to best solution components and high values make the values being select randomly.
#'
#' @param solution The Problem solution
#' @param mu The mixing intensity ratio, from 0 to 1. The mix intensity
#' controls de the probability of chosing a worst solutions
#'
#' @export
ees1.mating1<- function(solution, mu) {
  P<- function(p, x) { (p^x) }
  k<- length(solution[,1])
  #p<- c( rep(1,(k/2)), P(mu,((k/2)+1):k) )
  p<- c( P(mu,1:k ) )

  solution[sample(1:k, size = (k/2), prob = p),]
}

#' @title ees1.recombination
#'
#' @description Performs the recombination on solution
#'
#' @param solution The Problem solution
#' @param mates The mixed parents
#'
#' @export
ees1.recombination<- function(solution, mates) {
  m<- length(solution[,1])
  n<- length(solution[1,])

  solute<- apply(getSolution(mates),2,gm.mean)

  summatory<- sum(solution[,"fitness"])

  ssi<- getSolution(solution)
  stock<- c()
  for(i in 1:m) {
    si<- ssi[i,]

    fitness<- solution[i,"fitness"]
    weight<- fitness/summatory

    if(runif(1) < 1/5) {
      #stock<- rbind(stock, apply(rbind(si + si * runif(n,-1, 1), solute * weight),2,mean) )
      stock<- rbind(stock, apply(rbind(si, (si + solute) * weight),2,mean) )
      #stock<- rbind(stock, apply(rbind(si, si +solute),2,mean) )
    } else {
      stock<- rbind(stock, apply(rbind(si, solute),2,mean) )
      #stock<- rbind(stock,ees1.explore(si,fitness))
    }
  }
  as.data.frame(cbind(stock,fitness=getFitness(solution)))
}

#' @title ees1.mutation
#'
#' @description Performs the mutation on generated solution
#'
#' @param solution The Problem solution
#' @param mates The mixed parents
#' @param p The mutation probability
#'
#' @export
ees1.mutation<- function(solution, mates, p=0.01) {
  m<- length(solution[,1])
  n<- length(solution[1,])
  s<- getSolution(solution)
  for(i in 1:m) {
    s[i,]<- ees1.explore(s[i,],getFitness(solution, i))
  }
  s
}

#' @title ees1.explore
#'
#' @description Explore the solution space on the neighborhood of
#' solution 's' in order to find a new best.
#'
#' @param s The Problem solution
#' @param weight The exploration intensity
#' @param p The mutation probability
#'
#' @export
ees1.explore<- function(s, weight, p=0.01) {
  m<- length(s[,1])
  n<- length(s[1,])
  w<- Magnitude(weight)
  w<- ifelse(w < 0, 1/abs(w), 1)

  for(i in 1:m) {
    for(j in 1:n) {
      if(runif(1) < p) {
        mm<- Magnitude(s[i,j])
        e<- ifelse(mm < 0, 0, mm)
        p<- s[i,j]
        s[i,j]<- s[i,j] + runif(1,-10^e,10^e) * w
      }
    }
  }
  s
}

#' @title ees1.challenge
#'
#' @description Repeat the evalution of best solution to tacke with
#' variability.
#'
#' @param solution The Problem solution
#' @param objective The objective function
#'
#' @export
ees1.challenge<- function(solution, objective) {
  s<- getSolution(solution)
  ss0<- s + .Machine$double.eps^0.5
  solution1<- es.evaluate(objective, ss0)
  s0<- solution[with(solution,order(pset)),]
  s1<- solution1[with(solution1,order(pset)),]
  for(i in 1:length(s0[,1])) {
    if(getFitness(s1, i) > getFitness(s0, i)) {
      s0[i,]<- s1[i,]
    }
  }
  s0
}

#' @title ees.selection
#'
#' @description Select the elements with best fitness but accept uphill
#' moves with probability 'kkappa'.
#'
#' @param s0 The current best solution set
#' @param s1 The new solution
#' @param kkappa The selection pressure
#'
#' @export
ees1.selection<- function(s0, s1, kkappa) {
  assert(length(s0[1,]) == length(s1[1,]),"Invalid solution!")
  k<- length(getSolution(s0)[1,]) - 1
  m<- ifelse(length(s0[,1]) < length(s1[,1]), length(s0[,1]), length(s1[,1]))

  for(i in 1:m) {
    for(j in i:m) {
      if(s1[j,"fitness"] < s0[i,"fitness"]) {
        s0[i,]<- s1[j,]
        break
      } else {
        if(i > k && runif(1) < kkappa) {
          s0[i,]<- s1[j,]
        }
      }
    }
  }
  s0
}

## ##################################################################
##
## --- EvoPER Evolutionary Strategy 2
##
## ##################################################################

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
    mmean<- apply(s,2,mean)
    ssd<- apply(s,2,sd)
    interval<-  abs(mmax-mmin)/2

    for(k in parameters$name) {
      parameters[which(parameters$name == k),"min"]<- as.numeric((mmean-interval-runif(1))[k])
      parameters[which(parameters$name == k),"max"]<- as.numeric((mmean+interval+runif(1))[k])
    }



    s<- initSolution(parameters, N, "lhs")

    ## -- Evaluate
    s1<- es.evaluate(objective, s)

    ## -- Save the best iteration value
    estimates$addIterationBest(iteration, s1[iteration,])

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


## ##################################################################
##
## ---------------------- General functions -------------------------
##
## ##################################################################


#' @title initSolution
#'
#' @description Creates the initial Solution population
#' taking into account the lower an upper bounds of
#' provided experiment factors.
#'
#' @param parameters The Objective Function parameter list
#' @param N The size of Solution population
#' @param sampling The population sampling scheme (mcs|lhs|ffs)
#'
#' @return A random set of solutions
#'
#' @importFrom rrepast AoE.RandomSampling AoE.LatinHypercube
#' @export
initSolution<- function(parameters, N=20, sampling="mcs") {
  switch(sampling,
    lhs={
      rrepast::AoE.LatinHypercube(N,parameters)
    },

    mcs={
      rrepast::AoE.RandomSampling(N,parameters)
    },

    ffs={
      rrepast::AoE.FullFactorial(N,parameters)
    }
  )
}

#' @title partSolutionSpace
#'
#' @description Creates the initial Solution population
#' taking into account the lower an upper bounds of
#' provided experiment factors. This method works by
#' dividing the solution space into partitions of size 'd'
#' and then creating a full factorial combination of partitions.
#'
#' @param parameters The Objective Function parameter list
#' @param d The partition size. Default value 4.
#'
#' @return A set of solutions
#'
#' @export
partSolutionSpace<- function(parameters, d=4) {
  i<- 1
  m<- list()
  keys<- c()

  for(k in 1:length(parameters[,1])) {
    keys<- c(keys,p<- parameters[k, "name"])
    p<- parameters[k, ]
    mmin<- as.numeric(p$min)
    mmax<- as.numeric(p$max)
    step<- ceiling((mmax-mmin) / d)
    print(sprintf("a=%g, b=%g, c=%g", mmin+step/(d/2), mmax, ceiling((mmax-mmin) / d) ))
    s<- seq(mmin+step/(d/2), mmax, ceiling((mmax-mmin) / d) )
    m[i] <- list(round(s))
    i<- i + 1
  }
  g<- expand.grid(m)
  names(g)<-  keys
  g
}

#' @title sortSolution
#'
#' @description Sort solution by its respective fitness
#'
#' @param s Problem solution
#' @param f The function evaluation for s
#'
#' @export
sortSolution<- function(s, f) {
  solution<- cbind(s,f[,c("pset","fitness")])
  as.data.frame(solution[with(solution,order(fitness)),])
}

#' @title bestFitness
#'
#' @description Given a set S of N solutions created with sortSolution, this function
#' returns the fitness component fot the best solution.
#'
#' @param S The solution set
#'
#' @return The best fitness value
#'
#' @export
bestFitness<- function(S) {
  S[1,"fitness"]
}

#' @title getFitness
#'
#' @description Given a set S of N solutions created with sortSolution, this function
#' returns the solution component fot the best solution.
#'
#' @param S The solution set
#' @param i The fitness index, if null return the whole column.
#'
#' @return The selected fitness entry
#'
#' @export
getFitness<- function(S, i=NULL) {
  v<- S[,"fitness"]
  if(!is.null(i)) {
    v<- S[i,"fitness"]
  }
  v
}

#' @title bestSolution
#'
#' @description Given a set S of N solutions created with sortSolution, this function
#' returns the best solution found.
#'
#' @param S The solution set
#'
#' @return The best solution
#'
#' @export
bestSolution<- function(S) {
  S[,c("pset","fitness")]<- list(NULL)
  S[1,]
}

#' @title getSolution
#'
#' @description Given a set S of N solutions created with sortSolution, this function
#' returns the solution component. A solutions is a set of solutions and their associated
#' fitness
#'
#' @param S The solution set
#'
#' @return The solution set
#'
#' @export
getSolution<- function(S) {
  S[,c("pset","fitness")]<- list(NULL)
  S
}

#' @title lowerBound
#'
#' @description Checks if parameters is greater than the lower bounds
#' @param particles The particle set
#' @param factors the defined range for objective function parameters
#'
#' @return The particle greater than or equal to lower limit
#'
#'
#' @export
lowerBound<- function(particles, factors) {
  k<- rrepast::GetFactorsSize(factors)
  for(p in 1:nrow(particles)){
    for(i in 1:k) {
      lb<- as.numeric(factors[i,"min"])
      if(particles[p,i] < lb) { particles[p,i]<- lb}
    }
  }
  return(particles)
}

#' @title upperBound
#'
#' @description Checks if parameters is below the upper bounds
#' @param particles The particle set
#' @param factors the defined range for objective function parameters
#'
#' @return The particle inside the valid upper bound
#'
#' @export
upperBound<- function(particles, factors) {
  k<- rrepast::GetFactorsSize(factors)
  for(p in 1:nrow(particles)){
    for(i in 1:k) {
      ub<- as.numeric(factors[i,"max"])
      if(particles[p,i] > ub) { particles[p,i]<- ub}
    }
  }
  return(particles)
}

#' @title enforceBounds
#'
#' @description Checks if parameters fall within upper an lower bounds
#'
#' @param particles The particle set
#' @param factors the defined range for objective function parameters
#'
#' @return The particle inside the valid limits
#'
#' @export
enforceBounds<- function(particles, factors) {
  k<- rrepast::GetFactorsSize(factors)
  for(p in 1:nrow(particles)){
    for(i in 1:k) {
      lb<- as.numeric(factors[i,"min"]);
      ub<- as.numeric(factors[i,"max"])
      if( particles[p,i] < lb || particles[p,i] > ub || is.na(particles[p,i])) {
        particles[p,i]<- runif(1,lb,ub)
      }
    }
  }
  return(particles)
}

#' @title cbuf
#'
#' @description Simple implementation of a circular buffer.
#'
#' @param b The variable holding the current buffer content
#' @param v The new valued to be added to b
#' @param e The length of circular buffer
#'
#' @return The buffer b plus the element v minus the least recently added element
#'
#' @importFrom utils head
#' @export
cbuf<- function(b,v,e) {
  c(head(v,1),head(b,e))
}

#' @title slopes
#'
#' @description Calcule all slopes for the discrete x,y series
#'
#' @param x The x vector
#' @param y The y vector
#'
#' @return A vector with all slopes
#'
#' @export
slopes<- function(x, y) {
  sapply(1:length(x),function(i, x , y) { slope(x,y,i) }, x=x, y=y)
}

#' @title slope
#'
#' @description Simple function for calculate the slope on the ith element
#' position
#'
#' @param x The x vector
#' @param y The y vector
#' @param i The position
#'
#' @return The slope
#'
#' @export
slope<- function(x, y, i) {
  if(i > 1 && i < length(y)) {
    v<- 0.5 * ((y[i]-y[i-1])/(x[i]-x[i-1]) + (y[i+1] - y[i])/(x[i+1]-x[i]))
  } else {
    v<- ifelse(i == 1,(y[i+1] - y[i])/(x[i+1]-x[i]),(y[i]-y[i-1])/(x[i]-x[i-1]))
  }
  v
}

#' @title gm.mean
#'
#' @description Simple implementation for geometric mean
#'
#' @param x data
#' @return geometric mean for data
#'
#' @export
gm.mean<- function(x) {
  x<- as.vector(x)
  k<- length(x)
  P<- x[1]
  for(i in 2:k) {
    P<- P * x[i]
  }
  (abs(P)^(1/k)) * sign(P)
}

#' @title gm.sd
#'
#' @description Simple implementation for geometric standard deviation
#'
#' @param x data
#' @param mu The geometric mean. If not provided it is calculated.
#'
#' @return geometric standard deviation for data
#'
#' @export
gm.sd<- function(x, mu= NULL) {
  if(is.null(mu)) {
    mu<- gm.mean(x)
  }
  ssum<- 0
  k<- length(x)
  for(i in 1:k) {
    ssum<- ssum + (x[i]-mu)^2
  }
  exp(sqrt(ssum/k))
}

#' @title naiveperiod
#'
#' @description A naive approach for finding the period in a series
#' of data points
#'
#' @param d The data to search period
#'
#' @return A list with the average period and amplitude
#' @export
naiveperiod<- function(d) {
  s0<- 0
  i<- 0
  sum.max<- 0
  sum.per<- 0
  n<- 0
  distance<- 0

  for(i in 1:(length(d)-1)) {
    ## Checking data consistency
    if(is.na(d[i]) || is.na(d[i+1])) { break }

    ## Evaluate if current point is growing or decreasing
    if(d[i] < d[i+1]) { signal<- 1 }
    if(d[i] > d[i+1]) { signal<- -1 }
    if(d[i] == d[i+1]) { signal<- 0 }

    ## Increment distance between peaks
    distance<- distance + 1

    ## Detecting inflection points
    if(s0 != signal) {

      if(signal < s0) {
        ## Discard the first peak
        if(i>1) {
          sum.per= sum.per + distance
          n<- n + 1
          sum.max<- sum.max + d[i]
        }
        elog.debug("iteration=%d, amplitude=%g, period=%g",i, d[i], distance)
        distance<- 0
        i<- i + 1
      }
      s0<- signal
    }
  }
  elog.debug("summatory=%g, n=%g",sum.per, n)
  v<- list(period= (sum.per/n), value=(sum.max/n))
  v
}

#' @title OptionsFactory
#'
#' @description Instantiate the Options class required for
#' the specific metaheuristic method.
#'
#' @param type The metaheuristic method
#' @param v The options object
#'
#' @return Options object
#'
#' @export
OptionsFactory<- function(type, v=NULL) {
  if(is.null(v)) {
    switch(type,
      pso = { v<- OptionsPSO$new() },
      saa = { v<- OptionsSAA$new() },
      acor= { v<- OptionsACOR$new() },
      ees1= { v<- OptionsEES1$new() },
      ees2= { v<- OptionsEES2$new() },
      { stop("Invalid optimization function!") }
    )
  }

  if(v$getType() != type) stop(paste("Invalid option of type [", v$getType(),"]"))
  v
}

#' @title Magnitude
#'
#' @description Calculates the magnitude order for a given value
#'
#' @param v The numerical value
#'
#' @return The magnitude order
#'
#' @export
Magnitude<- function(v) {
  if(v == 0) {
    m<- 0
  } else {
    m<- trunc(log(abs(v),10))
  }
  ## Maybe a R bug
  m<- ifelse(m == -0, 0, m)
  m
}

#' @title xmeanci1
#'
#' @description Calculates confidence interval of mean for provided
#' data with desired confidence level. This functions uses bootstrap
#' resampling scheme for estimanting the CI.
#'
#' @param x The data set for which CI will be calculated
#' @param alpha The confidence level. The default value is 0.95 (95\%)
#'
#' @return The confidence interval for the mean calculated using 'boot.ci'
#'
#' @importFrom boot boot boot.ci
#' @importFrom stats sd
#' @export
xmeanci1<- function(x, alpha=0.95) {
  f<- function(x,i) { c( mean(x[i]), (sd(x[i])/sqrt(length(x[i])))^2 ) }

  b<- boot(x, f, R = 1000)
  boot.ci(b, conf = alpha, type = c("norm", "basic", "perc", "stud", "bca"))
}

#' @title xmeanci2
#'
#' @description Calculates confidence interval of mean for provided
#' data with desired confidence level.
#'
#' @param x The data set for which CI will be calculated
#' @param alpha The confidence level. The default value is 0.95 (95\%)
#'
#' @return The confidence interval for the mean
#' @importFrom stats sd qt
#' @export
xmeanci2<- function(x, alpha=0.95) {
  n<- length(x)
  avg<- mean(x)
  se<- sd(x)/sqrt(n);
  e<- se * qt(1-((1-alpha)/2), df=(n-1))
  c(avg-e, avg+e)
}


#' @title fixdfcolumns
#'
#' @description Coerce dataframe columns to a specic type.
#'
#' @param df The data frame.
#' @param cols The dataframe columns to be skiped or included.
#' @param skip If TRUE the column names in 'cols' are skiped. When FALSE logic is inverted.
#' @param type The type for which data frame columns must be converted.
#'
#' @return The data frame with converted column types.
#'
#' @export
fixdfcolumns<- function(df, cols= c(), skip=TRUE, type=as.numeric) {
  df<- as.data.frame(df)
  assert(is.data.frame(df),"Invalid data")
  df<- df[,!(names(df) %in% c("V1"))]

  for(k in names(df)) {
    if(skip) {
      if(k %in% cols) next
    } else {
      if(!k %in% cols) next
    }
    df[,k]<- type(df[,k])
  }
  df
}


#' @title xyplothelper
#'
#' @description Simple helper for ploting xy dispersion points.
#'
#' @param d A data frame.
#' @param x A string with the dataframe column name for x axis
#' @param y A string with the dataframe column name for y axis
#' @param title The plot title
#'
#' @return A ggplot2 plot object
#'
#' @importFrom ggplot2 aes_string labs geom_point ggtitle theme element_text
#' @export
xyplothelper<- function(d, x, y, title=NULL) {
  d<- as.data.frame(d)
  d<- fixdfcolumns(d,cols = c(x,y),skip = FALSE)

  p<- ggplot(d, with(d,aes_string(x = x, y = y)))
  p<- p + labs(y = y)
  p<- p + labs(x = x)
  p<- p + geom_point(size = 1.5)
  if(!is.null(title)) {
    p<- p + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  }
  p
}

#' @title histplothelper
#'
#' @description Simple helper for ploting histograms
#'
#' @param d A data frame.
#' @param x A string with the dataframe column name for histogram
#' @param title The plot title
#'
#' @return A ggplot2 plot object
#' @importFrom ggplot2 ggplot aes_string geom_histogram ggtitle theme element_text
#' @export
histplothelper<- function(d, x, title=NULL) {
  d<- as.data.frame(d)
  d<- fixdfcolumns(d,cols = c(x),skip = FALSE)


  p<- ggplot(d, aes_string(x=x))
  p<- p + geom_histogram(aes_string(y="..ncount.."), bins = 4, colour="black", fill="orange")
  if(!is.null(title)) {
    p<- p + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  }
  p
}

#' @title scatterplotlothelper
#'
#' @description Simple helper for ploting 3d scaterplots
#'
#' @param d A data frame.
#' @param x A string with the dataframe column name for x axis
#' @param y A string with the dataframe column name for y axis
#' @param z A string with the dataframe column name for z axis
#' @param title The plot title
#'
#' @return A scatter3D plot
#' @import plot3D
#' @export
scatterplotlothelper<- function(d, x, y, z, title=NULL) {
  d<- as.data.frame(d)
  d<- fixdfcolumns(d,cols = c(x, y, z),skip = FALSE)
  p<- scatter3D(x=d[,x],y=d[,y],z=d[,z],main=title,phi= 0,  bty = "g", pch = 20, cex = 1.5, ticktype = "simple", xlab=x, ylab=y,zlab=z)
  p
}


#' @title contourplothelper
#'
#' @description Simple helper for ploting histograms
#'
#' @param d A data frame.
#' @param x A string with the dataframe column name for x axis
#' @param y A string with the dataframe column name for y axis
#' @param z A string with the dataframe column name for z axis
#'
#' @importFrom ggplot2 ggplot aes_string geom_density2d stat_density2d stat_contour
#' @importFrom ggplot2 ggtitle theme element_text scale_fill_gradient facet_wrap stat_summary2d
#' @importFrom ggplot2 stat_summary_2d scale_alpha facet_grid xlim ylim
#' @export
contourplothelper<- function(d, x, y, z) {
  d<- as.data.frame(d)
  d<- fixdfcolumns(d,cols = c(x, y, z),skip = FALSE)

  ## --- Pruning duplicated entries
  d<- d[!duplicated(d[,c(x,y)]),]
  d<- d[!duplicated(d[,c(z)]),]

  d$title<- sprintf("solution landscape for %s ~ %s, %s", z, x, y)
  p<- with(d,ggplot(d, aes_string(x=x, y=y, z=z)))
  p<- p + stat_density2d(data = d, aes_string(fill = "..level..", alpha = "..level.."), size = 0.01, bins = 8, geom = 'polygon', n = 300)
  p<- p + scale_fill_gradient(low = "green", high = "red")
  p<- p + scale_alpha(range = c(0.15, 0.30), guide = FALSE)
  p<- p + facet_grid(. ~ title) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  p
}

## ##################################################################
##
## ----------------- Functions for testing package -----------------
##
## ##################################################################


##
# rm(list=ls())
# f<- PlainFunction$new(f0.rosenbrock2)
# f$Parameter(name="x1",min=-100,max=100)
# f$Parameter(name="x2",min=-100,max=100)
# extremize("acor", f)
##


#' @title f0.test
#'
#' @description Simple test function f(1,2,3,4) = 0
#'
#' @param x1 Parameter 1
#' @param x2 Parameter 2
#' @param x3 Parameter 3
#' @param x4 Parameter 4
#'
#' @export
f0.test<- function(x1,x2,x3,x4) { 10 * (x1 - 1)^2 + 20 * (x2 - 2)^2 + 30 * (x3 - 3)^2 + 40 * (x4 - 4)^2 }

#' @title f1.test
#'
#' @description Simple test function f(c(1,2,3,4)) = 0
#'
#' @param x Parameter vector
#'
#' @export
f1.test<- function(x) { f0.test(x[1],x[2],x[3],x[4]) }

#' @title f0.rosenbrock2
#'
#' @description Two variable Rosenbrock function, where f(1,1) = 0
#'
#' @param x1 Parameter 1
#' @param x2 Parameter 2
#'
#' @export
f0.rosenbrock2<- function(x1, x2) { (1 - x1)^2 + 100 * (x2 - x1^2)^2 }

#' @title f1.rosenbrock2
#'
#' @description Two variable Rosenbrock function, where f(c(1,1)) = 0
#'
#' @param x Parameter vector
#'
#' @export
f1.rosenbrock2<- function(x) { f0.rosenbrock2(x[1], x[2]) }

#' @title f0.adtn.rosenbrock2
#'
#' @description Two variable Rosenbrock function with random additive noise.
#'
#' @param x1 Parameter 1
#' @param x2 Parameter 2
#'
#' @export
f0.adtn.rosenbrock2<- function(x1, x2) { (1 - x1)^2 + 100 * (x2 - x1^2)^2 + runif(1,0,0.5)}

#' @title f1.adtn.rosenbrock2
#'
#' @description Two variable Rosenbrock function with random additive noise.
#'
#' @param x Parameter vector
#'
#' @export
f1.adtn.rosenbrock2<- function(x) { f0.adtn.rosenbrock2(x[1], x[2]) }


#' @title f0.nlnn.rosenbrock2
#'
#' @description Two variable Rosenbrock function with random additive noise.
#'
#' @param x1 Parameter 1
#' @param x2 Parameter 2
#'
#' @export
f0.nlnn.rosenbrock2<- function(x1, x2) {
  v<- (1 - x1)^2 + 100 * (x2 - x1^2)^2 + runif(1,0,0.5)
  if(runif(1) < 0.5) {
    v<- v * stats::rnorm(1,1,0.1)
  } else {
    v<- v + stats::rnorm(1,v,v * 0.1)
  }
  v
}

#' @title f1.nlnn.rosenbrock2
#'
#' @description Two variable Rosenbrock function with random additive noise.
#'
#' @param x Parameter vector
#'
#' @export
f1.nlnn.rosenbrock2<- function(x) { f0.nlnn.rosenbrock2(x[1], x[2]) }

#' @title f0.cigar
#'
#' @description The Cigar function of N variables for testing optimization
#' methods. The global optima for the function is given by
#' xi = 0, forall i E {1...N}, f(x) = 0.
#'
#' @param ... The variadic list of function variables.
#'
#' @return The function value
#'
#' @references
#'
#' http://deap.gel.ulaval.ca/doc/dev/api/benchmarks.html
#'
#' @export
f0.cigar<- function(...) {
  x<- list(...)
  f1.cigar(unlist(x))
}

#' @title f1.cigar
#'
#' @description The Cigar function of N variables for testing optimization
#' methods. The global optima for the function is given by
#' xi = 0, forall i E {1...N}, f(x) = 0.
#'
#' @param x The vector of function variables.
#'
#' @return The function value
#'
#' @references
#'
#' http://deap.gel.ulaval.ca/doc/dev/api/benchmarks.html
#'
#' @export
f1.cigar<- function(x) {
  ssum<- function(x) {
    s<- 0
    for(i in 2:length(x)) {
      s<- s + x[i]^2
    }
    s
  }
  x[1]^2 + 10^6 * ssum(x)
}

#' @title f0.cigar4
#'
#' @description The Cigar function of four variables for testing optimization
#' methods. The global optima for the function is given by
#' xi = 0, forall i E {1...N}, f(x) = 0.
#'
#' @param x1 The first function variable
#' @param x2 The second function variable
#' @param x3 The third function variable
#' @param x4 The fourth function variable
#'
#' @return The function value
#'
#' @export
f0.cigar4<- function(x1, x2, x3, x4) {
  f1.cigar(c(x1, x2, x3, x4))
}

#' @title predatorprey
#'
#' @description The solver for Lotka-Volterra differential equation.
#'
#' @param x1 The growth rate of prey
#' @param x2 The decay rate of predator
#' @param x3 The predating effect on prey
#' @param x4 The predating effecto on predator
#'
#' @return The ODE solution
#'
#' @export
predatorprey<- function(x1, x2, x3, x4) {
  value<- c(x = 12, y = 12)
  parameters<- c(c1 = x1, c2 = x2, c3 = x3, c4 = x4)
  time<- seq(0, 96, by = 1)

  ## Predator-Prey ODE function
  f.diffequation<- function (t, y, parms) {
    with(as.list(c(y, parms)), {
      dx = c1 * x - c3 * x * y
      dy = -c2 * y + c4 * x * y
      return(list(c(dx, dy)))
    })
  }

  deSolve::ode(func = f.diffequation, y = value, parms = parameters, times = time,  method = "radau")
}

#' @title predatorprey.plot
#'
#' @description Generate a plot for the predator-prey ODE output.
#'
#' @param x1 The growth rate of prey
#' @param x2 The decay rate of predator
#' @param x3 The predating effect on prey
#' @param x4 The predating effecto on predator
#'
#' @return An ggplot2 object
#'
#' @examples \dontrun{
#'  predatorprey.plot(1.351888, 1.439185, 1.337083, 0.9079049)
#' }
#'
#' @importFrom reshape melt
#' @importFrom ggplot2 ggplot aes geom_line ggtitle
#' @export
predatorprey.plot<- function(x1, x2, x3, x4) {
  v<- as.data.frame(predatorprey(x1, x2, x3, x4))
  v.data<- melt(v, id.vars="time", value.name="value", variable_name="species")
  p<- ggplot(data= v.data, with(v.data,aes(x=time, y=value, group = species, colour = species)))
  p + geom_line() + ggtitle("Predator/Prey period")
}

#' @title Period tuning for Predator-Prey
#'
#' @description This function is an example on how EvoPER can be
#' used for estimating the parameter values in order to produce
#' oscilations with the desired period.
#'
#' @param x1 The growth rate of prey
#' @param x2 The decay rate of predator
#' @param x3 The predating effect on prey
#' @param x4 The predating effecto on predator
#'
#' @return The solution fitness cost
#'
#' @examples \dontrun{
#'	rm(list=ls())
#'	set.seed(-27262565)
#'	f<- PlainFunction$new(f0.periodtuningpp)
#'	f$Parameter(name="x1",min=0.5,max=2)
#'	f$Parameter(name="x2",min=0.5,max=2)
#'	f$Parameter(name="x3",min=0.5,max=2)
#'	f$Parameter(name="x4",min=0.5,max=2)
#'	extremize("pso", f)
#' }
#'
#' @export
f0.periodtuningpp<- function(x1, x2, x3, x4) {
  v<- predatorprey(x1, x2, x3, x4)
  rrepast::AoE.NRMSD(naiveperiod(v[,"y"])$period,24)
}
