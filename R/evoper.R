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
    parameters = 'ANY',
    value = 'ANY',
    tolerance = 'ANY',
    converged = 'ANY',
    maximize = 'ANY',
    counter = 'ANY'),

  methods = list(
    initialize = function(funct) {
      object<<- NULL
      objective<<- funct
      parameters<<- NULL
      value<<- NULL
      tolerance<<- .Machine$double.eps^0.30
      converged<<- FALSE
      maximize<<- FALSE
      counter<<- 0
    },

    stats = function() {
      cbind(total_evals=counter,converged=converged)
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
      parameters[which(parameters[,"name"] == key),]
    },

    getParameterNames = function() {
      parameters[,"name"]
    },

    getParameterValue = function(key, name) {
      parameters[which(parameters[,"name"] == key),name]
    },

    Parameters = function() {
      return(parameters)
    },

    ParametersSize = function() {
      length(parameters[,1])
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
      object<<- RunExperiment(model,r=1,p, objective)

      # Sum the objective output and change the column name
      object$output<<- col.sum(object$output)
      n<- names(object$output)
      names(object$output)<<- replace(n, which(n == "total"),c("fitness"))

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
    iteration.best = 'ANY',
    overall.best = 'ANY'
  ),

  methods = list(
    initialize = function() {
      ci<<- NULL
      iteration.best<<- c()
    },

    setBest = function(v) {
      overall.best<<- v
    },

    getBest = function() {
      overall.best
    },

    addPartialBest = function(iteration, solution) {
      iteration.best<<- rbind(iteration.best, c(iteration,solution))
    },

    getPartialBest = function() {
      iteration.best
    },

    CI = function() {
      if(is.null(ci)) {
        v<- c()
        F<- function(x, i) {
          c( mean(x[i]), (sd(x[i])/sqrt(length(x[i])))^2 )
        }
        x<- as.matrix(iteration.best[Magnitude( unlist(iteration.best[,"fitness"]) ) == Magnitude( unlist(overall.best["fitness"])),])
        for(i in 2:(length(x[1,])-2)) {
          print(unlist(x[,i]))
          b<- boot(unlist(x[,i]), F, R = 10000)
          ###colnames(x)[i],
          v<- rbind(v, boot.ci(b, conf = 0.99, type = c("norm", "basic", "perc", "stud")))
        }
        ##b<- boot(x, F, R = 1000)
        ##ci<<- boot.ci(b, conf = 0.99)
        ci<<- v
      }
      ##ci
      ci
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
      setValue("iterations", 50)
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
      setValue("t0", 1)
      setValue("t.min", 10^-10)
      setValue("L", 144)
      setValue("d", 0.05)
      setValue("max.accept",32)
      setValue("max.reject",250)
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
      setValue("n.ants", 32)  ## The number of simulated ants
      setValue("k", 16)       ## The archive size
      setValue("q", 0.2)      ## Locality of the search process
      setValue("Xi", 0.5)     ## Equivalent to evaporation rate, higher Xi reduce convergence speed
    }
  )
)

#' @title OptionsSDA
#'
#' @description Options for Serial Dilutions method
#'
#' @field dilutions The desired dilutions
#'
#' @importFrom methods new
#' @export OptionsSDA
#' @exportClass OptionsSDA
OptionsSDA<- setRefClass("OptionsSDA", contains = "Options",
  fields = list(
  ),

  methods = list(
    initialize = function() {
      callSuper()
      setType("sda")
      setValue("iterations", 100)
      setValue("mu", 0.7641)    ## shaking ratio
      setValue("kkappa", 0.1)   ## Dilution factor
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

    sda={
      optimization.fun<- abm.sda
    },

    acor={
      optimization.fun<- abm.acor
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


  Pg<- Pi<- x<- initSolution(objective$parameters,N)
  lbest<- pbest<- objective$EvalFitness(x)
  Pi<- x

  for(i in 1:nrow(x)){
    lbest[i,]<- pso.lbest(i,pbest,f.neighborhood)
    Pg[i,]<- Pi[lbest[i,"pset"],]
  }

  vi<- Pi * 0.1

  for(index in 1:iterations) {
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
      if(objective$isConverged(gbest[1,"fitness"])) break
    }

    ## --- Storing the best of this iteration
    estimates$addPartialBest(index, pso.best(lbest, Pg))
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
  T0<- options$getValue("t0")
  TMIN<- options$getValue("t.min")
  L<- options$getValue("L")
  d<- options$getValue("d")
  max.accept<- options$getValue("max.accept")
  max.reject<- options$getValue("max.reject")
  f.neighborhood<- options$neighborhoodFunction()
  f.temp<- options$getTemperatureF()

  ## Generates an initial solution
  S0<- S<- initSolution(objective$parameters,1)
  f0<- (f<- objective$EvalFitness(S0))
  C<- C0<- f[1,"fitness"]

  v.accept<- 0
  v.reject<- 0

  for(k in 1:iterations) {
    t<- f.temp(T0, k)
    for(l in 1:L) {
      ## Evaluate a neighbor of S
      S1<- f.neighborhood(objective,S,d)
      f1<- objective$EvalFitness(S1)
      delta<- (C1<- f1[1,"fitness"]) - C

      if(delta < 0) {
        v.reject<- 0
        v.accept<- v.accept + 1
        S<- S1; f<- f1; C<- C1
        if(C < C0) { S0<- S; f0<- f; C0<- C }
      } else {
        if(runif(1,0,1) < exp(-delta/t)) { v.accept<- v.accept + 1; S<- S1; f<- f1; C<- C1 }
        else {v.reject<- v.reject + 1}
      }
      if(v.accept > max.accept) {v.accept<- 0; break}
      if(v.reject > max.reject) { S<- S0; C<- C0; v.reject<- 0}
    }

    ## --- Storing the best of this iteration
    estimates$addPartialBest(k, merge(S0,f0[1,c("pset","fitness")]))

    if(t <= TMIN) break;
    if(objective$isConverged(f0[1,"fitness"])) break
  }

  estimates$setBest(merge(S0,f0[1,c("pset","fitness")]))
  estimates
  #v<- merge(S0,f0[1,c("pset","fitness")])
  #return(v)
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
    newS[,i]<- newS[,i] + newS[,i] * runif(1,-1,1) * distance
    ####newS[,i]<- newS[,i] + 0.01 * f.range(k) * rnorm(1)
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
  S<- initSolution(objective$parameters,n.ants)

  C<- objective$EvalFitness(S)
  T<- acor.archive(S, C, W, k)

  for(index in 1:iterations) {
    s.s<- acor.S(T)
    s.sd<- acor.sigma(Xi, k, T)

    ## --- AntBasedSolutionConstruction
    S<- acor.updateants(S, n.ants, W, s.s, s.sd)
    C<- objective$EvalFitness(S)

    ## --- PheromoneUpdate
    T<- acor.archive(S, C, W, k)

    ## --- Storing the best of this iteration
    estimates$addPartialBest(index, T[1,])

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
  T$fitness
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
  T$w
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
## --- Simulated Dilution Approximation -
##
## ##################################################################


#' @title Simulated Dilution Approximation
#'
#' @description This function tries to provide a better approximation to the
#' best solution when no information is available on the correct range of
#' input parameters for the objective function. Basically the function search
#' for the best cost n replications and adjust the range based on that value
#' of best particles and a 10^-1 of initially provided range.
#'
#' @param objective An instance of ObjectiveFunction (or subclass) class \link{ObjectiveFunction}
#' @param options An apropiate instance from a sublclass of \link{Options} class
#'
#' @export
abm.sda<- function(objective, options= NULL) {
  ## Handling the heuristic specific options
  options<- OptionsFactory("sda", options)

  ## --- Creating the estimation object for returning results
  estimates<- Estimates$new()

  iterations<- options$getValue("iterations")
  k<- objective$ParametersSize()
  mu<- options$getValue("mu")           ## shaking ratio
  kkappa<- options$getValue("kkappa")   ## Dilution factor

  ## --- Initialize the solution
  s<- initSolution(objective$parameters, 4*k, "lhs")
  c<- objective$EvalFitness(s)
  s0<- sda.solution(s, c)

  for(index in 1:iterations) {
    s<- s0[,1:k]

    ##print(length(s0[1,]))

    mix<- sda.shaking1(s, c, mu)
    solute<- sda.solute(mix[,1:k])
    s<- sda.mixing(s, c, solute, kkappa)
    s<- enforceBounds(s, objective$parameters)

    ##print(apply(mixing[,1:length(s[1,])],2,gm.mean))
    #s<- sda.stock(s, c, mu, kkappa)

    c<- objective$EvalFitness(s)
    s1<- sda.solution(s,c)

    elog.debug("iteration best [%g]",s1[1,"fitness"])
    s0<- sda.dilute(s0, s1, 0.01)

    ## --- Storing the best of this iteration
    estimates$addPartialBest(index, s0[1,])

    ## Check for algorithm convergence
    if(objective$isConverged(s0[1, "fitness"])) break
  }

  estimates$setBest(s0[1,])
  estimates
}

#' @title sda.solution
#'
#' @description Sort solutions by its cost
#'
#' @param s Problem solution
#' @param f The function evaluation for s
#'
#' @export
sda.solution<- function(s, f) {
  solution<- cbind(s,f)
  as.matrix(solution[with(solution,order(fitness)),])
}

#' @title sda.shaking1
#'
#' @description This function 'mix' the elements present in the solution. The
#' parameter 'mu' controls the intensity of mixing. Low values give preference
#' to best solutions and high values make the values being select randomly.
#'
#' @param s The Problem solution
#' @param f The function evaluation for s
#' @param mu The mixing intensity ratio, from 0 to 1. The shaking intensity
#' controls de the probability of chosing a 'heavier' values.
#'
#' @export
sda.shaking1<- function(s, f, mu) {
  P<- function(p, x) { (p^x) }
  k<- length(s[,1])
  solution<- sda.solution(s, f)
  p<- c( rep(1,(k/2)), P(mu,((k/2)+1):k) )

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

#' @title sda.shaking2
#'
#' @description This function 'mix' the elements present in the solution. The
#' parameter 'mu' controls the intensity of mixing. Low values give preference
#' to best solutions and high values make the values being select randomly.
#'
#' @param s The Problem solution
#' @param f The function evaluation for s
#' @param mu The mixing intensity ratio, from 0 to 1. The shaking intensity
#' controls de the probability of chosing a 'heavier' values.
#'
#' @export
sda.shaking2<- function(s, f, mu) {
  P<- function(p, x) { (p^x) }
  k<- length(s[,1])
  solution<- sda.solution(s, f)
  p<- c( rep(1,(k/2)), P(mu,((k/2)+1):k) )

  solution[sample(1:k, size = (k/2), prob = p),]
}

#' @title sda.solute
#'
#' @description This function mimics the pipetting a solute
#' for multiplicatively diluting it in a new solution. Basically,
#' it is calculating the geometric mean of problem parameters in
#' the mixed solute.
#'
#' @param mix The solution mix
#'
#' @return The 'solute' The geometric mean of mix columns
#' @export
sda.solute<- function(mix) {
  n<- length(mix[1,])
  apply(mix[,1:n],2,gm.mean)
}

#' @title sda.mixing
#'
#' @description Mix solute
#'
#' @param s The Problem solution
#' @param f The function evaluation for s
#' @param solute The solute generated with sda.solute.
#' @param kkappa The dilution factor
#'
#' @export
sda.mixing<- function(s, f, solute, kkappa) {
  randomize<- function(x, kkappa) { stats::rnorm(1, x, exp(abs(x * kkappa))) }
  m<- length(s[,1])
  n<- length(s[1,])

  solution<- cbind(s,f)
  summatory<- sum(solution[,"fitness"])

  b<- which(solution[,"fitness"] == min(solution[,"fitness"]),arr.ind = TRUE)

  ##print(apply(s[,],2,sd)/apply(s[,],2,mean))
  print(apply(s[,],2,gm.mean))

  stock<- c()
  for(i in 1:m) {
    si<- matrix(s[i,],1,n)
    sf<- 10^trunc(log(solution[i,"fitness"],10))
    ##if(i == b) {
    if(i == b) {
      stock<- rbind(stock,si + si * runif(n,-1, 1))
      #stock<- rbind(stock, apply(rbind(si + si * runif(n,-1/kkappa, 1/kkappa), solute * solution[i,"fitness"]/summatory),2,gm.mean) )
    } else {
      stock<- rbind(stock, apply(rbind(si + runif(n,-1/kkappa, 1/kkappa), solute * solution[i,"fitness"]/summatory),2,mean) )
      #stock<- rbind(stock, apply(rbind(apply(si, 2, randomize, kkappa=0.10), solute * solution[i,"fitness"]/summatory),2,gm.mean) )
    }
  }
  as.data.frame(stock)
}

#' @title sda.dilute
#'
#' @description Dilute the old solute in the new solution
#'
#' @param s0 The old solution
#' @param s1 The new solution
#' @param kkappa The dilution factor
#'
#' @export
sda.dilute<- function(s0, s1, kkappa) {
  assert(length(s0[1,]) == length(s1[1,]),"Invalid solution!")
  k<- length(s0[,1])
  for(i in 1:k) {
    ##print(sprintf("[i=%g] f0=%g, f1=%g", i, s0[i,"fitness"], s1[i,"fitness"]))
    if(s1[i,"fitness"]< s0[i,"fitness"]) {
        s0[i,]<- s1[i,]
    } else {
      if(i > 1 && runif(1) < kkappa) {
        s0[i,]<- s1[i,]
      }
    }
  }
  s0
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
#' @param sampling The population sampling scheme (random|lhs)
#'
#' @return A random set of solutions
#'
#' @importFrom rrepast AoE.RandomSampling AoE.LatinHypercube
#' @export
initSolution<- function(parameters, N=20, sampling="random") {
  switch(sampling,
    lhs={
      rrepast::AoE.LatinHypercube(N,parameters)
    },

    {
      rrepast::AoE.RandomSampling(N,parameters)
    }
  )
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
      sda = { v<- OptionsSDA$new() },
      acor= { v<- OptionsACOR$new() },
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
  trunc(log(v,10))
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
