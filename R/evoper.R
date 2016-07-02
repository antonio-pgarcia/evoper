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
    bestS = 'ANY',
    bestF = 'ANY',
    counter = 'ANY'),

  methods = list(
    initialize = function(funct) {
      object<<- NULL
      objective<<- funct
      parameters<<- NULL
      value<<- NULL
      tolerance<<- .Machine$double.eps^0.30
      converged<<- FALSE
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

    Value = function(v = NULL) {
      if(!is.null(v)) {
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

    EvaluateV = function(swarm) {
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

##
## ----- Options class
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

    setNeighborhoodF = function(v) {
      neighborhood<<- v
    },

    getNeighborhoodF = function() {
      neighborhood
    },

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
      setValue("N",16)
      setValue("phi1",1.193)
      setValue("phi2",1.193)
      setValue("W",0.721)
      setNeighborhoodF(pso.neighborhood.KN)
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
      setValue("T0", 1)
      setValue("TMIN", 10^-15)
      setValue("L", 144)
      setValue("d", 0.05)
      setValue("max.accept",32)
      setValue("max.reject",250)
      setNeighborhoodF(saa.neighborhoodN)
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
    optf = 'ANY'
  ),

  methods = list(
    initialize = function() {
      callSuper()
      setType("sda")
      setValue("iterations", 12)
      setValue("dilutions", 8)
      setOptF(abm.pso)
    },

    setOptF = function(v) {
      optf<<- v
    },

    getOptF = function() {
      optf
    },

    getOptions = function() {
      if(identical(optf,abm.pso)) {
        opt<- OptionsPSO$new()
      }

      if(identical(optf,abm.saa)) {
        opt<- OptionsSAA$new()
      }

      opt$setValue("iterations",getValue("iterations"))
      return(opt)
    }

  )
)

##
## ----- Entry point function
##

#' @title extremize
#'
#' @description Entry point for optimization functions
#'
#' @param type The optimization method (pso,saa,sda,aco)
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

    aco={
    },

    {
      stop("Invalid optimization function!")
    }
  )

  optimization.fun(objective, options)
}

##
## ----- Particle Swarm Optimization functions
##


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

  if(is.null(options)) {
    options<- OptionsPSO$new()
  } else {
    if(options$getType() != "pso") stop(paste("Invalid option of type [", options$getType(),"]"))
  }

  iterations<- options$getValue("iterations")
  N<- options$getValue("N")
  phi1<- options$getValue("phi1")
  phi2<- options$getValue("phi2")
  W<- options$getValue("W")
  f.neighborhood<- options$getNeighborhoodF()


  Pg<- Pi<- x<- initSolution(objective$parameters,N)
  lbest<- pbest<- objective$EvaluateV(x)
  Pi<- x

  for(i in 1:nrow(x)){
    lbest[i,]<- pso.lbest(i,pbest,f.neighborhood)
    Pg[i,]<- Pi[lbest[i,"pset"],]
  }

  vi<- Pi * 0.1

  for(index in 1:iterations) {
    vi<- pso.Velocity(W,vi,phi1,phi2,Pi,Pg,x)
    x<- enforceBounds((x + vi), objective$parameters)

    f1<- objective$EvaluateV(x)

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
 }
 return(pso.best(lbest, Pg))
}

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

  if(is.null(options)) {
    options<- OptionsSDA$new()
  } else {
    if(options$getType() != "sda") stop(paste("Invalid option of type [", options$getType(),"]"))
  }

  iterations<- options$getValue("iterations")
  steps<- options$getValue("dilutions")

  o<- objective$copy()

  f.optimization<- options$getOptF()

  v<- v0<- f.optimization(o, options$getOptions())
  f0<- v0$fitness

  for(s in 1:steps) {

    o1<- o$copy()
    for(k in o$getParameterNames()) {
      min0<- as.numeric(o$getParameterValue(k,"min"))
      max0<- as.numeric(o$getParameterValue(k,"max"))

      #max1<- max0 * .90^s * runif(1, .6, 1)
      #min1<- min0 * .90^s * runif(1, .6, 1)
      #print(paste("k=", k, "v",v[1,k]))

      #U<- runif(1, .9, 1)
      U<- 1
      rmax1<-  (max0 * .90^s * U)
      rmin1<-  (min0 * .90^s * U)
      max1<- v[1,k] + (rmax1 - rmin1)/exp(1)
      min1<- v[1,k] - (rmax1 - rmin1)/exp(1)

      o1$parameters[ o1$parameters[,"name"] == k,"min"]<- ifelse(min1 > min0, min1, (runif(1) * min0))
      o1$parameters[ o1$parameters[,"name"] == k,"max"]<- ifelse(max1 < max0, max1, (runif(1) * max0))
    }

    print("************************************************************")
    print(o1$parameters)
    print("------------------------------------------------------------")
    print(paste("[",s,"]","f0=",f0,"f1=",v$fitness))
    print("************************************************************")

    #v<- abm.pso(o1, iterations=13)
    #v<- extremize("saa", o1, iterations=13)
    v<- f.optimization(o1, options$getOptions())

    if(v$fitness < f0) {
      f0<- v$fitness
      v0<- v
      o<- o1
    }
  }

  eval.parent(substitute(objective<-o))
  return(v0)
}

#' @export
abm.sda1<- function(objective, options= NULL) {

  n.sample<- function(p) {
    (nrow(p) * 5 + nrow(p))
  }

  x.sample<- function(n,p) {
    S<- rrepast::AoE.LatinHypercube(n, p)
    f<- objective$EvaluateV(S)
    v<- merge(S,f,by=0)
    v[with(v,order(fitness)),]
  }

  dilutions<- 10
  c0<- n.sample(objective$parameters)
  c1<- floor((c0-1)/2)

  x1<- c()
  x0<- NULL
  print(paste("parameters: ", dilutions, "c0=", c0, "c1=", c1))

  for(i in 1:dilutions) {
    x<- x.sample(c0,objective$parameters)

    if(is.null(x0)){
      x0<- x
    } else{
      if(x[1,"fitness"]< x0[1,"fitness"]) {
        x0<- x[1,]
      }
    }

    fsum<- sum(x[,"fitness"])
    #P<- x[ (U0<- runif(1,1,c0)) ,"fitness"]/fsum
    P<- sum(x[1:c1,"fitness"])/fsum

    print(paste("***** generation best S for dilution: ",dilutions," *****"))
    print(x[1,])

    for(k in objective$getParameterNames()) {
      objective$parameters[objective$parameters[,"name"] == k,"min"]<- min(x[1:c1,k])
      objective$parameters[objective$parameters[,"name"] == k,"max"]<- max(x[1:c1,k])
    }
  }

  #print(x1[with(x1,order(fitness)),])
  #return(x0[1,])
  return(x0)
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

#' @title initSolution
#'
#' @description Creates the initial Solution population
#' taking into account the lower an upper bounds of
#' provided experiment factors.
#'
#' @param parameters The Objective Function parameter list
#' @param N The size of Solution population
#'
#' @return A random set of solutions
#'
#' @importFrom rrepast AoE.RandomSampling
#' @export
initSolution<- function(parameters, N=20) {
  rrepast::AoE.RandomSampling(N,parameters)
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

enforceBounds1<- function(particles, factors) {
  k<- rrepast::GetFactorsSize(factors)

  bounds<- function(i, p, f) {
    if(nrow(p) > 1) {
      x<- p[,i]
    } else {
      x<- p[i]
    }

    lb<- as.numeric(f[i,"min"]);
    ub<- as.numeric(f[i,"max"]);

    for(j in 1:length(x)) {
      if( x[j] < lb || x[j] > ub) {
        x[j]<- runif(1,lb,ub)
      }
    }
    return(x)
  }
  v<- as.data.frame(sapply(1:k,bounds,p=particles,f=factors))
  names(v)<- factors[,"name"]
  return(v)
}


##
## ----- Simulated Annealing Algorithm
##


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
#abm.saa<- function(objective, iterations= 1000, T0=1, TMIN=10^-12,  L=144, d=0.16, f.neighborhood=saa.neighborhoodN, f.temp=saa.bolt) {
abm.saa<- function(objective, options= NULL) {

  if(is.null(options)) {
    options<- OptionsSAA$new()
  } else {
    if(options$getType() != "saa") stop(paste("Invalid option of type [", options$getType(),"]"))
  }

  iterations<- options$getValue("iterations")
  T0<- options$getValue("T0")
  TMIN<- options$getValue("TMIN")
  L<- options$getValue("L")
  d<- options$getValue("d")
  max.accept<- options$getValue("max.accept")
  max.reject<- options$getValue("max.reject")
  f.neighborhood<- options$getNeighborhoodF()
  f.temp<- options$getTemperatureF()

  ## Generates an initial solution
  S0<- S<- initSolution(objective$parameters,1)
  f0<- (f<- objective$EvaluateV(S0))
  C<- C0<- f[1,"fitness"]

  v.accept<- 0
  v.reject<- 0

  for(k in 1:iterations) {
    t<- f.temp(T0, k)
    for(l in 1:L) {
      ## Evaluate some neighbor of S
      S1<- f.neighborhood(objective,S,d)
      f1<- objective$EvaluateV(S1)
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
    #print(merge(S0,f0[1,c("pset","fitness")]))
    #print(merge(S1,f1[1,c("pset","fitness")]))
    #print(paste("rejects=", v.reject))
    if(t <= TMIN) break;
    if(objective$isConverged(f0[1,"fitness"])) break
  }

  v<- merge(S0,f0[1,c("pset","fitness")])
  return(v)
}

#' @export
abm.saa.1<- function(objective, options= NULL) {

  if(is.null(options)) {
    options<- OptionsSAA$new()
  } else {
    if(options$getType() != "saa") stop(paste("Invalid option of type [", options$getType(),"]"))
  }

  t<- options$getValue("T0")
  TMIN<- options$getValue("TMIN")
  L<- options$getValue("L")
  d<- options$getValue("d")
  max.accept<- options$getValue("max.accept")
  max.reject<- options$getValue("max.reject")
  f.neighborhood<- options$getNeighborhoodF()
  f.temp<- options$getTemperatureF()

  max.run<- 150
  max.accept<- 15
  max.reject<- 250
  e.norm<- 1e-2

  c.run<- 0
  c.accept<- 0
  rejects<- 0

  ## Generates an initial solution
  S0<- S<- initSolution(objective$parameters,1)
  C0<- C<- (f0<- (f<- objective$EvaluateV(S0)))[1,"fitness"]

  while(t > TMIN) {
    if ( c.accept > max.accept || c.run > max.run ) {
      t<- f.temp(t, 0)
      c.accept<- 0
      c.run<- 0
    }

    S<- f.neighborhood(objective,S0,d)
    f<- objective$EvaluateV(S)
    delta<- (C<- f[1,"fitness"]) - C0

    if(-delta < e.norm) {
      f0<- f; C0<- C; S0<- S
      c.accept<- c.accept + 1
    } else {
      if(exp(-delta/t) > runif(1,0,1)) {
        f0<- f; C0<- C; S0<- S
        c.accept<- c.accept + 1
      }
    }

    if(rejects > options$getValue("max.reject")) { break }
  }

  v<- merge(S0,f0[1,c("pset","fitness")])
  return(v)
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
    #newS[,i]<- newS[,i] + newS[,i] * runif(1,-1,1) * distance
    newS[,i]<- newS[,i] + newS[,i] * rnorm(1) * distance
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


##
## ----- Ant colony optimization
##
abm.aco<- function() {
  ##Still under test
}


##
## ----- Test functions
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
slope<- function(x, y,i) {
  if(i > 1 && i < length(y)) {
    v<- 0.5 * ((y[i]-y[i-1])/(x[i]-x[i-1]) + (y[i+1] - y[i])/(x[i+1]-x[i]))
  } else {
    v<- ifelse(i == 1,(y[i+1] - y[i])/(x[i+1]-x[i]),(y[i]-y[i-1])/(x[i]-x[i-1]))
  }
  v
}
