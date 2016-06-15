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
#' @export ObjectiveFunction
#' @exportClass ObjectiveFunction
ObjectiveFunction<- setRefClass("ObjectiveFunction",
  fields = list(
    object = 'ANY',
    objective = 'function',
    parameters = 'ANY',
    value = 'ANY'),

  methods = list(
    initialize = function(funct) {
      object<<- NULL
      objective<<- funct
      parameters<<- NULL
      value<<- NULL
    },

    Parameter = function(name, min, max) {
      if(is.null(parameters)) {
        parameters<<- rrepast::AddFactor(c(), name= name, min= min, max= max)
      } else {
        parameters<<- rrepast::AddFactor(parameters, name= name, min= min, max= max)
      }
    },

    GetParameter = function(key) {
      parameters[which(parameters[,"name"] == key),]
    },

    GetParameterValue = function(key, name) {
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
    }

  )
)

#' @title PlainFunction
#'
#' @description PlainFunction Class
#'
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

      v<- c()
      for(i in 1:nrow(swarm)){
        p<- sapply(1:length(swarm[i,]),function(j,d) {d[,j]}, d = swarm[i,] )
        f<- do.call(objective,as.list(p))
        v<- rbind(v,cbind(pset=i,fitness=f))
      }
      v<- Value(as.data.frame(v))
    }
  )

)

#' @title RepastFunction
#'
#' @description RepastFunction class
#'
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
## ----- Particle Swarm Optimization functions
##


#' @title abm.pso
#'
#' @description An implementaion of Particle Swarm Optimization
#' method for parameter estimation of Individual-based models.
#'
#' @param objective An instance of ObjectiveFunction (or subclass) class \link{ObjectiveFunction}
#' @param iterations The total number of interactions
#' @param iwc The maximun number of iteractions without changes in the output
#' @param N The Particle Swarm size
#' @param phi1 Acceleration coefficient toward the previous best
#' @param phi2 Acceleration coefficient toward the global best
#' @param W Inertia weight or Constriction coefficient
#'
#' @examples \dontrun{
#'  f<- RepastFunction$new("c:/usr/models/BactoSim(HaldaneEngine-1.0)","ds::Output",300)
#'
#'  f$addFactor(name="cyclePoint",min=0,max=90)
#'  f$addFactor(name="conjugationCost",min=0,max=100)
#'  f$addFactor(name="pilusExpressionCost",min=0,max=100)
#'  f$addFactor(name="gamma0",min=1,max=10)
#'
#'  abm.pso(f, iwc= 10, iterations=100, N=20, phi1, phi2)
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
#' @importFrom rrepast col.sum
#' @export
abm.pso<- function(objective, iterations=100, iwc=10, N=20, phi1=2.0, phi2=2.0, W=0) {

  particles<- initSolution(objective$parameters,N)
  objective$Evaluate(particles)

  Pi<- particles
  f0<- objective$Value()
  pbest<- f0
  lbest<- f0
  Pg<- Pi

  for(i in nrow(f0)){
    lbest[i,]<- pso.lbest(i,pbest,pso.neighborhood.K4)
    Pg[i,]<- Pi[lbest[i,"pset"],]
  }

  vi<- Pi * 0.2

  fitness.v<- cbuf(c(),(pso.best(lbest, Pg))[1,"fitness"], iwc)

  for(index in 1:iterations) {

    vi<- pso.Velocity(W,vi,phi1,phi2,Pi,Pg,particles)
    particles<- enforceBounds((particles + vi), objective$parameters)


    objective$Evaluate(particles)
    f1<- objective$Value()

    for(i in nrow(f1)){
      if(f1[i,"fitness"] < pbest[i,"fitness"]) {
        pbest[i,"fitness"]<- f1[i,"fitness"]
        Pi[i,]<- particles[i,]
      }

      gbest<- pso.lbest(i,f1,pso.neighborhood.K4)
      if(gbest[1,"fitness"] < lbest[i,"fitness"]) {
        lbest[i,]<- gbest
        Pg[i,]<- particles[lbest[i,"pset"],]
      }
    }


    ## Exit if the last iwc iteractions do no produce changes
    ## in the output.
    f.v<- (pso.best(lbest, Pg))[1,"fitness"]
    if(index > iwc && abs(mean(fitness.v)-f.v) < 0.0001) break;
    fitness.v<- cbuf(fitness.v, f.v, iwc)
  }

  #-- print(paste("index=", index, "fitness=", fitness.v))
  return(pso.best(lbest, Pg))
}

#' @title rangesearch.pso
#'
#' @description This function tries to provide a better approximation to the
#' best solution when no information is available on the correct range of
#' input parameters for the objective function. Basically the function search
#' for the best cost n replications and adjust the range based on that value
#' of best particles and a 10% of initially provided range.
#'
#' @param aproximations The number of aproximations
#' @param replications The number of repetitions of each aproximation
#' @param objective An instance of ObjectiveFunction (or subclass) class \link{ObjectiveFunction}
#' @param iterations The total number of interactions
#' @param iwc The maximun number of iteractions without changes in the output
#' @param N The Particle Swarm size
#' @param phi1 Acceleration coefficient toward the previous best
#' @param phi2 Acceleration coefficient toward the global best
#' @param W Inertia weight or Constriction coefficient

#'
#' @export
rangesearch.pso<- function(aproximations=10, replications=4, objective, iterations=30, iwc=10, N=10, phi1=2.0, phi2=2.0, W=0) {

  o<- objective$copy()

  for(a in 1:aproximations) {
    v<- c()

    for(r in 1:replications) {
      v<- rbind(v,abm.pso(o, iterations, iwc, N, phi1, phi2,W))
    }

    indz<- which(v == min(v$fitness), arr.ind = TRUE)
    indz<- as.integer(indz[,"row"])
    vv<- v[indz,]

    for(k in o$parameters[,"name"]) {
      value<- vv[,k]
      p.min<- as.numeric(objective$parameters[ objective$parameters[,"name"] == k,"min"])
      p.max<- as.numeric(objective$parameters[ objective$parameters[,"name"] == k,"max"])

      v.max<- value + (p.max - p.min)*0.1
      v.min<- value - (p.max - p.min)*0.1

      o$parameters[ o$parameters[,"name"] == k,"min"]<- ifelse(v.min > p.min, v.min, p.min)
      o$parameters[ o$parameters[,"name"] == k,"max"]<- ifelse(v.max < p.max, v.max, p.max)
    }
  }
  return(v)
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
  m<- matrix(seq(1,n),nrow=(ceiling(sqrt(n))))
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
pso.lowerBound<- function(particles, factors) {
  k<- rrepast::GetFactorsSize(factors)
  v<- as.data.frame(sapply(1:k, function(p) {ifelse(particles[,p] > as.numeric(factors[p,"min"]),particles[,p],as.numeric(factors[p,"min"]))}))
  names(v)<- factors[,"name"]
  return(v)
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
  v<- as.data.frame(sapply(1:k, function(p) {ifelse( particles[,p] < as.numeric(factors[p,"max"]),particles[,p],as.numeric(factors[p,"max"]))}))
  names(v)<- factors[,"name"]
  return(v)
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
#' @param T0 The initial temperature
#' @param TMIN The final temperature
#' @param L The temperature length
#' @param alpha The cooling ratio
#' @param d The neighborhood distance. The default value is 0.1 (10\%) of provided parameter range.
#'
#' @return The best solution. The first row is the best of all
#' solutions and the second row the current best solution.
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
#' @export
abm.saa<- function(objective, T0, TMIN,  L, alpha, d=0.1) {

  ## Initial solution
  S<- initSolution(objective$parameters,1)
  objective$Evaluate(S)
  f<- objective$Value()
  C<- f[1,"fitness"]

  ## The best solution found
  SS<- S
  ff<- f

  ## Initialize temperature t
  t<- T0

  ## Loop until t greater than minimun temperature
  while(t > TMIN) {

    ## Temperature Length (L) loop
    for(l in 1:L) {
      ## Evaluate some neighbor of S
      S1<- saa.neighborhood.t1(objective,S,d,1)
      objective$Evaluate(S1)
      f1<- objective$Value()
      C1<- f1[1,"fitness"]

      DELTA<- C1 - C

      ## Downhill
      if(DELTA <= 0) {
        S<- S1
        f<- f1
        C<- C1

        SS<- S1
        ff<- f1
      } else {
        ## Uphill with P = exp^-DELTA/t
        if(runif(1,0,1) < exp(-DELTA/t)) {
          S<- S1
          f<- f1
          C<- C1
        }
      }
    }
    ## Cooling scheme
    t<- alpha * t
  }
  v<- c()
  v<- rbind(v,merge(SS,ff[1,c("pset","fitness")]))
  v<- rbind(v,merge(S,f[1,c("pset","fitness")]))
  return(v)
}

#' @title saa.neighborhood.t1
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
saa.neighborhood.t1<- function(f, S, d, n) {
  assert(n > 0 && n <= ncol(S),"Invalid number of parameters to be perturbed!")
  newS<- S
  for(i in sample(1:ncol(S),n)) {
    k<- colnames(S)[i]
    distance<- (as.numeric(f$GetParameterValue(k,"max")) - as.numeric(f$GetParameterValue(k,"min"))) * d
    newS[,i]<- newS[,i] + newS[,i] * runif(1,-1,1) * distance
  }
  enforceBounds(as.data.frame(newS), f$parameters)
}
