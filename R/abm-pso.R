##================================================================================
## This file is part of the evoper package - EvoPER
##
## (C)2016, 2017 Antonio Prestes Garcia <@>
## For license terms see DESCRIPTION and/or LICENSE
##
## @file: abm-pso.R
##
## This file contains the Particle Swarm Optimization (PSO) metaheuristic and its
## associated auxiliary functions.
##================================================================================


#' @title abm.pso
#'
#' @description An implementaion of Particle Swarm Optimization
#' method for parameter estimation of Individual-based models.
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
#'  extremize("pso", f)
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


## ##################################################################
##
## --- Particle Swarm Optimization helper functions -
##
## ##################################################################


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
