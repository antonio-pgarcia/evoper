##================================================================================
## This file is part of the evoper package - EvoPER
##
## (C)2016, 2017 Antonio Prestes Garcia <@>
## For license terms see DESCRIPTION and/or LICENSE
##
## @file: abm-ees1.R
##
## This file contains the EvoPER Evolutionary Strategy 1 (EES1) metaheuristic
## and its associated auxiliary functions.
##================================================================================



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


## ##################################################################
##
## --- EvoPER Evolutionary Strategy 1 helper functions -
##
## ##################################################################


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
