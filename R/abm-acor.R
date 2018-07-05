##================================================================================
## This file is part of the evoper package - EvoPER
##
## (C)2016, 2017 Antonio Prestes Garcia <@>
## For license terms see DESCRIPTION and/or LICENSE
##
## @file: abm-acor.R
##
## This file contains the Ant colony optimization for continuous domains (ACOr)
## metaheuristic and its associated auxiliary functions.
##================================================================================


#' @title Ant colony optimization for continuous domains
#'
#' @description An implementation of Ant Colony Optimization algorithm
#' for continuous variables.
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
#'  extremize("acor", f)
#' }
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
  elog.info("Options(%s): %s", options$getType(), options$toString())

  ## --- Adjusting parameter types
  parameterz<- paramconverter(objective$parameters, options$isDiscrete(), options$getLevels())
  objective$parameters<- parameterz

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

  tmp<- T[1,-which(names(T[1,]) %in% c("w"))]
  tmp$pset<- row.names(tmp)
  estimates$setBest(tmp)
  estimates
}



## ##################################################################
##
## --- Ant colony optimization helper functions -
##
## ##################################################################


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
