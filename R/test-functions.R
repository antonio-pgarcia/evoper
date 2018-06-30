##================================================================================
## This file is part of the evoper package - EvoPER
##
## (C)2016, 2017 Antonio Prestes Garcia <@>
## For license terms see DESCRIPTION and/or LICENSE
##
## @file: test-functions.R
##
## This file contains the standard functions for benchmarking optimization
## metaheuristics.
##================================================================================



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


#' @title f0.rosenbrockn
#'
#' @description The rosenbrock function of N variables for testing
#' optimization methods. The global optima for the function is given by
#' xi = 1, forall i E {1...N}, f(x) = 0.
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
f0.rosenbrockn<- function(...) {
  x<- list(...)
  f1.rosenbrockn(unlist(x))
}

#' @title f1.rosenbrockn
#'
#' @description The rosenbrock function of N variables for testing
#' optimization methods. The global optima for the function is given by
#' xi = 1, forall i E {1...N}, f(x) = 0.
#'
#' @param x The vector of function parameters
#'
#' @return The function value
#'
#' @references
#'
#' http://deap.gel.ulaval.ca/doc/dev/api/benchmarks.html
#'
#' @export
f1.rosenbrockn<- function(x) {
  ssum<- 0
  for(i in 1:(length(x)-1)) {
    ssum<- ssum + (1 - x[i])^2 + 100 * (x[i+1] - x[i]^2)^2
  }
  ssum
}

#' @title f0.rosenbrock4
#'
#' @description The rosenbrock function of 4 variables for testing
#' optimization methods. The global optima for the function is given by
#' xi = 1, forall i E {1...N}, f(x) = 0.
#'
#' @param x1 The first function variable
#' @param x2 The second function variable
#' @param x3 The third function variable
#' @param x4 The fourth function variable
#'
#' @return The function value
#'
#' @export
f0.rosenbrock4<- function(x1, x2, x3, x4) {
  f1.rosenbrockn(c(x1, x2, x3, x4))
}

#' @title f0.schwefel
#'
#' @description The schwefel function of N variables for testing optimization
#' methods. The global optima for the function is given by
#' xi = 420.96874636, forall i E {1...N}, f(x) = 0. The range of xi is [-500,500]
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
f0.schwefel<- function(...) {
  x<- list(...)
  f1.schwefel(unlist(x))
}

#' @title f1.schwefel
#'
#' @description The schwefel function of N variables for testing optimization
#' methods. The global optima for the function is given by
#' xi = 420.96874636, forall i E {1...N}, f(x) = 0. The range of xi is [-500,500]
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
f1.schwefel<- function(x) {
  N<- length(x)
  ssum<- 0
  for(i in 1:N) {
    ssum<- ssum + x[i] * sin(sqrt(abs(x[i])))
  }
  abs(round((418.9828872724339 * N - ssum), sqrt(.Machine$double.eps)))
}

#' @title f0.schwefel4
#'
#' @description The schwefel function of N variables for testing optimization
#' methods. The global optima for the function is given by
#' xi = 420.96874636, forall i E {1...N}, f(x) = 0. The range of xi is [-500,500]
#'
#' @param x1 The first function variable
#' @param x2 The second function variable
#' @param x3 The third function variable
#' @param x4 The fourth function variable
#'
#' @return The function value
#'
#' @export
f0.schwefel4<- function(x1, x2, x3, x4) {
  f1.schwefel(c(x1, x2, x3, x4))
}

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

#' @title f0.schaffer
#'
#' @description The schaffer function of N variables for testing
#' optimization methods. The global optima for the function is given by
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
f0.schaffer<- function(...) {
  x<- list(...)
  f1.schaffer(unlist(x))
}

#' @title f1.schaffer
#'
#' @description The schaffer function of N variables for testing
#' optimization methods. The global optima for the function is given by
#' xi = 0, forall i E {1...N}, f(x) = 0.
#'
#' @param x The vector of function parameters
#'
#' @return The function value
#'
#' @references
#'
#' http://deap.gel.ulaval.ca/doc/dev/api/benchmarks.html
#'
#' @export
f1.schaffer<- function(x) {
  ssum<- 0
  for(i in 1:(length(x)-1)) {
    ssum<- ssum + (x[i]^2 + x[i+1]^2)^0.25 * (sin(50 * (x[i]^2 + x[i+1]^2)^0.1)^2 + 1)
  }
  ssum
}

#' @title f0.schaffer4
#'
#' @description The Schaffer function of four variables for testing optimization
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
f0.schaffer4<- function(x1, x2, x3, x4) {
  f0.schaffer(c(x1, x2, x3, x4))
}

#' @title f0.bohachevsky
#'
#' @description The Bohachevsky function of N variables for testing
#' optimization methods. The global optima for the function is given by
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
f0.bohachevsky<- function(...) {
  x<- list(...)
  f1.bohachevsky(unlist(x))
}

#' @title f1.bohachevsky
#'
#' @description The Bohachevsky function of N variables for testing
#' optimization methods. The global optima for the function is given by
#' xi = 0, forall i E {1...N}, f(x) = 0.
#'
#' @param x The vector of function parameters
#'
#' @return The function value
#'
#' @references
#'
#' http://deap.gel.ulaval.ca/doc/dev/api/benchmarks.html
#'
#' @export
f1.bohachevsky<- function(x) {
  ssum<- 0
  for(i in 1:(length(x)-1)) {
    ssum<- ssum + ( x[i]^2 + x[i+1]^2 - 0.3 * cos(3 * pi * x[i]) - 0.4 * cos(4 * pi * x[i+1]) + 0.7 )
  }
  ssum
}

#' @title f0.bohachevsky4
#'
#' @description The Bohachevsky function of four variables for testing optimization
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
f0.bohachevsky4<- function(x1, x2, x3, x4) {
  f0.bohachevsky(c(x1, x2, x3, x4))
}

#' @title f0.griewank
#'
#' @description The griewank function of N variables for testing
#' optimization methods. The global optima for the function is given by
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
f0.griewank<- function(...) {
  x<- list(...)
  f1.griewank(unlist(x))
}

#' @title f1.griewank
#'
#' @description The griewank function of N variables for testing
#' optimization methods. The global optima for the function is given by
#' xi = 0, forall i E {1...N}, f(x) = 0.
#'
#' @param x The vector of function parameters
#'
#' @return The function value
#'
#' @references
#'
#' http://deap.gel.ulaval.ca/doc/dev/api/benchmarks.html
#'
#' @export
f1.griewank<- function(x) {
  ssum<- 0
  prod<- 1
  for(i in 1:(length(x))) {
    ssum<- ssum + x[i]^2
    prod<- prod * cos(x[i]/sqrt(1))
  }
  (1/4000 * ssum - prod + 1)
}

#' @title f0.griewank4
#'
#' @description The griewank function of four variables for testing optimization
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
f0.griewank4<- function(x1, x2, x3, x4) {
  f0.griewank(c(x1, x2, x3, x4))
}

#' @title f0.ackley
#'
#' @description The ackley function of N variables for testing
#' optimization methods. The global optima for the function is given by
#' xi = 0, forall i E {1...N}, f(x) = 0.
#' Domain xi E [-32.768, 32.768], for all i = 1, ..., d
#'
#' @param ... The variadic list of function variables.
#'
#' @return The function value
#'
#' @references
#'
#' https://www.sfu.ca/~ssurjano/ackley.html
#'
#' @export
f0.ackley<- function(...) {
  x<- list(...)
  f1.ackley(unlist(x))
}

#' @title f1.ackley
#'
#' @description The ackley function of N variables for testing
#' optimization methods. The global optima for the function is given by
#' xi = 0, forall i E {1...N}, f(x) = 0.
#' Domain xi E [-32.768, 32.768], for all i = 1, ..., d
#'
#' @param x The vector of function parameters
#'
#' @return The function value
#'
#' @references
#'
#' https://www.sfu.ca/~ssurjano/ackley.html
#'
#' @export
f1.ackley<- function(x) {
  a<- 20.0
  b<- 0.2
  c<- 2.0*pi
  d<- length(x)

  v<- (-a * exp(-b * sqrt( (1/d) * sum(x^2))) - exp((1/d) * sum(cos(c*x))) + a + exp(1))
  ifelse(v < 10^-15, 0, v)
}

#' @title f0.ackley4
#'
#' @description The ackley function of four variables for testing optimization
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
f0.ackley4<- function(x1, x2, x3, x4) {
  f0.ackley(c(x1, x2, x3, x4))
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
  value<- c(x = 120, y = 12)
  parameters<- c(c1 = x1, c2 = x2, c3 = x3, c4 = x4)
  time<- seq(0, 120, by = 1)

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

#' @title predatorprey.plot0
#'
#' @description Generate a plot for the predator-prey ODE output.
#'
#' @param x1 The growth rate of prey
#' @param x2 The decay rate of predator
#' @param x3 The predating effect on prey
#' @param x4 The predating effect on predator
#' @param title The optional plot title. May be omited.
#'
#' @return An ggplot2 object
#'
#' @examples \dontrun{
#'  predatorprey.plot0(1.351888, 1.439185, 1.337083, 0.9079049)
#' }
#'
#' @importFrom reshape melt
#' @importFrom ggplot2 ggplot aes geom_line ggtitle
#' @export
predatorprey.plot0<- function(x1, x2, x3, x4, title = NULL) {
  v<- as.data.frame(predatorprey(x1, x2, x3, x4))
  v.data<- melt(v, id.vars="time", value.name="value", variable_name="species")
  p<- ggplot(data= v.data, with(v.data,aes(x=time, y=value, group = species, colour = species)))
  p + geom_line() + ggtitle(title)
}

#' @title predatorprey.plot1
#'
#' @description Simple wrapper for 'predatorprey.plot0' accepting the parameters as a list.
#'
#' @param x A list containing the values of predator/prey parameters c1, c2, c3 and c4
#' denoting respectivelly the growth rate of prey, the decay rate of predator,
#' the predating effect on prey and the predating effect on predator
#' @param title The optional plot title. May be omited.
#'
#' @return An ggplot2 object
#'
#' @examples \dontrun{
#'  rm(list=ls())
#   set.seed(27262565)
#   f<- PlainFunction$new(f0.periodtuningpp24)
#   f$Parameter(name="x1",min=0.5,max=2)
#   f$Parameter(name="x2",min=0.5,max=2)
#   f$Parameter(name="x3",min=0.5,max=2)
#   f$Parameter(name="x4",min=0.5,max=2)
#   v<- extremize("acor", f)
#'  predatorprey.plot1(v$getBest()[1:4])
#' }
#'
#' @export
predatorprey.plot1<- function(x, title = NULL) {
  predatorprey.plot0(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[3]),as.numeric(x[4]), title)
}


#' @title Period tuning for Predator-Prey base
#'
#' @description This function is an example on how EvoPER can be
#' used for estimating the parameter values in order to produce
#' oscilations with the desired period. It is not intended to be
#' used directelly, the provided wrappers should be instead.
#'
#' @param x1 The growth rate of prey
#' @param x2 The decay rate of predator
#' @param x3 The predating effect on prey
#' @param x4 The predating effecto on predator
#' @param period The desired oscilation period
#'
#' @return The solution fitness cost
#'
#' @export
f0.periodtuningpp<- function(x1, x2, x3, x4, period) {
  v<- predatorprey(x1, x2, x3, x4)
  rrepast::AoE.NRMSD(naiveperiod(ifelse(v[,"y"] < 0, .Machine$double.xmax, v[,"y"]))$period, period)
}

#' @title Period tuning of 12 time units for Predator-Prey
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
#'	f<- PlainFunction$new(f0.periodtuningpp12)
#'	f$Parameter(name="x1",min=0.5,max=2)
#'	f$Parameter(name="x2",min=0.5,max=2)
#'	f$Parameter(name="x3",min=0.5,max=2)
#'	f$Parameter(name="x4",min=0.5,max=2)
#'	extremize("pso", f)
#' }
#'
#' @export
f0.periodtuningpp12<- function(x1, x2, x3, x4) {
  f0.periodtuningpp(x1, x2, x3, x4, 12)
}

#' @title Period tuning of 24 time units for Predator-Prey
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
#'	f<- PlainFunction$new(f0.periodtuningpp24)
#'	f$Parameter(name="x1",min=0.5,max=2)
#'	f$Parameter(name="x2",min=0.5,max=2)
#'	f$Parameter(name="x3",min=0.5,max=2)
#'	f$Parameter(name="x4",min=0.5,max=2)
#'	extremize("pso", f)
#' }
#'
#' @export
f0.periodtuningpp24<- function(x1, x2, x3, x4) {
  f0.periodtuningpp(x1, x2, x3, x4, 24)
}

#' @title Period tuning of 48 time units for Predator-Prey
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
#'	f<- PlainFunction$new(f0.periodtuningpp24)
#'	f$Parameter(name="x1",min=0.5,max=2)
#'	f$Parameter(name="x2",min=0.5,max=2)
#'	f$Parameter(name="x3",min=0.5,max=2)
#'	f$Parameter(name="x4",min=0.5,max=2)
#'	extremize("pso", f)
#' }
#'
#' @export
f0.periodtuningpp48<- function(x1, x2, x3, x4) {
  f0.periodtuningpp(x1, x2, x3, x4, 48)
}

#' @title Period tuning of 72 time units for Predator-Prey
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
#'	f<- PlainFunction$new(f0.periodtuningpp24)
#'	f$Parameter(name="x1",min=0.5,max=2)
#'	f$Parameter(name="x2",min=0.5,max=2)
#'	f$Parameter(name="x3",min=0.5,max=2)
#'	f$Parameter(name="x4",min=0.5,max=2)
#'	extremize("pso", f)
#' }
#'
#' @export
f0.periodtuningpp72<- function(x1, x2, x3, x4) {
  f0.periodtuningpp(x1, x2, x3, x4, 72)
}

