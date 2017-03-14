##================================================================================
## This file is part of the evoper package - EvoPER
##
## (C)2016, 2017 Antonio Prestes Garcia <@>
## For license terms see DESCRIPTION and/or LICENSE
##
## @file: helper-functions.R
##
## This file contains a collection of auxiliary functions for general tasks
## required by the package main functionality as well as for plot the
## metaheuristics results.
##================================================================================

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


## ##################################################################
##
## ---------------------- General functions -------------------------
##
## ##################################################################


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

#' @title random.whell
#'
#' @description A simple randon seed generator
#'
#' @return A random number for seeding
#' @export
random.wheel<- function() {
  set.seed(as.numeric(format(Sys.time(), "%OS3")) * 1000)
  wheel<- runif(1,1,10^6)
  for(i in 1:wheel) {
    runif(1)
    if(runif(1) < runif(1)) {
      v<- trunc(runif(1,1,.Machine$integer.max * 10^-3))
      break
    }
  }
  v
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


## ##################################################################
##
## -------------- Helper functions for ploting --------------
##
## ##################################################################


#' @title xyplothelper
#'
#' @description Simple helper for ploting xy dispersion points.
#'
#' @param d A data frame.
#' @param x A string with the dataframe column name for x axis
#' @param y A string with the dataframe column name for y axis
#' @param title The optional plot title. May be omited.
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
#' @param title The optional plot title. May be omited.
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
#' @description Simple helper function for countour plots
#'
#' @param d A data frame.
#' @param x A string with the dataframe column name for x axis.
#' @param y A string with the dataframe column name for y axis.
#' @param z A string with the dataframe column name for z axis.
#' @param nbins The number bins. The default is 32.
#' @param binwidth The binwidths for 'kde2d'. Can be an scalar or a vector.
#' @param points The number of grid points. Can be an scalar or a vector.
#' @param title The optional plot title. May be omited.
#'
#' @importFrom ggplot2 ggplot aes_string geom_density2d stat_density2d stat_contour
#' @importFrom ggplot2 ggtitle theme element_text scale_fill_gradient facet_wrap stat_summary2d
#' @importFrom ggplot2 stat_summary_2d scale_alpha facet_grid xlim ylim
#' @export
contourplothelper<- function(d, x, y, z, nbins=32, binwidth=c(10,10), points=c(300,300), title=NULL) {
  d<- as.data.frame(d)
  d<- fixdfcolumns(d,cols = c(x, y, z),skip = FALSE)

  ## --- Pruning duplicated entries
  d<- d[!duplicated(d[,c(x,y)]),]
  d<- d[!duplicated(d[,c(z)]),]

  if(is.null(title)) {
    d$title<- sprintf("%s ~ %s, %s", z, x, y)
  } else {
    d$title<- title
  }
  p<- with(d,ggplot(d, aes_string(x=x, y=y, z=z)))
  p<- p + stat_density2d(data = d, aes_string(fill = "..level..", alpha = "..level.."), size = 0.01, bins = nbins, geom = 'polygon', h = binwidth, n = points)
  p<- p + scale_fill_gradient(low = "green", high = "red")
  p<- p + scale_alpha(range = c(0.15, 0.30), guide = FALSE)
  p<- p + facet_grid(. ~ title) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  p
}

