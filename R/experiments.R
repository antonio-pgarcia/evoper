##================================================================================
## This file is part of the evoper package - EvoPER
##
## (C)2016 Antonio Prestes Garcia <@>
## For license terms see DESCRIPTION and/or LICENSE
##
## $Id$
##================================================================================

#' @title compare.algorithms1
#'
#' @description Compare the number of function evalutions and convergence for the
#' following optimization algorithms, ("saa","pso","acor","ees1").
#'
#' @param F The function to be tested
#' @param seeds The random seeds which will be used for testing algorithms
#'
#' @examples \dontrun{
#'  rm(list=ls())
#'  d.cigar4<- compare.algorithms1(f0.cigar4)
#'  d.schaffer4<- compare.algorithms1(f0.schaffer4)
#'  d.griewank4<- compare.algorithms1(f0.griewank4)
#'  d.bohachevsky4<- compare.algorithms1(f0.bohachevsky4)
#'  d.rosenbrock4<- compare.algorithms1(f0.rosenbrock4)
#' }
#'
#' @export
compare.algorithms1<- function(F, seeds= c(27, 2718282, 36190727, 3141593, -91190721, -140743, 1321)) {
  algorithms<- c("saa","pso","acor","ees1")

  mydata<- c()
  for(algorithm in algorithms) {
    for(seed in seeds) {
      set.seed(seed)
      f<- PlainFunction$new(F)
      f$setTolerance(10^-1)
      f$Parameter(name="x1",min=-100,max=100)
      f$Parameter(name="x2",min=-100,max=100)
      f$Parameter(name="x3",min=-100,max=100)
      f$Parameter(name="x4",min=-100,max=100)
      v<- extremize(algorithm, f)

      myrow<- cbind(algorithm, seed, f$stats(), v$getBest())
      mydata<- rbind(mydata, myrow)
    }
  }
  as.data.frame(mydata)
}

#' @title summarize.comp1
#'
#' @description Provides as summary with averged values of experimental setup
#'
#' @param mydata The data frame generated with 'compare.algorithms1'
#'
#' @return The summarized data
#'
#' @import plyr
#' @export
summarize.comp1<- function(mydata) {
  ddply(mydata, .(algorithm), summarize,  evals=mean(total_evals), convergence=mean(converged), fitness=mean(fitness))
}

#' @title show.comp1
#'
#' @description Generates a barplot comparing the number of evalutions for
#' algorithms ("saa","pso","acor","ees1").
#'
#' @param mydata The data generated with 'summarize.comp1'
#' @param what The name of variable to plot on 'y' axis
#' @param title the plot title
#'
#' @examples \dontrun{
#' 	p.a<- show.comp1(d.cigar4,"evals","(a) Cigar function")
#' 	p.b<- show.comp1(d.schaffer4,"evals","(b) Schafer function")
#' 	p.c<- show.comp1(d.griewank4,"evals","(c) Griewank function")
#' 	p.d<- show.comp1(d.bohachevsky4,"evals","(d) Bohachevsky function")
#' }
#'
#' @importFrom ggplot2 geom_bar
#' @export
show.comp1<- function(mydata, what, title=NULL) {
  mydata$title<- sprintf("%s", title)
  p<- ggplot(data= mydata, with( mydata, aes_string(x="algorithm", y=what)) )
  p<- p + geom_bar(stat="identity", fill=ifelse(mydata$convergence < 0.6,"gray", "steelblue"))
  p<- p + facet_grid(. ~ title)
  p
}

