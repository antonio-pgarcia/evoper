% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test-functions.R
\name{f0.periodtuningpp48}
\alias{f0.periodtuningpp48}
\title{Period tuning of 48 time units for Predator-Prey}
\usage{
f0.periodtuningpp48(x1, x2, x3, x4)
}
\arguments{
\item{x1}{The growth rate of prey}

\item{x2}{The decay rate of predator}

\item{x3}{The predating effect on prey}

\item{x4}{The predating effecto on predator}
}
\value{
The solution fitness cost
}
\description{
This function is an example on how EvoPER can be
used for estimating the parameter values in order to produce
oscilations with the desired period.
}
\examples{
\dontrun{
rm(list=ls())
set.seed(-27262565)
f<- PlainFunction$new(f0.periodtuningpp24)
f$Parameter(name="x1",min=0.5,max=2)
f$Parameter(name="x2",min=0.5,max=2)
f$Parameter(name="x3",min=0.5,max=2)
f$Parameter(name="x4",min=0.5,max=2)
extremize("pso", f)
}

}
