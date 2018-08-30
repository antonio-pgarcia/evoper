##================================================================================
## This file is part of EvoPER 
##
## Example of parameter estimation for NetLogo models
## 
## (C) 2016, 2017, 2018 - Antonio Prestes Garcia <@>
## For license terms see DESCRIPTION and/or LICENSE
##
## @file: netlogo-predator-prey.R
##
## This file contains R script for the parameter estimation of predator prey
## NetLogo model
##================================================================================

## Initialization
rm(list=ls())
library(evoper)
set.seed(2718282)

netlogo<- "C:/Program Files/NetLogo 6.0.4/app"
model<- file.path(netlogo, "models", "Sample Models", "Biology", "Wolf Sheep Predation.nlogo")
output<- c("count sheep", "count wolves")

default<- c(`initial-number-sheep`=100, `initial-number-wolves`=50)

## Objective function
my.objectivefn<- function(params, results) {
  	p.sheep<- naiveperiod(results$`count sheep`)
  	p.wolves<- naiveperiod(results$`count wolves`)
	dsheep<- rrepast::AoE.NRMSD(p.sheep$period, 30)
	dwolves<- rrepast::AoE.NRMSD(p.sheep$period, 30)   
   	criteria<- cbind(dsheep,dwolves)
   	return(criteria)
}


## Objective function instantiation
objective<- NetLogoFunction$new(netlogo, model, output, 100, my.objectivefn)

## Factors under study
objective$Parameter(name="grass-regrowth-time",min=30,max=100)
objective$Parameter(name="sheep-gain-from-food",min=1,max=25)
objective$Parameter(name="wolf-gain-from-food",min=1,max=50)
objective$Parameter(name="sheep-reproduce",min=1,max=10)
objective$Parameter(name="wolf-reproduce",min=1,max=10)

## Setting tolerance and default parameters
objective$setTolerance(0.1)
objective$defaults(default)

## Optimization method
results<- extremize("ees2", objective)

