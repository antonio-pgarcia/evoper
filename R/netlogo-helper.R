##================================================================================
## This file is part of the evoper package - EvoPER
##
## (C)2016, 2017, 2018 Antonio Prestes Garcia <@>
## For license terms see DESCRIPTION and/or LICENSE
##
## @file: netlogo-helper.R
##
## This file contains a simple integration wrapper for RNetLogo API
##================================================================================


#' @title NLWrapper.FindJar
#'
#' @description Search for the netlogo jar file on the provided path
#'
#' @param path The base path for searching
#'
#' @return The path for NetLogo jar file
#'
#' @export
NLWrapper.FindJar<- function(path) {
  list.files(path = path, pattern = "netlogo.*\\.jar")
}

#' @title NLWrapper.Model
#'
#' @description This wrapper prepares the environment and instantiates the model
#'
#' @param netlogodir The base path of NetLogo installation
#' @param modelfile The absolute path for NetLogo model file
#' @param maxtime The total number of iterations
#' @param dataset The names of model variables
#'
#' @examples \dontrun{
#'    rm(list=ls())
#'    p<- "C:/Program Files/NetLogo 6.1.1/app"
#'    output<- c("count sheep", "count wolves")
#'    m<- file.path(p, "models", "Sample Models", "Biology", "Wolf Sheep Predation.nlogo")
#'    o<- NLWrapper.Model(p, m, output, 150)
#' }
#'
#' @export
NLWrapper.Model<- function(netlogodir, modelfile, dataset, maxtime) {
  ppath<- getwd()
  setwd(netlogodir)
  Sys.setenv(NOAWT=1)
  netlogojar<- NLWrapper.FindJar(netlogodir)

  ## Instantiating
  obj<- RNetLogo::NLStart(netlogodir, gui=FALSE, nl.jarname=netlogojar)
  RNetLogo::NLLoadModel(modelfile, obj)
  v<- list(nl.obj=obj, nl.path=netlogodir, nl.model=modelfile, ppath=ppath, dataset=dataset, maxtime=maxtime )
  v
}

#' @title NLWrapper.Shutdown
#'
#' @description This wrapper terminates RNetLogo execution environment
#'
#' @param obj The object retuned by \link{NLWrapper.Model}
#'
#' @export
NLWrapper.Shutdown<- function(obj) {
  RNetLogo::NLQuit()
  setwd(obj$ppath)
}

#' @title NLWrapper.GetParameter
#'
#' @description Gets the value of a model parameter
#'
#' @param obj The object retuned by \link{NLWrapper.Model}
#' @param name The parameter name string or the collection of parameter names
#'
#' @return The parameter values
#'
#' @examples \dontrun{
#'    rm(list=ls())
#'    p<- "C:/Program Files/NetLogo 6.1.1/app"
#'    m<- file.path(nlpath, "models", "Sample Models", "Biology", "Wolf Sheep Predation.nlogo")
#'    o<- NLWrapper.Model(p, m)
#'    v<- NLWrapper.GetParameter(o, c("initial-number-sheep"))
#'
#'    or
#'
#'    v<- NLWrapper.GetParameter(o, c("initial-number-sheep","initial-number-wolves")))
#' }
#'
#' @export
NLWrapper.GetParameter<- function(obj, name) {
  v<- RNetLogo::NLReport(name, obj$nl.obj)
  v<- as.data.frame(v)
  names(v)<- name
  v
}

#' @title NLWrapper.SetParameter
#'
#' @description Set parameter values
#'
#' @param obj The object retuned by \link{NLWrapper.Model}
#' @param parameters The data frame containing the paramters
#'
#' @examples \dontrun{
#'    rm(list=ls())
#'    p<- "C:/Program Files/NetLogo 6.1.1/app"
#'    m<- file.path(nlpath, "models", "Sample Models", "Biology", "Wolf Sheep Predation.nlogo")
#'    o<- NLWrapper.Model(p, m)
#' }
#'
#' @export
NLWrapper.SetParameter<- function(obj, parameters) {
  for(name in colnames(parameters)) {
    RNetLogo::NLCommand(paste("set",name,parameters[1,name]))
  }
}

#' @title NLWrapper.SetRandomSeed
#'
#' @description Configures the random seed
#'
#' @param obj The object retuned by \link{NLWrapper.Model}
#' @param seed The new random seed
#'
#' @export
NLWrapper.SetRandomSeed<- function(obj, seed) {
  RNetLogo::NLCommand("random-seed", seed)
}

#' @title NLWrapper.Run
#'
#' @description Executes a NetLogo Model using rNetLogo
#'
#' @param obj The object retuned by \link{NLWrapper.Model}
#' @param r The number of replications
#' @param seed The collection of random seeds
#'
#' @examples \dontrun{
#' p<- "C:/Program Files/NetLogo 6.1.1/app"
#' m<- file.path(p, "models", "Sample Models", "Biology", "Wolf Sheep Predation.nlogo")
#' output<- c("count sheep", "count wolves")
#' o<- NLWrapper.Model(p, m, output, 150)
#' v<- NLWrapper.Run(o)
#' NLWrapper.Shutdown(o)}
#'
#'
#' @export
NLWrapper.Run<- function(obj, r=1, seed=c()) {
  if(length(seed) == 0) {
    seed= runif(r,-10^8,10^8)
  } else if(length(seed) != r) {
    stop("The provided set of random numbers doesn't match replications!")
  }

  output<- c()
  myreporter<- c("ticks",obj$dataset)
  ## --- Execute the replications
  for(run in 1:r) {
    random_seed<- seed[run]
    NLWrapper.SetRandomSeed(obj, seed[run])
    RNetLogo::NLCommand("setup")
    v<- RNetLogo::NLDoReport(iterations= obj$maxtime, command= "go", reporter= myreporter, as.data.frame= TRUE, df.col.names= c())
    names(v)<- c(c("time"),obj$dataset)
    output<- rbind(output,cbind(random_seed,run,v))
  }
  output
}

#' @title NLWrapper.RunExperiment
#'
#' @description Executes a NetLogo Model using rNetLogo
#'
#' @param obj The object retuned by \link{NLWrapper.Model}
#' @param r The number of replications
#' @param design The desing matrix holding parameter sampling
#' @param FUN THe calibration function.
#'
#' @return A list containing the the parameters, the calibration functio output and the whole resultset
#'
#' @examples \dontrun{
#'    rm(list=ls())
#'    objectivefn<- function(params, results) { 0 }
#'
#'    f<- AddFactor(name="initial-number-sheep",min=100,max=250)
#'    f<- AddFactor(factors=f, name="initial-number-wolves",min=50,max=150)
#'    f<- AddFactor(factors=f, name="grass-regrowth-time",min=30,max=100)
#'    f<- AddFactor(factors=f, name="sheep-gain-from-food",min=1,max=50)
#'    f<- AddFactor(factors=f, name="wolf-gain-from-food",min=1,max=100)
#'    f<- AddFactor(factors=f, name="sheep-reproduce",min=1,max=20)
#'    f<- AddFactor(factors=f, name="wolf-reproduce",min=1,max=20)
#'
#'    design<- AoE.LatinHypercube(factors=f)
#'
#'    p<- "C:/Program Files/NetLogo 6.1.1/app"
#'    m<- file.path(p, "models", "Sample Models", "Biology", "Wolf Sheep Predation.nlogo")
#'    output<- c("count sheep", "count wolves")
#'    o<- NLWrapper.Model(p, m, output, 150)
#'    v<- RunExperiment(o, r=1, design, objectivefn)
#'    NLWrapper.Shutdown(o)
#' }
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
NLWrapper.RunExperiment<- function(obj, r=1, design, FUN) {
  paramset<- c()
  output<- c()
  dataset<- c()

  psets<- nrow(design)
  myactivity<- txtProgressBar(min = 0, max = psets, style = 3)

  for(pset in 1:psets) {
    d<- design[pset,]

    NLWrapper.SetParameter(obj, d)

    # -- Run model with current parameter set
    results<- evoper::NLWrapper.Run(obj, r)

    # -- The user provided calibration function.
    # -- Calibration function must return 0 for perfect fit between
    # -- observed and experimental data.
    calibration<- FUN(d, results)

    if(is.null(calibration)) {
      stop("Invalid user provided calibration function!")
    }

    setTxtProgressBar(myactivity, pset)

    paramset<- rbind(paramset,cbind(pset,d))
    output<- rbind(output,cbind(pset,calibration))
    dataset<- rbind(dataset,cbind(pset,results))
    #paramset<- rbindlist(list(paramset, c(pset,d)), use.names=TRUE)
    #output<- rbindlist(list(output, c(pset,setNames(calibration,"calibration"))), use.names=TRUE)
    #dataset<- rbindlist(list(dataset, c(pset,results)), use.names=TRUE)
  }
  close(myactivity)
  return(list(paramset=paramset, output=output, dataset=dataset))
}

