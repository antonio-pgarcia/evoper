##================================================================================
## This file is part of the evoper package - EvoPER
##
## (C)2016, 2017 Antonio Prestes Garcia <@>
## For license terms see DESCRIPTION and/or LICENSE
##
## @file: class-options.R
##
## This file contains the classes abstracting the options for metaheuristics
##================================================================================


#' @title Options
#'
#' @description The base class for the options for the optimization metaheuristics
#'
#' @field type The configuration type
#' @field neighborhood The neighborhood function for population methods
#' @field discrete Flag indicating that and specific algorithm is discrete or continuous
#' @field nlevelz Default value for generating parameter levels when range is provided, default value is 5
#' @field container The object holding the configuration otions
#'
#' @importFrom methods new
#' @export Options
#' @exportClass Options
Options<- setRefClass("Options",

  fields = list(
    type = 'ANY',
    discrete = 'ANY',
    nlevelz = 'ANY',
    neighborhood = 'ANY',
    container = 'list'
  ),


  methods = list(

    initialize = function() {
      type<<- 'none'
      discrete<<- FALSE
      nlevelz<<- 5
      container<<- list(iterations=500, trace=FALSE)
    },

    ## Set/Get the neighborhood function
    neighborhoodFunction = function(f=NULL) {
      if(!is.null(f)) {
        neighborhood<<- f
      }
      neighborhood
    },


    setType = function(v) {
      type<<- v
    },

    getType = function() {
      type
    },

    setDiscrete = function(v) {
      discrete<<- v
    },

    isDiscrete = function() {
      discrete
    },

    setLevels = function(v) {
      nlevelz<<- v
    },

    getLevels = function() {
      nlevelz
    },

    setValue = function(k , v) {
      container[k]<<- v
    },

    getValue = function(k) {
      (container[k])[[1]]
    },

    Keys = function() {
      names(container)
    },

    toString = function() {
      sstring<- c()
      values<- unlist(container)
      for(k in names(values)) { sstring<- paste0(sstring,k,"=",values[k],sep='\t') }
      sstring
    }

  )
)


#' @title OptionsPSO
#'
#' @description Options for PSO optimization metaheuristic
#'
#' @importFrom methods new
#' @export OptionsPSO
#' @exportClass OptionsPSO
OptionsPSO<- setRefClass("OptionsPSO", contains = "Options",

  methods = list(

    initialize = function() {

      callSuper()
      setType("pso")
      setValue("N",16)
      setValue("phi1",1.193)
      setValue("phi2",1.193)
      setValue("W",0.721)
      neighborhoodFunction(pso.neighborhood.K4)
    }

  )
)

#' @title OptionsSAA
#'
#' @description Options for SAA method
#'
#' @field temperature The temperature dacay function
#'
#' @importFrom methods new
#' @export OptionsSAA
#' @exportClass OptionsSAA
OptionsSAA<- setRefClass("OptionsSAA", contains = "Options",

  fields = list(
    temperature = 'ANY'
  ),


  methods = list(

    initialize = function() {

      callSuper()
      setType("saa")
      setValue("iterations", 500)
      setValue("t0", 8*10^12)
      setValue("t.min", 10^-5)
      setValue("L", 10)
      setValue("d", 0.5)
      neighborhoodFunction(saa.neighborhoodN)
      setTemperatureF(saa.tcte)
    },

    setTemperatureF = function(v) {
        temperature<<- v
    },

    getTemperatureF = function() {
      temperature
    }

  )
)


#' @title OptionsACOR
#'
#' @description Options for ACOR method
#'
#' @importFrom methods new
#' @export OptionsACOR
#' @exportClass OptionsACOR
OptionsACOR<- setRefClass("OptionsACOR", contains = "Options",

  methods = list(

    initialize = function() {

      callSuper()
      setType("acor")
      setValue("n.ants", 64)  ## The number of simulated ants
      setValue("k", 32)       ## The archive size
      setValue("q", 0.2)      ## Locality of the search process
      setValue("Xi", 0.5)     ## Equivalent to evaporation rate, higher Xi reduce convergence speed
    }

  )
)


#' @title OptionsEES1
#'
#' @description Options for EvoPER Evolutionary Stratety 1
#'
#'
#' @importFrom methods new
#' @export OptionsEES1
#' @exportClass OptionsEES1
OptionsEES1<- setRefClass("OptionsEES1", contains = "Options",

  fields = list(
  ),


  methods = list(

    initialize = function() {

      callSuper()
      setType("ees1")
      setValue("N", 10)             ## Solution size
      setValue("mu", 0.3)           ## Fitness preference strenght
      setValue("rho", 0.01)         ## Mutation probability
      setValue("kkappa", 0.2)       ## Selective pressure
      setValue("iterations", 100)   ## Total number of iterations
    }

  )
)


#' @title OptionsEES2
#'
#' @description Options for Serial Dilutions method
#'
#' @field dilutions The desired dilutions
#'
#' @importFrom methods new
#' @export OptionsEES2
#' @exportClass OptionsEES2
OptionsEES2<- setRefClass("OptionsEES2", contains = "Options",

  fields = list(
  ),


  methods = list(

    initialize = function() {
      callSuper()
      setType("ees2")
      setValue("N", 20)            ## Solution size 100
      setValue("rho", 0.25)        ## Solution size 0.05
      setValue("iterations", 30)   ## Total number of iterations
    }

  )
)

#' @title OptionsTS
#'
#' @description Options for Tabu search optimization metaheuristic
#'
#' @importFrom methods new
#' @export OptionsTS
#' @exportClass OptionsTS
OptionsTS<- setRefClass("OptionsTS", contains = "Options",

  methods = list(

    initialize = function() {

      callSuper()
      setType("tabu")
      setDiscrete(TRUE)
      setValue("N", 25)             ## Neighbor solution size
      setValue("tabu_size", 16)     ## Tabu size
      setValue("iterations", 400)   ## Total number of iterations
    }

  )
)

#' @title OptionsGA
#'
#' @description Options for Genetic Algorithm optimization metaheuristic
#'
#' @importFrom methods new
#' @export OptionsGA
#' @exportClass OptionsGA
OptionsGA<- setRefClass("OptionsGA", contains = "Options",

  methods = list(

    initialize = function() {

      callSuper()
      setType("ga")
      setDiscrete(TRUE)
      setValue("N", 25)             ## Population size
      setValue("iterations", 400)   ## Total number of iterations
      setValue("tournament_size", 10)
      setValue("crossover_rate", 0.60)
      setValue("mutation_rate", 0.10)
    }

  )
)

#' @title OptionsFactory
#'
#' @description Instantiate the Options class required for
#' the specific metaheuristic method.
#'
#' @param type The metaheuristic method
#' @param v The options object
#'
#' @return Options object
#'
#' @export
OptionsFactory<- function(type, v=NULL) {

  if(is.null(v)) {
    switch(type,
      pso = { v<- OptionsPSO$new() },
      saa = { v<- OptionsSAA$new() },
      acor= { v<- OptionsACOR$new() },
      ees1= { v<- OptionsEES1$new() },
      ees2= { v<- OptionsEES2$new() },
      tabu= { v<- OptionsTS$new() },
      ga  = { v<- OptionsGA$new() },
      { stop("Invalid optimization function!") }
    )
  }

  if(v$getType() != type) stop(paste("Invalid option of type [", v$getType(),"]"))
  v
}
