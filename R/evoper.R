##================================================================================
## This file is part of the evoper package - EvoPER
##
## (C)2016 Antonio Prestes Garcia <@>
## For license terms see DESCRIPTION and/or LICENSE
##
## $Id$
##================================================================================


# ------------------------------------------------------------
# .onLoad, Hook for loading package namespace
#
# ------------------------------------------------------------
.onLoad<- function(libname, pkgname) {
  assign("pkg.globals", new.env(), envir=parent.env(environment()))

  # The default Random Seed. You may want to change this.
  set.seed(exp(1)*10^6)

  # Internal variables
  assign("pkg.basedir", NA, pkg.globals)

  # Set default package log log leve
  assign("pkg.LogLevel", 'ERROR', pkg.globals)
}

#' @title assert
#'
#' @description The assert function stop the execution if the logical
#' expression given by the parameter \code{expresion} is false.
#'
#' @param expresion Some logical expression
#' @param string The text message to show if expression does not hold
#'
#' @export
assert<- function(expresion, string) {
  if(!expresion) stop(string, call. = FALSE)
}


## ##################################################################
##
## ----- Wrapping the logging system
##
## ##################################################################


#' @title GetLogLevel
#'
#' @description Get the current log level
#'
#' @return The log level
#'
#' @export
GetLogLevel<- function() {
  return(get("pkg.LogLevel", pkg.globals))
}

#' @title SetLogLevel
#'
#' @description Configure the current log level for the package
#'
#' @param level The log level (ERROR|WARN|INFO|DEBUG)
#'
#' @export
SetLogLevel<- function(level) {
  assign("pkg.LogLevel", level, pkg.globals)
}

#' @title elog.level
#'
#' @description Configure the current log level
#'
#' @param level The log level (ERROR|WARN|INFO|DEBUG)
#'
#' @return The log level
#'
#' @importFrom logging basicConfig
#' @export
elog.level<- function(level=NULL) {
  if(!is.null(level)) {
    #flog.threshold(level)
    SetLogLevel(level)
    basicConfig(level = level)
  }
  GetLogLevel()
}

#' @title elog.error
#'
#' @description Wrapper for logging error messages.
#'
#' @param ... Variable number of arguments including a format string.
#'
#' @importFrom logging logerror
#' @export
elog.error<- function(...) {
  #flog.error(...)
  logerror(...)
}

#' @title elog.info
#'
#' @description Wrapper for logging info messages.
#'
#' @param ... Variable number of arguments including a format string.
#'
#' @importFrom logging loginfo
#' @export
elog.info<- function(...) {
  #log.info(...)
  loginfo(...)
}

#' @title elog.debug
#'
#' @description Wrapper for logging debug messages.
#'
#' @param ... Variable number of arguments including a format string.
#'
#' @importFrom logging logdebug
#' @export
elog.debug<- function(...) {
  #log.debug(...)
  logdebug(...)
}



## ##################################################################
##
## --- Entry point function -
##
## ##################################################################


#' @title extremize
#'
#' @description Entry point for optimization functions
#'
#' @param type The optimization method (aco,pso,saa,sda)
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
#' @export
extremize<- function(type, objective, options = NULL) {

  switch(type,
    pso={
      optimization.fun<- abm.pso
    },

    saa={
      optimization.fun<- abm.saa
    },

    acor={
      optimization.fun<- abm.acor
    },

    ees1={
      optimization.fun<- abm.ees1
    },

    ees2={
      optimization.fun<- abm.ees2
    },

    tabu={
      optimization.fun<- abm.tabu
    },

    ga={
      # Temporal for fixing futile.logging
      optimization.fun<- abm.tabu
    },

    {
      stop("Invalid optimization function!")
    }
  )

  optimization.fun(objective, options)
}


