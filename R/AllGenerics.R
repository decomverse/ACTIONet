#' Get rowNets
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowNets
setGeneric("rowNets", function(x, ...) standardGeneric("rowNets"))

#' Get colNets
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colNets
#' @export
setGeneric("colNets", function(x, ...) standardGeneric("colNets"))

#' Get rowFactors
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowFactors
#' @export
setGeneric("rowFactors", function(x, ...) standardGeneric("rowFactors"))

#' Get colFactors
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colFactors
#' @export
setGeneric("colFactors", function(x, ...) standardGeneric("colFactors"))

#' Set rowNets
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowNets<-
setGeneric("rowNets<-", function(x, ...) standardGeneric("rowNets<-"))

#' Set colNets
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colNets<-
setGeneric("colNets<-", function(x, ...) standardGeneric("colNets<-"))

#' Set rowFactors
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowFactors<-
setGeneric("rowFactors<-", function(x, ...) standardGeneric("rowFactors<-"))

#' Set colFactors
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colFactors<-
setGeneric("colFactors<-", function(x, ...) standardGeneric("colFactors<-"))



