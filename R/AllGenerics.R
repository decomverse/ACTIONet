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

#' Get rowMaps
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowMaps
#' @export
setGeneric("rowMaps", function(x, ...) standardGeneric("rowMaps"))

#' Get colMaps
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colMaps
#' @export
setGeneric("colMaps", function(x, ...) standardGeneric("colMaps"))

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

#' Set rowMaps
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowMaps<-
setGeneric("rowMaps<-", function(x, ...) standardGeneric("rowMaps<-"))

#' Set colMaps
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colMaps<-
setGeneric("colMaps<-", function(x, ...) standardGeneric("colMaps<-"))



########################### Type functions ################################
#' Get rowMapTypes
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowMapTypes
setGeneric("rowMapTypes", function(x, ...) standardGeneric("rowMapTypes"))

#' Get colMapTypes
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colMapTypes
setGeneric("colMapTypes", function(x, ...) standardGeneric("colMapTypes"))


#' Set rowMapTypes
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowMapTypes<-
setGeneric("rowMapTypes<-", function(x, ...) standardGeneric("rowMapTypes<-"))

#' Set colMapTypes
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colMapTypes<-
setGeneric("colMapTypes<-", function(x, ...) standardGeneric("colMapTypes<-"))

########################### Metadata functions ##############################

#' Get rowMapMeta
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowMapMeta
setGeneric("rowMapMeta", function(x, ...) standardGeneric("rowMapMeta"))

#' Get colMapTypes
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colMapTypes
setGeneric("colMapMeta", function(x, ...) standardGeneric("colMapMeta"))


#' Set rowMapMeta
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowMapMeta<-
setGeneric("rowMapMeta<-", function(x, ...) standardGeneric("rowMapMeta<-"))

#' Set colMapMeta
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colMapMeta<-
setGeneric("colMapMeta<-", function(x, ...) standardGeneric("colMapMeta<-"))


########################### For compatibility ##############################
#' Get counts assay
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod counts
setGeneric("counts", function(x, ...) standardGeneric("counts"))



#' Set counts assay
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod counts<-
setGeneric("counts<-", function(x, ...) standardGeneric("counts<-"))

#' Get logcounts assay
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod logcounts
setGeneric("logcounts", function(x, ...) standardGeneric("logcounts"))



#' Set logcounts assay
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod logcounts<-
setGeneric("logcounts<-", function(x, ...) standardGeneric("logcounts<-"))



#' Get normcounts assay
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod normcounts
setGeneric("normcounts", function(x, ...) standardGeneric("normcounts"))



#' Set normcounts assay
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod normcounts<-
setGeneric("normcounts<-", function(x, ...) standardGeneric("normcounts<-"))




#' Get colMaps() as reducedDims()
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod reducedDims
setGeneric("reducedDims", function(x, ...) standardGeneric("reducedDims"))




#' Get colMaps() as reducedDims()
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod reducedDims
setGeneric("reducedDims", function(x, ...) standardGeneric("reducedDims"))


#' Get names(colMaps()) as reducedDimNames()
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod reducedDimNames
setGeneric("reducedDimNames", function(x, ...) standardGeneric("reducedDimNames"))
