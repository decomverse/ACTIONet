#' Get rowNets
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowNets
setGeneric("rowNets", function(object, ...) standardGeneric("rowNets"))

#' Get colNets
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colNets
#' @export
setGeneric("colNets", function(object, ...) standardGeneric("colNets"))

#' Get rowMaps
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowMaps
#' @export
setGeneric("rowMaps", function(object, ...) standardGeneric("rowMaps"))

#' Get colMaps
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colMaps
#' @export
setGeneric("colMaps", function(object, ...) standardGeneric("colMaps"))

#' Set rowNets
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowNets<-
setGeneric("rowNets<-", function(object, ...) standardGeneric("rowNets<-"))

#' Set colNets
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colNets<-
setGeneric("colNets<-", function(object, ...) standardGeneric("colNets<-"))

#' Set rowMaps
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowMaps<-
setGeneric("rowMaps<-", function(object, ...) standardGeneric("rowMaps<-"))

#' Set colMaps
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colMaps<-
setGeneric("colMaps<-", function(object, ...) standardGeneric("colMaps<-"))

########################### Type functions ################################
#' Get rowMapTypes
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowMapTypes
setGeneric("rowMapTypes", function(object, ...) standardGeneric("rowMapTypes"))

#' Get colMapTypes
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colMapTypes
setGeneric("colMapTypes", function(object, ...) standardGeneric("colMapTypes"))


#' Set rowMapTypes
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowMapTypes<-
setGeneric("rowMapTypes<-", function(object, ...) standardGeneric("rowMapTypes<-"))

#' Set colMapTypes
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colMapTypes<-
setGeneric("colMapTypes<-", function(object, ...) standardGeneric("colMapTypes<-"))

########################### Metadata functions ##############################

#' Get rowMapMeta
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowMapMeta
setGeneric("rowMapMeta", function(object, ...) standardGeneric("rowMapMeta"))

#' Get colMapTypes
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colMapTypes
setGeneric("colMapMeta", function(object, ...) standardGeneric("colMapMeta"))


#' Set rowMapMeta
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowMapMeta<-
setGeneric("rowMapMeta<-", function(object, ...) standardGeneric("rowMapMeta<-"))

#' Set colMapMeta
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colMapMeta<-
setGeneric("colMapMeta<-", function(object, ...) standardGeneric("colMapMeta<-"))

#' Get row embeddings
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowEmbeddings
setGeneric("rowEmbeddings", function(object, ...) standardGeneric("rowEmbeddings"))

#' Get column embeddings
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colEmbeddings
setGeneric("colEmbeddings", function(object, ...) standardGeneric("colEmbeddings"))

#' Set row embeddings
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowEmbeddings
setGeneric("rowEmbeddings<-", function(object, ...) standardGeneric("rowEmbeddings<-"))

#' Set column embeddings
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colEmbeddings
setGeneric("colEmbeddings<-", function(object, ...) standardGeneric("colEmbeddings<-"))


#' Get row reductions
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowReductions
setGeneric("rowReductions", function(object, ...) standardGeneric("rowReductions"))

#' Get column reductions
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colReductions
setGeneric("colReductions", function(object, ...) standardGeneric("colReductions"))

#' Set column reductions
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colReductions<-
setGeneric("colReductions<-", function(object, ...) standardGeneric("colReductions<-"))

#' Set row reductions
#'
#' @param object ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowReductions<-
setGeneric("rowReductions<-", function(object, ..., value) standardGeneric("rowReductions<-"))

#' @export
setGeneric("logcounts", function(object, ...) standardGeneric("logcounts"))

#' @export
setGeneric("logcounts<-", function(object, ...) standardGeneric("logcounts<-"))

#' @export
setGeneric("normcounts", function(object, ...) standardGeneric("normcounts"))

#' @export
setGeneric("normcounts<-", function(object, ..., value) standardGeneric("normcounts<-"))

#' @export
setGeneric("reducedDimNames", function(object) standardGeneric("reducedDimNames"))

#' @export
setGeneric("reducedDimNames<-", function(object, value) standardGeneric("reducedDimNames<-"))

#' @export
setGeneric("reducedDims", function(object, ...) standardGeneric("reducedDims"))

#' @export
setGeneric("reducedDims<-", function(object, value) standardGeneric("reducedDims<-"))

#' @export
setGeneric("rownames<-", function(object, ..., value) standardGeneric("rownames<-"))

#' @export
setGeneric("colnames<-", function(object, ..., value) standardGeneric("colnames<-"))
