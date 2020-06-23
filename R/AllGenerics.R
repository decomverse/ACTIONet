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




#' Get row embeddings
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowEmbeddings
setGeneric("rowEmbeddings", function(x, ...) standardGeneric("rowEmbeddings"))



#' Get column embeddings
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colEmbeddings
setGeneric("colEmbeddings", function(x, ...) standardGeneric("colEmbeddings"))

#' Set row embeddings
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowEmbeddings
setGeneric("rowEmbeddings<-", function(x, ...) standardGeneric("rowEmbeddings<-"))



#' Set column embeddings
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colEmbeddings
setGeneric("colEmbeddings<-", function(x, ...) standardGeneric("colEmbeddings<-"))


#' Get row reductions
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowReductions
setGeneric("rowReductions", function(x, ...) standardGeneric("rowReductions"))



#' Get column reductions
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colReductions
setGeneric("colReductions", function(x, ...) standardGeneric("colReductions"))

#' Set column reductions
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod colReductions<-
setGeneric("colReductions<-", function(x, ...) standardGeneric("colReductions<-"))

#' Set row reductions
#'
#' @param x ACTIONetExperiment object
#' @param ... other parameters
#'
#' @exportMethod rowReductions<-
setGeneric("rowReductions<-", function(x, ...) standardGeneric("rowReductions<-"))
