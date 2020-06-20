#' Set row-associated networks
#'
#' @return List of adjacency matrices
#'
#' @rdname rowNets
setReplaceMethod("rowNets", "ACTIONetExperiment", function(x, value) {
    x@rowNets <- value
    validObject(x)
    x
})


#' Set column-associated networks
#'
#' @return List of adjacency matrices
#'
#' @rdname colNets
setReplaceMethod("colNets", "ACTIONetExperiment", function(x, value) {
    x@colNets <- value
    validObject(x)
    x
})


#' Set row-associated factors
#'
#' @return List of matrices
#'
#' @rdname rowMaps
setReplaceMethod("rowMaps", "ACTIONetExperiment", function(x, value) {
	SEs = lapply(value, function(X) {
		SE = SummarizedExperiment(assays=list(X=X))
		metadata(SE)$type = "internal"
		return(SE)
	}
	x@rowMaps <- SEs

    validObject(x)
    x
})


#' Set column-associated factors
#'
#' @return List of matrices
#'
#' @rdname colMaps
setReplaceMethod("colMaps", "ACTIONetExperiment", function(x, value) {
	x@colMaps <- value
    validObject(x)
	x
})
