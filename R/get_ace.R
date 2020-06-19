#' Get row-associated networks
#'
#' @return List of adjacency matrices
#'
#' @rdname rowNets
setMethod("rowNets", "ACTIONetExperiment", function(x) {
    out <- x@rowNets

    out
})


#' Get column-associated networks
#'
#' @return List of adjacency matrices
#'
#' @rdname colNets
setMethod("colNets", "ACTIONetExperiment", function(x) {
    out <- x@colNets

    out
})


#' Get row-associated factors
#'
#' @return List of matrices
#'
#' @rdname rowMaps
setMethod("rowMaps", "ACTIONetExperiment", function(x, all=TRUE) {
    out <- x@rowMaps
    if((all == FALSE) & ("rowMapsAnnot" %in% slotNames(x))) {
		mask = sapply(x@rowMapsAnnot, function(x) x$type != "internal")
		out = out[mask]
	}    
	
    out
})

#' Get column-associated factors
#'
#' @return List of matrices
#'
#' @rdname colMaps
setMethod("colMaps", "ACTIONetExperiment", function(x, all=TRUE) {
    out <- x@colMaps
    if((all == FALSE) & ("colMapsAnnot" %in% slotNames(x))) {
		mask = sapply(x@colMapsAnnot, function(x) x$type != "internal")
		out = out[mask]
	}

    out
})



GET_FUN <- function(exprs_values, ...) {
    (exprs_values) # To ensure evaluation
    function(object, ...) {
        assay(object, i=exprs_values, ...)
    }
}

SET_FUN <- function(exprs_values, ...) {
    (exprs_values) # To ensure evaluation
    function(object, ..., value) {
        assay(object, i=exprs_values, ...) <- value
        object
    }
}

#' @export
#' @importFrom BiocGenerics counts
setMethod("counts", "ACTIONetExperiment", GET_FUN("counts"))

#' @export
#' @importFrom BiocGenerics "counts<-"
setReplaceMethod("counts", c("ACTIONetExperiment", "ANY"), SET_FUN("counts"))

#' @export
setMethod("logcounts", "ACTIONetExperiment", GET_FUN("logcounts"))

#' @export
setReplaceMethod("logcounts", c("ACTIONetExperiment", "ANY"), SET_FUN("logcounts"))

#' @export
setMethod("normcounts", "ACTIONetExperiment", GET_FUN("normcounts"))

#' @export
setReplaceMethod("normcounts", c("ACTIONetExperiment", "ANY"), SET_FUN("normcounts"))
