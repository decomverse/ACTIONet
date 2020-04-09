#' Get row-associated networks
#'
#' @return List of adjacency matrices
#'
#' @rdname rowNets
setMethod("rowNets", "ACTIONetExperiment", function(x, withDimnames=TRUE) {
    out <- x@rowNets
    if (withDimnames & (length(out) > 0)) {
        for (i in 1:length(out)) {
            rownames(out[[i]]) <- colnames(out[[i]]) <- rownames(x)
        }                        
	}
    out
})


#' Get column-associated networks
#'
#' @return List of adjacency matrices
#'
#' @rdname colNets
setMethod("colNets", "ACTIONetExperiment", function(x, withDimnames=TRUE) {
    out <- x@colNets
    if (withDimnames & (length(out) > 0)) {
        for (i in 1:length(out)) {
            rownames(out[[i]]) <- colnames(out[[i]]) <- colnames(x)
        }                        
	}
    out
})


#' Get row-associated factors
#'
#' @return List of matrices
#'
#' @rdname rowFactors
setMethod("rowFactors", "ACTIONetExperiment", function(x, withDimnames=TRUE) {
    out <- x@rowFactors
    if (withDimnames & (length(out) > 0)) {
        for (i in 1:length(out)) {
            rownames(out[[i]]) <- rownames(x)
        }                        
	}
    out
})

#' Get column-associated factors
#'
#' @return List of matrices
#'
#' @rdname colFactors
setMethod("colFactors", "ACTIONetExperiment", function(x, withDimnames=TRUE) {
    out <- x@colFactors
    if (withDimnames & (length(out) > 0)) {
        for (i in 1:length(out)) {
            colnames(out[[i]]) <- colnames(x)
        }                        
	}
    out
})
