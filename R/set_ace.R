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
#' @rdname rowFactors
setReplaceMethod("rowFactors", "ACTIONetExperiment", function(x, value) {
    x@rowFactors <- value
    validObject(x)
    x
})

#' Set column-associated factors
#'
#' @return List of matrices
#'
#' @rdname colFactors
setReplaceMethod("colFactors", "ACTIONetExperiment", function(x, value) {
    x@colFactors <- value
    validObject(x)
    x
})
