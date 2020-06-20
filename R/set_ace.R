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
  x <- .check_new_map(x)
	x@rowMaps <- value
  validObject(x)
  x
})


#' Set column-associated factors
#'
#' @return List of matrices
#'
#' @rdname colMaps
setReplaceMethod("colMaps", "ACTIONetExperiment", function(x, value) {
  x <- .check_new_map(x)
	x@colMaps <- value
  validObject(x)
	x
})

.check_new_map <- function(x){
  err = sprintf("Mapping must be a named list.\n")
  if(!(class(x) %in% c("list", "SimpleList")))
    stop(err)
  if(is.null(names(x)))
    stop(err)

  x = as(x, "SimpleList")
  return(x)
}
