#' Coerces a SingleCellExperiment method as an ACTIONetExperiment (ACE) object
#'
#' @param from SingleCellExperiment object
#'
#' @examples
#' @exportMethods coerce
setAs("SingleCellExperiment", "ACTIONetExperiment", function(from) {
    new("ACTIONetExperiment", from, 
    rowNets=S4Vectors::SimpleList(), 
    colNets=S4Vectors::SimpleList(), 
    rowFactors=S4Vectors::SimpleList(), 
    colFactors=S4Vectors::SimpleList())
})
