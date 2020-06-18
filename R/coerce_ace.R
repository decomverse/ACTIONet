#' Coerces a RangedSummarizedExperiment method as an ACTIONetExperiment (ACE) object
#'
#' @param from RangedSummarizedExperiment object
#'
#' @exportMethod coerce
setAs("RangedSummarizedExperiment", "ACTIONetExperiment", function(from) {
    new("ACTIONetExperiment", from, 
    rowNets=S4Vectors::SimpleList(), 
    colNets=S4Vectors::SimpleList(), 
    rowFactors=S4Vectors::SimpleList(), 
    colFactors=S4Vectors::SimpleList())
})

setAs("ACTIONetExperiment", "SingleCellExperiment", function(from) {	
    new("SingleCellExperiment", as(from, "RangedSummarizedExperiment"), 
    reducedDims=colFactors(from))
})

setAs("SingleCellExperiment", "ACTIONetExperiment", function(from) {	
    new("ACTIONetExperiment", as(from, "RangedSummarizedExperiment"), 
    colFactors=reducedDims(from))
})

