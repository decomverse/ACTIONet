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
	sce = as(as(from, "RangedSummarizedExperiment"), "SingleCellExperiment")
    #sce = as(from, "SingleCellExperiment")
    reducedDims(sce)=SimpleList(lapply(colFactors(from), function(x) Matrix::t(x)))
    return(sce)
})

setAs("SingleCellExperiment", "ACTIONetExperiment", function(from) {	
	ace = as(as(from, "RangedSummarizedExperiment"), "ACTIONetExperiment")
    #ace = as(from, "ACTIONetExperiment")
    
    colFactors(ace)=SimpleList(lapply(reducedDims(from), function(x) Matrix::t(x)))
    return(ace)
})

