#' Coerces a SummarizedExperiment method as an ACTIONetExperiment (ACE) object
#'
#' @param from SummarizedExperiment object
#'
#' @exportMethod coerce
setAs("SummarizedExperiment", "ACTIONetExperiment", function(from) {
    new("ACTIONetExperiment", from, 
    rowNets=S4Vectors::SimpleList(), 
    colNets=S4Vectors::SimpleList(), 
    rowFactors=S4Vectors::SimpleList(), 
    colFactors=S4Vectors::SimpleList())
})

setAs("ACTIONetExperiment", "SingleCellExperiment", function(from) {	
	sce = as(as(from, "SummarizedExperiment"), "SingleCellExperiment")
    #sce = as(from, "SingleCellExperiment")
    reducedDims(sce)=SimpleList(lapply(colFactors(from), function(x) Matrix::t(x)))
    return(sce)
})

setAs("SingleCellExperiment", "ACTIONetExperiment", function(from) {	
	ace = as(as(from, "SummarizedExperiment"), "ACTIONetExperiment")
    #ace = as(from, "ACTIONetExperiment")
    
    colFactors(ace)=SimpleList(lapply(reducedDims(from), function(x) Matrix::t(x)))
    return(ace)
})

