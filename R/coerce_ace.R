#' Coerces a SummarizedExperiment method as an ACTIONetExperiment (ACE) object
#'
#' @param from SummarizedExperiment object
#'
#' @exportMethod coerce
setAs("SummarizedExperiment", "ACTIONetExperiment", function(from) {
    new("ACTIONetExperiment", from, 
    rowNets=S4Vectors::SimpleList(), 
    colNets=S4Vectors::SimpleList(), 
    rowMaps=S4Vectors::SimpleList(), 
    colMaps=S4Vectors::SimpleList())
})

setAs("ACTIONetExperiment", "SingleCellExperiment", function(from) {	
	sce = as(as(from, "SummarizedExperiment"), "SingleCellExperiment")
    transposed_factors = SimpleList(lapply(colMaps(from), function(x) Matrix::t(x)))
    reducedDims(sce) = transposed_factors
    return(sce)
})

setAs("SingleCellExperiment", "ACTIONetExperiment", function(from) {	
	ace = as(as(from, "SummarizedExperiment"), "ACTIONetExperiment")
    #ace = as(from, "ACTIONetExperiment")
    transposed_factors = SimpleList(lapply(reducedDims(from), function(x) Matrix::t(x)))
    colMaps(ace) = transposed_factors
    return(ace)
})

reconstruct_ace <- function(from) {	
	ace = as(as(from, "SummarizedExperiment"), "ACTIONetExperiment")
	
    rowNets(ace)=rowNets(from)
    colNets(ace)=colNets(from)
    rowMaps(ace)=rowMaps(from)
    colMaps(ace)=colMaps(from)

	return(ace)
}
