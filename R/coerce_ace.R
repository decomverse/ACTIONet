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
    sce = new("SingleCellExperiment", as(from, "RangedSummarizedExperiment")) 
    reducedDims(sce)=colFactors(from)
    return(sce)
})

setAs("SingleCellExperiment", "ACTIONetExperiment", function(from) {	
    ace = new("ACTIONetExperiment", as(from, "RangedSummarizedExperiment"))
    colFactors(ace)=reducedDims(from)
    return(ace)
})

