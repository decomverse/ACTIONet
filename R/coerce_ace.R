#' Coerces a SummarizedExperiment method as an ACTIONetExperiment (ACE) object
#'
#' @param from SummarizedExperiment object
#'
#' @exportMethod coerce
setAs("SummarizedExperiment", "ACTIONetExperiment", function(from) {
    ace = .ACTIONetExperiment(from, rowNets=S4Vectors::SimpleList(), 
    colNets=S4Vectors::SimpleList(), 
    rowMaps=S4Vectors::SimpleList(), 
    colMaps=S4Vectors::SimpleList())

    #rowData(ace) = DataFrame(as.data.frame(rowData(ace)))
    #colData(ace) = DataFrame(as.data.frame(colData(ace)))

    return(ace)  
})

#' Coerces a ACTIONetExperiment (ACE) method as an SingleCellExperiment (SCE) object
#'
#' @param from ACTIONetExperiment object
#'
#' @exportMethod coerce
setAs("ACTIONetExperiment", "SingleCellExperiment", function(from) {	
	SE = as(from, "SummarizedExperiment")	
	sce = as(SE, "SingleCellExperiment")
	
	Xs = colMaps(from)
	Xs = Xs[colMapTypes(from) != "internal"]
	
    transposed_factors = as(lapply(Xs, function(X) Matrix::t(X)), "SimpleList")
    reducedDims(sce) = transposed_factors
    
    #rowData(sce) = DataFrame(as.data.frame(rowData(from)))
    #colData(sce) = DataFrame(as.data.frame(colData(from)))
        
    return(sce)
})

#' Coerces a SingleCellExperiment (SCE) method as an ACTIONetExperiment (ACE) object
#'
#' @param from SingleCellExperiment object
#'
#' @exportMethod coerce
setAs("SingleCellExperiment", "ACTIONetExperiment", function(from) {	
	SE = as(from, "SummarizedExperiment")
	rownames(SE) = rownames(from)
	rowData(SE) = rowData(from)
	
	
	ace = as(SE, "ACTIONetExperiment")
    #ace = as(from, "ACTIONetExperiment")
    
    transposed_factors = as(lapply(reducedDims(from), function(x) SummarizedExperiment(assays = list(X = Matrix::t(x)))), "SimpleList")
    colMaps(ace) = transposed_factors
    
    rowData(ace) = DataFrame(as.data.frame(rowData(ace)))
    colData(ace) = DataFrame(as.data.frame(colData(ace)))
    
    return(ace)
})

