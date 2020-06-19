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

setAs("ACTIONetExperiment", "SingleCellExperiment", function(from) {	
	SE = as(from, "SummarizedExperiment")	
	sce = as(SE, "SingleCellExperiment")
	
    transposed_factors = SimpleList(lapply(colMaps(from), function(x) Matrix::t(x)))
    reducedDims(sce) = transposed_factors
    
    #rowData(sce) = DataFrame(as.data.frame(rowData(from)))
    #colData(sce) = DataFrame(as.data.frame(colData(from)))
        
    return(sce)
})

setAs("SingleCellExperiment", "ACTIONetExperiment", function(from) {	
	SE = as(from, "SummarizedExperiment")
	rownames(SE) = rownames(from)
	rowData(SE) = rowData(from)
	
	
	ace = as(SE, "ACTIONetExperiment")
    #ace = as(from, "ACTIONetExperiment")
    
    transposed_factors = SimpleList(lapply(reducedDims(from), function(x) Matrix::t(x)))
    colMaps(ace) = transposed_factors
    
    rowData(ace) = DataFrame(as.data.frame(rowData(ace)))
    colData(ace) = DataFrame(as.data.frame(colData(ace)))
    
    return(ace)
})

reconstruct_ace <- function(from) {	
	SE = as(from, "SummarizedExperiment")	
	ace = as(SE, "ACTIONetExperiment")
	
    rowData(ace) = DataFrame(as.data.frame(rowData(ace)))
    colData(ace) = DataFrame(as.data.frame(colData(ace)))
    
    
    rowNets(ace)=rowNets(from)
    colNets(ace)=colNets(from)
    if('rowMaps' %in% slotNames(SE)) {
		rowMaps(ace)=rowMaps(from)
		colMaps(ace)=colMaps(from)
	} else {
		rowMaps(ace)=from@rowFactors
		colMaps(ace)=from@colFactors
	}

    
	return(ace)
}
