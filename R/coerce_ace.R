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
	
    transposed_factors = SimpleList(lapply(colMaps(from), function(SE) Matrix::t(assays(SE)$X)))
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
    
    transposed_factors = SimpleList(lapply(reducedDims(from), function(x) SummarizedExperiment(assays = list(X = Matrix::t(x)))))
    colMaps(ace) = transposed_factors
    
    rowData(ace) = DataFrame(as.data.frame(rowData(ace)))
    colData(ace) = DataFrame(as.data.frame(colData(ace)))
    
    return(ace)
})

# reconstruct_ace <- function(from) {	
# 	SE = as(from, "SummarizedExperiment")	
# 	ace = as(SE, "ACTIONetExperiment")
# 	
#     rowData(ace) = DataFrame(as.data.frame(rowData(ace)))
#     colData(ace) = DataFrame(as.data.frame(colData(ace)))
#     
#     
#     rowNets(ace)=rowNets(from)
#     colNets(ace)=colNets(from)
#     
#     if(.hasSlot(from, "rowFactors")) {
# 		rowMaps(ace)=from@rowFactors
# 		colMaps(ace)=from@colFactors		
# 	} else {
# 		rowMaps(ace)=rowMaps(from)
# 		colMaps(ace)=colMaps(from)
# 	}
# 
# 	if( length(ace@rowMapsAnnot) == 0 ) {
# 		nn = names(rowMaps(ace))
# 		Annot = lapply(1:length(nn), function(i) list(type="reduction"))
# 		names(Annot) = nn
# 		
# 		for(i in grep("^C_|^H_", nn)) {
# 			Annot[[i]]$type = "internal"
# 		}
# 
# 		for(i in which(sapply(rowMaps(ace), ncol) <= 3)) {
# 			Annot[[i]]$type = "embedding"
# 		}
# 
# 		ace@rowMapsAnnot = SimpleList(Annot)
# 	}
# 
# 	if( length(ace@colMapsAnnot) == 0 ) {
# 		nn = names(colMaps(ace))
# 		Annot = lapply(1:length(nn), function(i) list(type="reduction"))
# 		names(Annot) = nn
# 		
# 		for(i in grep("^C_|^H_", nn)) {
# 			Annot[[i]]$type = "internal"
# 		}
# 
# 		for(i in which(sapply(colMaps(ace), ncol) <= 3)) {
# 			Annot[[i]]$type = "embedding"
# 		}
# 
# 		ace@colMapsAnnot = SimpleList(Annot)
# 	}
# 
#     
# 	return(ace)
# }

