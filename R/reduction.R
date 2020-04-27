#' Takes a `SingleCellExperiment` object and adds the reduced kernel matrix
#'
#' @param sce Input sce object
#' @param norm.method Normalization method to use. See normalize.sce() function (default:"default")
#' (used only if the sce object is not already normalized)
#' @param reduced_dim Dimension of SVD used for reducing kernel matrix
#' @param max.iter Number of SVD iterations
#' @param passphrase Passphrase for encrypting column names of the sce object for anonymization
#' 
#' @return ACTIONetExperiment object (ACE), derived from SingleCellExperiment (SCE), with added ReducedDims(sce)[["S_r"]]
#' 
#' @examples
#' sce = import.sce.from.10X(input_path)
#' sce = reduce.sce(sce)
reduce.sce <- function(sce, norm.method = "default", reduced_dim = 50, max.iter = 5, passphrase = NULL) {        
    sce.norm = sce
    
    if (is.null(rownames(sce.norm))) {
        rownames(sce.norm) = sapply(1:nrow(sce.norm), function(i) sprintf("Gene%d", i))
    }
    if (is.null(colnames(sce.norm))) {
        colnames(sce.norm) = sapply(1:ncol(sce.norm), function(i) sprintf("Cell%d", i))
    }
    
    
    if (!("logcounts" %in% names(SummarizedExperiment::assays(sce.norm)))) {
        print("Normalizing sce object")
        
        sce.norm = normalize.sce(sce.norm, norm.method)
    }
    SummarizedExperiment::assays(sce.norm) = lapply(SummarizedExperiment::assays(sce.norm), function(A) {
		rownames(A) = rownames(sce.norm)
		colnames(A) = colnames(sce.norm)
		return(A)
	})
        
    
    print("Running main reduction")
    # reduction_algorithm=ACTION, SVD_algorithm=Halko
    suppressWarnings({
        reduction.out = reduce_kernel(as(SummarizedExperiment::assays(sce.norm)$logcounts, "sparseMatrix"), reduced_dim = reduced_dim, iter = max.iter, seed = 0, reduction_algorithm = 1, SVD_algorithm = 1) 
    })
    
    S_r = t(reduction.out$S_r)
    rownames(S_r) = colnames(sce.norm)
    colnames(S_r) = sapply(1:ncol(S_r), function(i) sprintf("Dim%d", i))    
    reducedDim(sce.norm, "ACTION") <- S_r
    
    
    metadata(sce.norm)$reduction.time = Sys.time()
    return(sce.norm)
}

batch.correct.sce.Harmony <- function(sce, batch.vec) {
    require(harmony)
    reducedDims(sce)$S_r = harmony::HarmonyMatrix(reducedDims(sce)$S_r, batch.vec, do_pca = FALSE)
    return(sce)
}


#' Takes a `SingleCellExperiment` object, batch corrects, and adds the reduced kernel matrix
#' (It used Harmony for batch-correction: https://github.com/immunogenomics/harmony)
#'
#' @param sce Input sce object
#' @param norm.method Normalization method to use. See normalize.sce() function (default:"default")
#' (used only if the sce object is not already normalized)
#' @param batch.vec Vector of batches per sample
#' @param reduced_dim Dimension of SVD used for reducing kernel matrix
#' @param max.iter Number of SVD iterations
#' @param passphrase Passphrase for encrypting column names of the sce object for anonymization
#' 
#' @return Reduced sce object with added ReducedDims(sce)[["S_r"]]
#' 
#' @examples
#' sce = import.sce.from.10X(input_path)
#' batch.vec = sce$Batch # Assumes sample annotations are in the input_path with "Batch" attribute being provided
#' sce = reduce.and.batch.correct.sce.Harmony(sce)
reduce.and.batch.correct.sce.Harmony <- function(sce, batch.vec = NULL, norm.method = "default", reduced_dim = 50, max.iter = 5, passphrase = NULL) {
	if( !("harmony" %in% rownames(installed.packages())) ) {
		message("You need to install harmony (https://github.com/immunogenomics/harmony) first for batch-correction.")
		return
	} else {
		library(harmony)
	}
	
    if (is.null(batch.vec)) {
        print("You need to provide the batch vector/attr")
        return(sce)
    }
    
    sce = reduce.sce(sce, reduced_dim = reduced_dim, max.iter = max.iter, norm.method = norm.method, passphrase = passphrase)
    sce = batch.correct.sce.Harmony(sce, batch.vec)
    
    return(sce)
}
