#' Takes a `SingleCellExperiment` object and adds the reduced kernel matrix
#'
#' @param sce Input sce object
#' @param norm.method Normalization method to use. See normalize.sce() function (default:"default")
#' (used only if the sce object is not already normalized)
#' @param reduced_dim Dimension of SVD used for reducing kernel matrix
#' @param max.iter Number of SVD iterations
#' @param passphrase Passphrase for encrypting column names of the sce object for anonymization
#' 
#' @return Reduced sce object with added ReducedDims(sce)[["S_r"]]
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
    
    
    if (!("cell.hashtag" %in% colnames(colData(sce.norm)))) {
        print("tagging cells")
        time.tag = Sys.time()
        if(is.null(passphrase))
			passphrase = as.character(time.tag)
        
        encoded.ids = encode.ids(colnames(sce.norm), passphrase)        
        sce.norm$cell.hashtag = cell.hashtags
        
        colData(sce.norm)$original.colnames = colnames(sce.norm)
        colnames(sce.norm) = sce.norm$cell.hashtag
        
        metadata(sce.norm)$tagging.time = time.tag
        metadata(sce.norm)$passphrase = passphrase
    }

    if (!("logcounts" %in% names(SummarizedExperiment::assays(sce.norm)))) {
        print("Normalizing sce object")
        
        sce.norm = normalize.sce(sce.norm, norm.method)
		rownames(SummarizedExperiment::assays(sce.norm)$counts) = rownames(SummarizedExperiment::assays(sce.norm)$logcounts) = rownames(sce.norm)
		colnames(SummarizedExperiment::assays(sce.norm)$counts) = colnames(SummarizedExperiment::assays(sce.norm)$logcounts) = colnames(sce.norm)
    } else {
		rownames(SummarizedExperiment::assays(sce.norm)$logcounts) = rownames(sce.norm)
		colnames(SummarizedExperiment::assays(sce.norm)$logcounts) = colnames(sce.norm)
    }
        
    
    print("Running main reduction")
    # reduction_algorithm=ACTION, SVD_algorithm=Halko
    suppressWarnings({
        reduction.out = reduce_kernel(as(SummarizedExperiment::assays(sce.norm)$logcounts, "sparseMatrix"), reduced_dim = reduced_dim, iter = max.iter, seed = 0, reduction_algorithm = 1, SVD_algorithm = 1) 
    })
    
    S_r = t(reduction.out$S_r)
    rownames(S_r) = colnames(sce.norm)
    colnames(S_r) = sapply(1:ncol(S_r), function(i) sprintf("Dim%d", i))
    
    SingleCellExperiment::reducedDim(sce.norm, "S_r") <- S_r
    
    V = reduction.out$V
    rownames(V) = rownames(sce.norm)
    colnames(V) = sapply(1:dim(V)[2], function(i) sprintf("Dim%d", i))
    
    X = rowData(sce.norm)
    PC.idx = -grep("^Dim", colnames(X))
    if (length(PC.idx) > 0) 
        X = X[, PC.idx]
    rowData(sce.norm) = cbind(V, X)
    
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


#' Takes a `SingleCellExperiment` object, batch corrects, and adds the reduced kernel matrix
#' (It used Mutua-Nearest-Neighbors (MNN) for batch-correction)
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
reduce.and.batch.correct.sce.MNN <- function(sce, batch.vec = NULL, norm.method = "scran", reduced_dim = 50, MNN.k = 20) {
    require(scran)
    
    if (is.null(batch.vec)) {
        print("You need to provide the batch vector/attr")
        return(sce)
    }
    
    if (is.null(rownames(sce))) {
        rownames(sce) = sapply(1:nrow(sce), function(i) sprintf("Gene%d", i))
    }
    if (is.null(colnames(sce))) {
        colnames(sce) = sapply(1:ncol(sce), function(i) sprintf("Cell%d", i))
    }
    
    if (!("logcounts" %in% names(SummarizedExperiment::assays(sce.norm)))) {
        print("Normalizing sce object")
                
		sce = clearSpikes(sce)
		assays(sce)[["counts"]] = as(assays(sce)[["counts"]], "sparseMatrix")		
		
		IDX = split(1:ncol(sce), batch.vec)
		
		sce.list = lapply(IDX, function(idx) suppressWarnings(normalize.sce(sce[, idx], norm.method)))
		sce.list.norm = do.call(scran::multiBatchNorm, sce.list)		
		
		# Sort based on 'complexity'
		perm = order(sapply(sce.list.norm, function(sce) dim(sce)[2]), decreasing = TRUE)
		sce.list.norm = sce.list.norm[perm]
		
		
		sce.norm = do.call(SingleCellExperiment::cbind, sce.list.norm)
		SummarizedExperiment::assays(sce.norm)$logcounts = as(SummarizedExperiment::assays(sce.norm)[["logcounts"]], "sparseMatrix")
		
		metadata(sce.norm)$normalization.method = "multiBatchNorm"
		metadata(sce.norm)$normalization.time = Sys.time()
		
		rownames(SummarizedExperiment::assays(sce.norm)$counts) = rownames(SummarizedExperiment::assays(sce.norm)$logcounts) = rownames(sce.norm)
		colnames(SummarizedExperiment::assays(sce.norm)$counts) = colnames(SummarizedExperiment::assays(sce.norm)$logcounts) = colnames(sce.norm)		
    } else {
		rownames(SummarizedExperiment::assays(sce.norm)$logcounts) = rownames(sce.norm)
		colnames(SummarizedExperiment::assays(sce.norm)$logcounts) = colnames(sce.norm)
    }
        	
    if (!("cell.hashtag" %in% colnames(colData(sce.norm)))) {
        print("tagging cells")
        time.tag = Sys.time()
        if(is.null(passphrase))
			passphrase = as.character(time.tag)
        
        encoded.ids = encode.ids(colnames(sce.norm), passphrase)        
        sce.norm$cell.hashtag = cell.hashtags
        
        colData(sce.norm)$original.colnames = colnames(sce.norm)
        colnames(sce.norm) = sce.norm$cell.hashtag
        
        metadata(sce.norm)$tagging.time = time.tag
        metadata(sce.norm)$passphrase = passphrase
    }
   
    set.seed(0)
    mnn.out <- do.call(scran::fastMNN, c(sce.list.norm, list(k = MNN.k, d = reduced_dim, auto.order = FALSE, approximate = TRUE, cos.norm = FALSE)))
    
    S_r = mnn.out$corrected
    rownames(S_r) = colnames(sce.norm)
    colnames(S_r) = sapply(1:ncol(S_r), function(i) sprintf("Dim%d", i))
    
    
    SingleCellExperiment::reducedDims(sce.norm)[["S_r"]] <- S_r
    
    V = mnn.out$rotation
    rownames(V) = rownames(sce.norm)
    colnames(V) = sapply(1:dim(V)[2], function(i) sprintf("Dim%d", i))
    
    colnames(V) = sapply(1:dim(V)[2], function(i) sprintf("Dim%d", i))
    
    X = rowData(sce.norm)
    PC.idx = grep("^Dim", colnames(X))
    if (length(PC.idx) > 0) 
        X = X[, -PC.idx]
    
    rowData(sce.norm) = cbind(V, X)
    
    metadata(sce.norm)$reduction.time = Sys.time()
    return(sce.norm)
}


