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
reduce.sce <- function(sce, reduced_dim = 50, max.iter = 5, data.slot = "logcounts", normalization.method = "default", reduction.slot = "ACTION", seed = 0, SVD_algorithm = 0) {
    if (!(data.slot %in% names(assays(sce)))) {
		if(normalization.method != "none") {
			print("Normalizing sce object ... ");
			sce = normalize.sce(sce, normalization.method)
		} else {
			R.utils::printf("Slot %s not found. This can be potentially due to missing normalization step.", data.slot)
			return(sce)
		}
    }

    sce.norm = sce
    if (is.null(rownames(sce.norm))) {
        rownames(sce.norm) = sapply(1:nrow(sce.norm), function(i) sprintf("Gene%d", i))
    } else {
		rn = rownames(sce.norm)
		if(length(unique(rn)) < length(rn)) {
			rownames(sce.norm) = make.names(rn, unique = TRUE)
		}
	}

    if (is.null(colnames(sce.norm))) {
        colnames(sce.norm) = sapply(1:ncol(sce.norm), function(i) sprintf("Cell%d", i))
    } else {
		cn = colnames(sce.norm)
		if(length(unique(cn)) < length(cn)) {
			colnames(sce.norm) = make.names(cn, unique = TRUE)
		}
	}




    for(n in names(assays(sce.norm))) {
		rownames(assays(sce.norm)[[n]]) = rownames(sce.norm)
		colnames(assays(sce.norm)[[n]]) = colnames(sce.norm)
	}


    print("Running main reduction")
    # reduction_algorithm=ACTION (1), SVD_algorithm=IRLB (0)
    if(SVD_algorithm == 0)
		max.iter = max.iter * 100
		
    S = assays(sce.norm)[[data.slot]]
    if(is.matrix(S)) {
        reduction.out = reduce_kernel_full(S, reduced_dim = reduced_dim, iter = max.iter, seed = seed, reduction_algorithm = 1, SVD_algorithm = SVD_algorithm)
	} else {
        reduction.out = reduce_kernel(S, reduced_dim = reduced_dim, iter = max.iter, seed = seed, reduction_algorithm = 1, SVD_algorithm = SVD_algorithm)
    }

    S_r = t(reduction.out$S_r)
    rownames(S_r) = colnames(sce.norm)
    colnames(S_r) = sapply(1:ncol(S_r), function(i) sprintf("Dim%d", i))
    reducedDims(sce.norm)[[reduction.slot]] <- S_r


    metadata(sce.norm)$reduction.time = Sys.time()
    return(sce.norm)
}

batch.correct.sce.Harmony <- function(sce, batch.vec, reduction.slot = "ACTION") {
    require(harmony)
    reducedDims(sce)[[reduction.slot]] = harmony::HarmonyMatrix(reducedDims(sce)[[reduction.slot]], batch.vec, do_pca = FALSE)
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
reduce.and.batch.correct.sce.Harmony <- function(sce, batch.vec = NULL, reduced_dim = 50, max.iter = 5, data.slot = "logcounts", normalization.method = "default", reduction_name = "ACTION", seed = 0, SVD_algorithm = 0) {
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

    sce = reduce.sce(sce, reduced_dim = reduced_dim, max.iter = max.iter, normalization.method = normalization.method, data.slot = data.slot, reduction.slot = reduction.slot, seed = seed, SVD_algorithm = SVD_algorithm)
    sce = batch.correct.sce.Harmony(sce, batch.vec)

    return(sce)
}
