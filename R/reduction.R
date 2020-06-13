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
reduce.sce <- function(sce, reduced_dim = 50, max.iter = 5, data.slot = "logcounts", normalization.method = "default", reduction.slot = "ACTION", seed = 0, SVD_algorithm = 1) {
    if (!(data.slot %in% names(assays(sce)))) {
		if(normalization.method != "none") {
			msg = sprintf("Normalizing sce object ...\n")
      message(msg)
			sce = normalize.sce(sce, normalization.method)
		} else {
			err = sprintf("Slot %s not found. This can be potentially due to missing normalization step.", data.slot)
      stop(err)
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

    msg = sprintf("Running main reduction.\n")
    message(msg)
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
