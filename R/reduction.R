#' Takes a `ACTIONetExperiment` object and adds the reduced kernel matrix
#'
#' @param ace Input ace object
#' @param norm.method Normalization method to use. See normalize.ace() function (default:"default")
#' (used only if the ace object is not already normalized)
#' @param reduced_dim Dimension of SVD used for reducing kernel matrix
#' @param max.iter Number of SVD iterations
#' @param passphrase Passphrase for encrypting column names of the ace object for anonymization
#'
#' @return ACTIONetExperiment object (ACE), derived from SummarizedExperiment (ace), with added colMaps(ace)[["S_r"]]
#'
#' @examples
#' ace = import.ace.from.10X(input_path)
#' ace = reduce.ace(ace)
reduce.ace <- function(ace, reduced_dim = 50, max.iter = 5, data.slot = "logcounts", normalization.method = "default", reduction.slot = "ACTION", seed = 0, SVD_algorithm = 1, return_V = FALSE) {
    ace <- check_if_ace(ace)
    if (!(data.slot %in% names(assays(ace)))) {
		if(normalization.method != "none") {
			msg = sprintf("Normalizing ace object ...\n")
      message(msg)
			ace = normalize.ace(ace, normalization.method)
		} else {
			err = sprintf("Slot %s not found. This can be potentially due to missing normalization step.", data.slot)
      stop(err)
		}
    }

    ace.norm = ace
    if (is.null(rownames(ace.norm))) {
        rownames(ace.norm) = sapply(1:nrow(ace.norm), function(i) sprintf("Gene%d", i))
    } else {
		rn = rownames(ace.norm)
		if(length(unique(rn)) < length(rn)) {
			rownames(ace.norm) = make.names(rn, unique = TRUE)
		}
	}

    if (is.null(colnames(ace.norm))) {
        colnames(ace.norm) = sapply(1:ncol(ace.norm), function(i) sprintf("Cell%d", i))
    } else {
		cn = colnames(ace.norm)
		if(length(unique(cn)) < length(cn)) {
			colnames(ace.norm) = make.names(cn, unique = TRUE)
		}
	}

	for(n in names(assays(ace.norm))) {
		rownames(assays(ace.norm)[[n]]) = rownames(ace.norm)
		colnames(assays(ace.norm)[[n]]) = colnames(ace.norm)
	}

    msg = sprintf("Running main reduction.\n")
    message(msg)
    # reduction_algorithm=ACTION (1), SVD_algorithm=IRLB (0)
    if(SVD_algorithm == 0)
		max.iter = max.iter * 100

    S = assays(ace.norm)[[data.slot]]
    if(is.matrix(S)) {
        reduction.out = reduce_kernel_full(S, reduced_dim = reduced_dim, iter = max.iter, seed = seed, reduction_algorithm = 1, SVD_algorithm = SVD_algorithm)
	} else {
        reduction.out = reduce_kernel(S, reduced_dim = reduced_dim, iter = max.iter, seed = seed, reduction_algorithm = 1, SVD_algorithm = SVD_algorithm)
    }

    S_r = reduction.out$S_r
    colnames(S_r) = colnames(ace.norm)
    rownames(S_r) = sapply(1:nrow(S_r), function(i) sprintf("Dim%d", i))
    colMaps(ace.norm)[[reduction.slot]] <- S_r
	ace.norm@colMapsAnnot[[reduction.slot]] = list(type = "reduction")								
    

	if(return_V){
		V = reduction.out[["V"]]
		colnames(V) = sapply(1:dim(V)[2], function(i) sprintf("PC%d", i))
		rowMaps(ace.norm)[["rotation"]] = V
		ace.norm@rowMapsAnnot[["rotation"]] = list(type = "internal")								
	}

    metadata(ace.norm)$reduction.time = Sys.time()
    return(ace.norm)
}
