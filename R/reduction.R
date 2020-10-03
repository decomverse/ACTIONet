#' Takes a `ACTIONetExperiment` object and adds the reduced kernel matrix
#'
#' @param ace Input ace object
#' @param norm.method Normalization method to use. See normalize.ace() function (default:'default')
#' (used only if the ace object is not already normalized)
#' @param reduced_dim Dimension of SVD used for reducing kernel matrix
#' @param max_iter Number of SVD iterations
#' @param passphrase Passphrase for encrypting column names of the ace object for anonymization
#'
#' @return ACTIONetExperiment object (ACE), derived from SummarizedExperiment (ace), with added colMaps(ace)[['S_r']]
#'
#' @examples
#' ace = import.ace.from.10X(input_path)
#' ace = reduce.ace(ace)
#' @export
reduce.ace <- function(ace, reduced_dim = 50, max_iter = 5, data_slot = "logcounts",
    norm_method = "default", reduction_slot = "ACTION", seed = 0, SVD_algorithm = 0) {
      
    ace <- as(ace, "ACTIONetExperiment")
    if (!(data_slot %in% names(assays(ace)))) {
      err = sprintf("Slot %s not found.\n", data_slot)
      stop(err)
    }

    if (norm_method != "none") {
        msg = sprintf("Normalizing ace object ...\n")
        message(msg)
        ace = normalize.ace(ace, norm_method)
    }

    ace.norm = ace
    if (is.null(rownames(ace.norm))) {
        rownames(ace.norm) = sapply(1:nrow(ace.norm), function(i) sprintf("Gene%d",
            i))
    } else {
        rn = rownames(ace.norm)
        if (length(unique(rn)) < length(rn)) {
            rownames(ace.norm) = make.names(rn, unique = TRUE)
        }
    }

    if (is.null(colnames(ace.norm))) {
        colnames(ace.norm) = sapply(1:ncol(ace.norm), function(i) sprintf("Cell%d",
            i))
    } else {
        cn = colnames(ace.norm)
        if (length(unique(cn)) < length(cn)) {
            colnames(ace.norm) = make.names(cn, unique = TRUE)
        }
    }

    for (n in names(assays(ace.norm))) {
        rownames(assays(ace.norm)[[n]]) = rownames(ace.norm)
        colnames(assays(ace.norm)[[n]]) = colnames(ace.norm)
    }

    msg = sprintf("Running main reduction.\n")
    message(msg)
    # reduction_algorithm=ACTION (1), SVD_algorithm=IRLB (0)
    if (SVD_algorithm == 0)
        max_iter = 100*max_iter

    S = assays(ace.norm)[[data_slot]]
    if (is.matrix(S)) {
        reduction.out = reduce_kernel_full(S, reduced_dim = reduced_dim, iter = max_iter,
            seed = seed, SVD_algorithm = SVD_algorithm)
    } else {
        reduction.out = reduce_kernel(S, reduced_dim = reduced_dim, iter = max_iter,
            seed = seed, SVD_algorithm = SVD_algorithm)
    }

    S_r = reduction.out$S_r
    colnames(S_r) = colnames(ace.norm)
    rownames(S_r) = sapply(1:nrow(S_r), function(i) sprintf("Dim%d", i))
    colMaps(ace.norm)[[reduction_slot]] <- Matrix::t(S_r)
    colMapTypes(ace.norm)[[reduction_slot]] = "reduction"


	V = reduction.out$V
	colnames(V) = sapply(1:dim(V)[2], function(i) sprintf("V%d", i))
	rowMaps(ace.norm)[[sprintf("%s_V", reduction_slot)]] = V
	rowMapTypes(ace.norm)[[sprintf("%s_V", reduction_slot)]] = "internal"


	A = reduction.out$A
	colnames(A) = sapply(1:dim(A)[2], function(i) sprintf("A%d", i))
	rowMaps(ace.norm)[[sprintf("%s_A", reduction_slot)]] = A
	rowMapTypes(ace.norm)[[sprintf("%s_A", reduction_slot)]] = "internal"


	B = reduction.out$B
	colnames(B) = sapply(1:dim(B)[2], function(i) sprintf("B%d", i))
	colMaps(ace.norm)[[sprintf("%s_B", reduction_slot)]] = B
	colMapTypes(ace.norm)[[sprintf("%s_B", reduction_slot)]] = "internal"


	metadata(ace.norm)[[sprintf("%s_sigma", reduction_slot)]] = reduction.out$sigma

    # metadata(ace.norm)$reduction.time = Sys.time()
    return(ace.norm)
}
