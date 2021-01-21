#' Takes a `ACTIONetExperiment` object and adds the reduced kernel matrix
#'
#' @param ace ACTIONetExperiment object.
#' @param reduced_dim Dimension of SVD used for reducing kernel matrix
#' @param max_iter Number of SVD iterations
#' @param assay_name Name of assay to reduce.
#' @param reduction_slot Name of slot to store reduction.
#'
#' @return ACTIONetExperiment object with reduction in colMaps(ace).
#'
#' @examples
#' ace = import.ace.from.10X(input_path)
#' ace = reduce.ace(ace)
#' @export
reduce.ace <- function(ace, reduced_dim = 50, max_iter = 10, assay_name = "logcounts",
    reduction_slot = "ACTION", seed = 0, SVD_algorithm = 0) {

    ace <- as(ace, "ACTIONetExperiment")
    if (!(assay_name %in% names(assays(ace)))) {
        err = sprintf("Assay '%s' not found.\n", assay_name)
        stop(err)
    }

    if (is.null(rownames(ace))) {
        rownames(ace) = .default_rownames(NROW(ace))
    } else {
        rownames(ace) = make.unique(rownames(ace), sep = "_")
    }

    if (is.null(colnames(ace))) {
        colnames(ace) = .default_colnames(NCOL(ace))
    } else {
        colnames(ace) = make.unique(colnames(ace), sep = "_")
    }

    # for (n in names(assays(ace))) { rownames(assays(ace)[[n]]) = rownames(ace)
    # colnames(assays(ace)[[n]]) = colnames(ace) }

    # msg = sprintf('Running main reduction.\n') message(msg)
    # reduction_algorithm=ACTION (1), SVD_algorithm=IRLB (0)
    if (SVD_algorithm == 0)
        max_iter = 100 * max_iter

    S = assays(ace)[[assay_name]]
    if (is.matrix(S)) {
        reduction.out = reduce_kernel_full(S, reduced_dim = reduced_dim, iter = max_iter,
            seed = seed, SVD_algorithm = SVD_algorithm)
    } else {
        reduction.out = reduce_kernel(S, reduced_dim = reduced_dim, iter = max_iter,
            seed = seed, SVD_algorithm = SVD_algorithm)
    }

    S_r = reduction.out$S_r
    colnames(S_r) = colnames(ace)
    # rownames(S_r) = sapply(1:nrow(S_r), function(i) sprintf('dim_%d', i))
    rownames(S_r) = paste0("dim_", 1:NROW(S_r))
    colMaps(ace)[[reduction_slot]] <- Matrix::t(S_r)
    colMapTypes(ace)[[reduction_slot]] = "reduction"

    V = reduction.out$V
    # colnames(V) = sapply(1:dim(V)[2], function(i) sprintf('V%d', i))
    colnames(V) = paste0("V", 1:NCOL(V))
    rowMaps(ace)[[sprintf("%s_V", reduction_slot)]] = V
    rowMapTypes(ace)[[sprintf("%s_V", reduction_slot)]] = "internal"

    A = reduction.out$A
    # colnames(A) = sapply(1:dim(A)[2], function(i) sprintf('A%d', i))
    colnames(A) = paste0("A", 1:NCOL(A))
    rowMaps(ace)[[sprintf("%s_A", reduction_slot)]] = A
    rowMapTypes(ace)[[sprintf("%s_A", reduction_slot)]] = "internal"

    B = reduction.out$B
    # colnames(B) = sapply(1:dim(B)[2], function(i) sprintf('B%d', i))
    colnames(B) = paste0("B", 1:NCOL(B))
    colMaps(ace)[[sprintf("%s_B", reduction_slot)]] = B
    colMapTypes(ace)[[sprintf("%s_B", reduction_slot)]] = "internal"

    metadata(ace)[[sprintf("%s_sigma", reduction_slot)]] = reduction.out$sigma

    return(ace)
}
