#' Perform batch correction on `ACTIONetExperiment` and `SingleCellExperiment` objects.
#' @export
reduce.and.batch.correct.ace.fastMNN <- function(ace, batch_attr, reduced_dim = 50, MNN_k = 20, return_V = FALSE, reduction_slot = "MNN",
    BPPARAM = SerialParam()) {
    .check_and_load_package("batchelor")

    ace = .check_if_ace(ace)
    SummarizedExperiment::assays(ace)[["counts"]] = as(SummarizedExperiment::assays(ace)[["counts"]], "sparseMatrix")
    m_data = metadata(ace)
    IDX = .get_ace_split_IDX(ace, batch_attr)

    ace.list = lapply(IDX, function(idx) computeSumFactors(ace[, idx], BPPARAM = BPPARAM))
    ace.list.norm = do.call(batchelor::multiBatchNorm, list(ace.list, BPPARAM = BPPARAM))
    # Sort based on 'complexity'
    merge_order = order(sapply(ace.list.norm, function(ace) dim(ace)[2]), decreasing = TRUE)

    ace.norm = do.call(cbind, ace.list.norm)
    SummarizedExperiment::assays(ace.norm)[["logcounts"]] = as(SummarizedExperiment::assays(ace.norm)[["logcounts"]], "sparseMatrix")

    set.seed(0)
    mnn.out <- do.call(batchelor::fastMNN, c(ace.list.norm, list(k = MNN_k, d = reduced_dim, auto.merge = FALSE, merge.order = merge_order, cos.norm = FALSE,
        assay.type = "logcounts", BPPARAM = BPPARAM)))

    S_r = reducedDims(mnn.out)[["corrected"]]
    rownames(S_r) = colnames(ace.norm)
    colnames(S_r) = sapply(1:ncol(S_r), function(i) sprintf("PC%d", i))


    ACTIONet::colMaps(ace.norm)[[reduction_slot]] <- Matrix::t(S_r)
    colMapTypes(ace)[[reduction_slot]] = "reduction"


    if (return_V) {
        V = rowData(mnn.out)[["rotation"]]
        colnames(V) = sapply(1:dim(V)[2], function(i) sprintf("PC%d", i))
        rowMaps(ace.norm)[["rotation"]] = V
    }
    metadata(ace.norm) = m_data
    return(ace.norm)
}

#' (It used Harmony for batch-correction: https://github.com/immunogenomics/harmony)
#'
#' @param ace Input ace object
#' @param norm.method Normalization method to use. See normalize.ace() function (default:'default')
#' (used only if the ace object is not already normalized)
#' @param batch_attr Vector of batches per sample
#' @param reduced_dim Dimension of SVD used for reducing kernel matrix
#' @param max_iter Number of SVD iterations
#' @param passphrase Passphrase for encrypting column names of the ace object for anonymization
#'
#' @return Reduced ace object with added colMaps(ace)
#'
#' @examples
#' ace = import.ace.from.10X(input_path)
#' batch_attr = ace$Batch # Assumes sample annotations are in the input_path with 'Batch' attribute being provided
#' ace = reduce.and.batch.correct.ace.Harmony(ace)
reduce.and.batch.correct.ace.Harmony <- function(ace, batch_attr, reduced_dim = 50, max_iter = 5, data.slot = "logcounts", norm_method = c("default", "scran", "Linnorm"),
    reduction_slot = "ACTION", seed = 0, SVD_algorithm = 1) {
    if(!require(harmony)){
      err = sprintf("You need to install harmony (https://github.com/immunogenomics/harmony) first for batch-correction.\n")
      stop(err)
    }

    ace <- .check_if_ace(ace)

    ace = reduce.ace(ace, reduced_dim = reduced_dim, max_iter = max_iter, norm_method = norm_method, data.slot = data.slot,
        reduction_slot = reduction_slot, seed = seed, SVD_algorithm = SVD_algorithm)
    ace = batch.correct.ace.Harmony(ace, batch_attr, reduction_slot = reduction_slot)

    return(ace)
}

#' @export
batch.correct.ace.Harmony <- function(ace, batch_attr, reduction_slot = "ACTION") {
    if(!require(harmony)){
      err = sprintf("You need to install harmony (https://github.com/immunogenomics/harmony) first for batch-correction.\n")
      stop(err)
    }
    ace <- .check_if_ace(ace)
    batch_attr = .get_ace_split_IDX(ace, batch_attr, return_split_vec = TRUE)
    ACTIONet::colMaps(ace)[[reduction_slot]] = harmony::HarmonyMatrix(ACTIONet::colMaps(ace)[[reduction_slot]], meta_data = batch_attr, do_pca = FALSE)
    return(ace)
}
