#' Perform batch correction on `ACTIONetExperiment` and `SingleCellExperiment` objects.
#' @export
reduce.and.batch.correct.ace.fastMNN <- function(ace, batch_attr, reduced_dim = 50,
    MNN_k = 20, return_V = FALSE, reduction_slot = "MNN", V_slot = NULL, BPPARAM = SerialParam()) {
    .check_and_load_package(c("scran", "SingleCellExperiment", "batchelor", "BiocParallel"))

    ace = .check_and_convert_se_like(ace, "ACE")
    SummarizedExperiment::assays(ace)[["counts"]] = as(SummarizedExperiment::assays(ace)[["counts"]],
        "sparseMatrix")
    m_data = metadata(ace)
    IDX = .get_ace_split_IDX(ace, batch_attr)

    sce.list = lapply(IDX, function(idx) {
        sce = as(ace[, idx], "SingleCellExperiment")
        sce = computeSumFactors(sce, BPPARAM = BPPARAM)
    })
    sce.list.norm = do.call(batchelor::multiBatchNorm, list(sce.list, BPPARAM = BPPARAM))
    # Sort based on 'complexity'
    merge_order = order(sapply(sce.list.norm, function(sce) dim(sce)[2]), decreasing = TRUE)

    sce.norm = do.call(cbind, sce.list.norm)
    if (!all(colnames(ace) == colnames(sce.norm)))
        sce.norm = sce.norm[, match(colnames(ace), colnames(sce.norm))]

    assays(ace)[["logcounts"]] <- assays(sce.norm)[["logcounts"]]
    SummarizedExperiment::assays(ace)[["logcounts"]] = as(SummarizedExperiment::assays(ace)[["logcounts"]],
        "sparseMatrix")
    sizeFactors(ace) = sizeFactors(sce.norm)

    rm("sce.norm")
    invisible(gc())

    set.seed(0)
    mnn.out <- do.call(batchelor::fastMNN, c(sce.list.norm, list(k = MNN_k, d = reduced_dim,
        auto.merge = FALSE, merge.order = merge_order, cos.norm = FALSE, assay.type = "logcounts",
        BPPARAM = BPPARAM)))

    S_r = Matrix::t(SingleCellExperiment::reducedDims(mnn.out)[["corrected"]])
    if (!all(colnames(ace) == colnames(S_r)))
        S_r = S_r[, match(colnames(ace), colnames(S_r))]

    # rownames(S_r) = colnames(ace)
    rownames(S_r) = sapply(1:nrow(S_r), function(i) sprintf("PC%d", i))

    ACTIONet::colMaps(ace)[[reduction_slot]] <- Matrix::t(S_r)
    colMapTypes(ace)[[reduction_slot]] = "reduction"


    if (return_V) {
        V = rowData(mnn.out)[["rotation"]]
        colnames(V) = sapply(1:dim(V)[2], function(i) sprintf("PC%d", i))
        if (is.null(V_slot) | length(V_slot) > 1) {
            V_slot = paste(reduction_slot, "rotation", sep = "_")
        }
        rowMaps(ace)[[V_slot]] = V
    }
    rm("sce.list.norm", "mnn.out")
    invisible(gc())

    metadata(ace) = m_data
    return(ace)
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
reduce.and.batch.correct.ace.Harmony <- function(ace, batch_attr, reduced_dim = 50,
    max_iter = 5, data_slot = "logcounts", norm_method = c("default", "scran", "Linnorm"),
    reduction_slot = "ACTION", seed = 0, SVD_algorithm = 0) {
    if (!require(harmony)) {
        err = sprintf("You need to install harmony (https://github.com/immunogenomics/harmony) first for batch-correction.\n")
        stop(err)
    }

    ace = .check_and_convert_se_like(ace, "ACE")

    ace = reduce.ace(ace, reduced_dim = reduced_dim, max_iter = max_iter, norm_method = norm_method,
        data_slot = data_slot, reduction_slot = reduction_slot, seed = seed, SVD_algorithm = SVD_algorithm)
    ace = batch.correct.ace.Harmony(ace, batch_attr, reduction_slot = sprintf("%s_harmony", reduction_slot))

    return(ace)
}

#' @export
batch.correct.ace.Harmony <- function(ace, batch_attr, reduction_slot = "ACTION") {
    if (!require(harmony)) {
        err = sprintf("You need to install harmony (https://github.com/immunogenomics/harmony) first for batch-correction.\n")
        stop(err)
    }
    ace = .check_and_convert_se_like(ace, "ACE")
    batch_attr = .get_ace_split_IDX(ace, batch_attr, return_split_vec = TRUE)
    ACTIONet::colMaps(ace)[[reduction_slot]] = harmony::HarmonyMatrix(ACTIONet::colMaps(ace)[[reduction_slot]],
        meta_data = batch_attr, do_pca = FALSE)
    return(ace)
}
