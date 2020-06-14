#' Perform batch correction on `ACTIONetExperiment` and `SingleCellExperiment` objects.

reduce.and.batch.correct.sce.fastMNN <- function(sce, batch.attr = NULL, reduced_dim = 50, MNN.k = 20, return_V = FALSE, reduction_name = "MNN", BPPARAM = SerialParam()) {
  require(scran)
  require(batchelor)
  colData(sce) = droplevels(colData(sce))
  SummarizedExperiment::assays(sce)[["counts"]] = as(SummarizedExperiment::assays(sce)[["counts"]], 'sparseMatrix')

  # if( length(batch.attr) ==  1) {
  #   IDX = split(1:dim(sce)[2], colData(sce)[[batch.attr]])
  # } else {
  #   IDX = split(1:dim(sce)[2], batch.attr)
  # }
  IDX = .get_ace_split_IDX(ace, attr)

  sce.list = lapply(IDX, function(idx) computeSumFactors(sce[, idx], BPPARAM = BPPARAM ) )
  sce.list.norm = do.call(batchelor::multiBatchNorm, list(sce.list, BPPARAM = BPPARAM))

  # Sort based on 'complexity'
  merge_order = order(sapply(sce.list.norm, function(sce) dim(sce)[2]), decreasing = TRUE)

  sce.norm = do.call(cbind, sce.list.norm)
  SummarizedExperiment::assays(sce.norm)[["logcounts"]] = as(SummarizedExperiment::assays(sce.norm)[["logcounts"]], "sparseMatrix")

  set.seed(0)
  mnn.out <- do.call(batchelor::fastMNN, c(sce.list.norm, list(k = 20, d = 50, auto.merge = FALSE, merge.order = merge_order, cos.norm = FALSE, assay.type = "logcounts", BPPARAM = BPPARAM)))

  S_r = reducedDims(mnn.out)[["corrected"]]
  rownames(S_r) = colnames(sce.norm)
  colnames(S_r) = sapply(1:ncol(S_r), function(i) sprintf("PC%d", i))


  SingleCellExperiment::reducedDims(sce.norm)[[reduction_name]] <- S_r

  if(return_V){
    V = rowData(mnn.out)[["rotation"]]
    colnames(V) = sapply(1:dim(V)[2], function(i) sprintf("PC%d", i))
    rowData(sce.norm)[["rotation"]] = V
  }

  return(sce.norm)
}

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
#' @return Reduced sce object with added ReducedDims(sce)
#'
#' @examples
#' sce = import.sce.from.10X(input_path)
#' batch.vec = sce$Batch # Assumes sample annotations are in the input_path with "Batch" attribute being provided
#' sce = reduce.and.batch.correct.sce.Harmony(sce)
reduce.and.batch.correct.sce.Harmony <- function(sce, batch.vec = NULL, reduced_dim = 50, max.iter = 5, data.slot = "logcounts", normalization.method = "default", reduction.slot = "ACTION", seed = 0, SVD_algorithm = 0) {
	if( !("harmony" %in% rownames(installed.packages())) ) {
		message("You need to install harmony (https://github.com/immunogenomics/harmony) first for batch-correction.")
		return
	} else {
		require(harmony)
	}

    if (is.null(batch.vec)) {
        err = sprintf("You need to provide the batch vector/attr.")
        stop(err)
    }

    sce = reduce.sce(sce, reduced_dim = reduced_dim, max.iter = max.iter, normalization.method = normalization.method, data.slot = data.slot, reduction.slot = reduction.slot, seed = seed, SVD_algorithm = SVD_algorithm)
    sce = batch.correct.sce.Harmony(sce, batch.vec)

    return(sce)
}

batch.correct.sce.Harmony <- function(sce, batch.vec, reduction.slot = "ACTION") {
    require(harmony)
    reducedDims(sce)[[reduction.slot]] = harmony::HarmonyMatrix(reducedDims(sce)[[reduction.slot]], batch.vec, do_pca = FALSE)
    return(sce)
}
