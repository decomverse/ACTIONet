#' Perform batch correction on `ACTIONetExperiment` and `SingleCellExperiment` objects.
#' @export
reduce.and.batch.correct.ace.fastMNN <- function(ace, batch_attr = NULL, assay_name = "logcounts", reduced_dim = 50, MNN_k = 20, return_V = FALSE, reduction_slot = "MNN", V_slot = NULL, BPPARAM = SerialParam()) {
    .check_and_load_package(c("scran", "SingleCellExperiment", "batchelor", "BiocParallel"))

	if(is.null(batch_attr)) {
		warning("'batch_attr' must be provided")
		return(ace)
	}

    ace = .check_and_convert_se_like(ace, "ACE")
    m_data = metadata(ace)
    ace = normalize.ace(ace, norm_method = "multiBatchNorm", assay_name = assay_name, batch_attr = batch_attr, BPPARAM = BPPARAM)

    S = SummarizedExperiment::assays(ace)[[assay_name]]
    mnn_batch = .get_ace_split_IDX(ace, batch_attr, return_split_vec = TRUE)
    IDX = .get_ace_split_IDX(ace, batch_attr)
    merge_order = order(sapply(IDX, function(idx) length(idx)), decreasing = TRUE)

    set.seed(0)
    mnn.out <- batchelor::fastMNN(S, batch = mnn_batch, k = MNN_k, d = reduced_dim, auto.merge = FALSE, merge.order = merge_order, cos.norm = FALSE, BPPARAM = BPPARAM)

    S_r = SingleCellExperiment::reducedDims(mnn.out)[["corrected"]]
    rownames(S_r) = colnames(ace)
    colnames(S_r) = sapply(1:dim(S_r)[2], function(i) sprintf("PC%d", i))

    ACTIONet::colMaps(ace)[[reduction_slot]] <- S_r
    colMapTypes(ace)[[reduction_slot]] = "reduction"

    if (return_V) {
        V = rowData(mnn.out)[["rotation"]]
        colnames(V) = sapply(1:dim(V)[2], function(i) sprintf("PC%d", i))
        if (is.null(V_slot) | length(V_slot) > 1) {
            V_slot = paste(reduction_slot, "rotation", sep = "_")
        }
        rowMaps(ace)[[V_slot]] = V
    }
    invisible(gc())

    metadata(ace) = m_data
    return(ace)
}

#' (It uses Harmony for batch-correction: https://github.com/immunogenomics/harmony)
#'
#' @param ace ACTIONetExperiment object
#' @param norm_method Normalization method to use. See normalize.ace() function (default:'default')
#' (used only if the ace object is not already normalized)
#' @param batch_attr Vector of length ncol(ace) or column name of colData(ace) containing batch labels.
#' @param reduced_dim Dimension of SVD used for reducing kernel matrix
#' @param max_iter Number of SVD iterations
#'
#' @return Reduced ace object with added colMaps(ace)
#'
#' @examples
#' ace = import.ace.from.10X(input_path)
#' batch_attr = ace$Batch # Assumes sample annotations are in the input_path with 'Batch' attribute being provided
#' ace = reduce.and.batch.correct.ace.Harmony(ace)
#' @export
reduce.and.batch.correct.ace.Harmony <- function(ace, batch_attr, reduced_dim = 50, max_iter = 10, assay_name = "logcounts", norm_method = c("default", "scran", "multiBatchNorm", "Linnorm"), reduction_slot = "ACTION", seed = 0, SVD_algorithm = 0) {
  if (!require(harmony)) {
      err = sprintf("You need to install harmony (https://github.com/immunogenomics/harmony) first for batch-correction.\n")
      stop(err)
  }

  if(is.null(batch_attr)) {
		err = sprintf("'batch_attr' must be provided.\n")
		stop(err)
	}

  ace = .check_and_convert_se_like(ace, "ACE")
  norm_method = match.arg(norm_method)
  batch_attr = .get_ace_split_IDX(ace, batch_attr, return_split_vec = TRUE)

  ace = reduce.ace(ace, reduced_dim = reduced_dim, max_iter = max_iter, norm_method = norm_method,
      assay_name = assay_name, reduction_slot = reduction_slot, seed = seed, SVD_algorithm = SVD_algorithm)

  ace = batch.correct.ace.Harmony(ace, batch_attr, reduction_slot = reduction_slot)

  return(ace)
}

#' @export
batch.correct.ace.Harmony <- function(ace, batch_attr = NULL, reduction_slot = "ACTION") {
  if (!require(harmony)) {
    err = sprintf("You need to install harmony (https://github.com/immunogenomics/harmony) first for batch-correction.\n")
    stop(err)
  }

	if(is.null(batch_attr)) {
		err = sprintf("'batch_attr' must be provided.\n")
		stop(err)
	}

  ace = .check_and_convert_se_like(ace, "ACE")
  batch_attr = .get_ace_split_IDX(ace, batch_attr, return_split_vec = TRUE)
  ACTIONet::colMaps(ace)[[reduction_slot]] = harmony::HarmonyMatrix(ACTIONet::colMaps(ace)[[reduction_slot]],
      meta_data = batch_attr, do_pca = FALSE)
  return(ace)
}


#' @export
orthogonalize.ace.batch <- function(ace, design_mat, reduction_slot = "ACTION", assay_name = "logcounts") {
	S = assays(ace)[[assay_name]]
	S_r = colMaps(ace)[[sprintf("%s", reduction_slot)]]
	V = rowMaps(ace)[[sprintf("%s_V", reduction_slot)]]
	A = rowMaps(ace)[[sprintf("%s_A", reduction_slot)]]
	B = colMaps(ace)[[sprintf("%s_B", reduction_slot)]]
	sigma = metadata(ace)[[sprintf("%s_sigma", reduction_slot)]]

	reduction.out = orthogonalize_batch_effect(S = S, old_S_r = S_r, old_V = V, old_A = A, old_B = B, old_sigma = sigma, design = design_mat)

    S_r = reduction.out$S_r
    colnames(S_r) = colnames(ace)
    rownames(S_r) = sapply(1:nrow(S_r), function(i) sprintf("Dim%d", i))
    colMaps(ace)[[sprintf("%s", reduction_slot)]] <- Matrix::t(S_r)
    colMapTypes(ace)[[sprintf("%s", reduction_slot)]] = "reduction"


	V = reduction.out$V
	colnames(V) = sapply(1:dim(V)[2], function(i) sprintf("V%d", i))
	rowMaps(ace)[[sprintf("%s_V", reduction_slot)]] = V
	rowMapTypes(ace)[[sprintf("%s_V", reduction_slot)]] = "internal"


	A = reduction.out$A
	colnames(A) = sapply(1:dim(A)[2], function(i) sprintf("A%d", i))
	rowMaps(ace)[[sprintf("%s_A", reduction_slot)]] = A
	rowMapTypes(ace)[[sprintf("%s_A", reduction_slot)]] = "internal"


	B = reduction.out$B
	colnames(B) = sapply(1:dim(B)[2], function(i) sprintf("B%d", i))
	colMaps(ace)[[sprintf("%s_B", reduction_slot)]] = B
	colMapTypes(ace)[[sprintf("%s_B", reduction_slot)]] = "internal"


	metadata(ace)[[sprintf("%s_sigma", reduction_slot)]] = reduction.out$sigma

	return(ace)
}

#' @export
orthogonalize.ace.batch.simple <- function(ace, batch_attr, reduction_slot = "ACTION") {
  batch_attr = .get_ace_split_IDX(ace, attr = batch_attr, return_split_vec = TRUE) %>% as.factor
	design_mat = model.matrix(~batch_attr)

	ace = orthogonalize.ace.batch(ace, design_mat, reduction_slot = reduction_slot)
	return(ace)

}

#' @export
reduce.and.batch.orthogonalize.ace <- function (ace, design_mat, reduced_dim = 50, max_iter = 10, assay_name = "logcounts", norm_method = "default", reduction_slot = "ACTION", seed = 0, SVD_algorithm = 0) {

	if(!is.matrix(design_mat)) {
		err = sprintf("'design_mat' must be a matrix.\n")
		stop(err)
	}

	ace = reduce.ace(ace, reduced_dim = reduced_dim, max_iter = max_iter, assay_name = assay_name,
    norm_method = norm_method, reduction_slot = reduction_slot, seed = seed, SVD_algorithm = SVD_algorithm)

	ace = orthogonalize.ace.batch(ace, design_mat, reduction_slot = reduction_slot, assay_name = assay_name)
	return(ace)
}
