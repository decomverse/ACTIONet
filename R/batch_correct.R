#' Perform batch correction on `ACTIONetExperiment` and `SingleCellExperiment` objects.

#' @export
reduce.and.batch.correct.ace.fastMNN <- function(ace,
                                                 batch_attr = NULL,
                                                 assay_name = NULL,
                                                 reduced_dim = 50,
                                                 MNN_k = 20,
                                                 reduction_out = "MNN",
                                                 BPPARAM = SerialParam()) {
  ACTIONetExperiment:::.check_and_load_package(c("scran", "SingleCellExperiment", "batchelor", "BiocParallel"))

  if (is.null(batch_attr)) {
    warning("'batch_attr' must be provided")
    return(ace)
  }

  ace <- as(ace, "ACTIONetExperiment")
  m_data <- metadata(ace)

  if (is.null(assay_name)) {
    if ("default_assay" %in% names(metadata(ace))) {
      message(sprintf("Input assay_name is NULL. Setting assay_name to the metadata(ace)[['default_assay']]"))
      assay_name <- metadata(ace)[["default_assay"]]
    } else {
      message(sprintf("Input assay_name is NULL. Setting assay_name to logcounts"))
      assay_name <- "logcounts"
    }
  }
  .validate_assay(ace, assay_name = assay_name, return_elem = FALSE)


  S <- SummarizedExperiment::assays(ace)[[assay_name]]
  mnn_batch <- ACTIONetExperiment::get.data.or.split(ace, attr = batch_attr, to_return = "data")
  IDX <- ACTIONetExperiment::get.data.or.split(ace, attr = batch_attr, to_return = "split")
  merge_order <- order(sapply(IDX, function(idx) length(idx)), decreasing = TRUE)

  set.seed(0)
  mnn.out <- batchelor::fastMNN(
    S,
    batch = mnn_batch,
    k = MNN_k,
    d = reduced_dim,
    auto.merge = FALSE,
    merge.order = merge_order,
    cos.norm = FALSE,
    BPPARAM = BPPARAM
  )

  S_r <- SingleCellExperiment::reducedDims(mnn.out)[["corrected"]]
  rownames(S_r) <- colnames(ace)
  colnames(S_r) <- sapply(1:dim(S_r)[2], function(i) sprintf("PC%d", i))

  colMaps(ace)[[reduction_out]] <- S_r
  colMapTypes(ace)[[reduction_out]] <- "reduction"

  metadata(ace)[["default_reduction"]] <- reduction_out
  metadata(ace)[["default_assay"]] <- assay_name


  V <- rowData(mnn.out)[["rotation"]]
  colnames(V) <- paste0("V", 1:NCOL(V))
  rowMaps(ace)[[sprintf("%s_V", reduction_out)]] <- V
  rowMapTypes(ace)[[sprintf("%s_V", reduction_out)]] <- "internal"

  invisible(gc())

  metadata(ace) <- m_data
  return(ace)
}


#' @export
reduce.and.batch.correct.ace.Harmony <- function(ace,
                                                 batch_attr,
                                                 reduced_dim = 50,
                                                 max_iter = 1000,
                                                 assay_name = NULL,
                                                 reduction_out = "ACTION",
                                                 harmony_out = "Harmony",
                                                 seed = 0,
                                                 SVD_algorithm = 0,
                                                 harmony_clustering_algorithm = 2) {
  if (!require(harmony)) {
    err <- sprintf("You need to install harmony (https://github.com/immunogenomics/harmony).\n")
    stop(err)
  }

  if (is.null(batch_attr)) {
    err <- sprintf("'batch_attr' must be provided.\n")
    stop(err)
  }

  ace <- as(ace, "ACTIONetExperiment")

  if (is.null(assay_name)) {
    if ("default_assay" %in% names(metadata(ace))) {
      message(sprintf("Input assay_name is NULL. Setting assay_name to the metadata(ace)[['default_assay']]"))
      assay_name <- metadata(ace)[["default_assay"]]
    } else {
      message(sprintf("Input assay_name is NULL. Setting assay_name to logcounts"))
      assay_name <- "logcounts"
    }
  }
  .validate_assay(ace, assay_name = assay_name, return_elem = FALSE)


  batch_attr <- ACTIONetExperiment::get.data.or.split(ace, attr = batch_attr, to_return = "data")

  ace <- reduce.ace(
    ace = ace,
    reduced_dim = reduced_dim,
    max_iter = max_iter,
    assay_name = assay_name,
    reduction_out = reduction_out,
    seed = seed,
    SVD_algorithm = SVD_algorithm
  )

  metadata(ace)[["default_assay"]] <- assay_name

  ace <- batch.correct.ace.Harmony(
    ace = ace,
    batch_attr = batch_attr,
    reduction_slot = reduction_out,
    harmony_out = harmony_out,
    harmony_clustering_algorithm = harmony_clustering_algorithm
  )

  return(ace)
}

#' @export
batch.correct.ace.Harmony <- function(ace,
                                      batch_attr = NULL,
                                      reduction_slot = "ACTION",
                                      harmony_out = "Harmony",
                                      harmony_clustering_algorithm = 2) {
  if (!require(harmony)) {
    err <- sprintf("You need to install Harmony (https://github.com/immunogenomics/harmony) first for batch-correction.\n")
    stop(err)
  }

  if (is.null(batch_attr)) {
    err <- sprintf("'batch_attr' must be provided.\n")
    stop(err)
  }

  ace <- as(ace, "ACTIONetExperiment")

  batch_vec <- as.numeric(factor(ACTIONetExperiment::get.data.or.split(ace, attr = batch_attr, to_return = "data")))

  X <- Matrix::t(colMaps(ace)[[reduction_slot]])
  X.norm <- normalize_mat(X, 2)
  SPA_out <- run_SPA(X, k = min(ncol(X), nrow(X)))
  W0 <- X.norm[, SPA_out$selected_columns]

  X_corr <- run_harmony(X = X, W0 = W0, batch = batch_vec, clustering_algorithm = harmony_clustering_algorithm)
  colMaps(ace)[[harmony_out]] <- Matrix::t(X_corr)

  metadata(ace)[["default_reduction"]] <- harmony_out

  return(ace)
}


#' @export
orthogonalize.ace.batch <- function(ace,
                                    design_mat,
                                    reduction_slot = "ACTION",
                                    ortho_out = "ACTION_ortho",
                                    assay_name = NULL) {
  ace <- as(ace, "ACTIONetExperiment")

  if (is.null(assay_name)) {
    if ("default_assay" %in% names(metadata(ace))) {
      message(sprintf("Input assay_name is NULL. Setting assay_name to the metadata(ace)[['default_assay']]"))
      assay_name <- metadata(ace)[["default_assay"]]
    } else {
      message(sprintf("Input assay_name is NULL. Setting assay_name to logcounts"))
      assay_name <- "logcounts"
    }
  }
  .validate_assay(ace, assay_name = assay_name, return_elem = FALSE)


  S <- SummarizedExperiment::assays(ace)[[assay_name]]
  S_r <- colMaps(ace)[[sprintf("%s", reduction_slot)]]
  V <- rowMaps(ace)[[sprintf("%s_V", reduction_slot)]]
  A <- rowMaps(ace)[[sprintf("%s_A", reduction_slot)]]
  B <- colMaps(ace)[[sprintf("%s_B", reduction_slot)]]
  sigma <- S4Vectors::metadata(ace)[[sprintf("%s_sigma", reduction_slot)]]

  if (is.matrix(S)) {
    reduction.out <- orthogonalize_batch_effect_full(
      S = S,
      old_S_r = S_r,
      old_V = V,
      old_A = A,
      old_B = B,
      old_sigma = sigma,
      design = design_mat
    )
  } else {
    reduction.out <- orthogonalize_batch_effect(
      S = S,
      old_S_r = S_r,
      old_V = V,
      old_A = A,
      old_B = B,
      old_sigma = sigma,
      design = design_mat
    )
  }
  S_r <- reduction.out$S_r
  colnames(S_r) <- colnames(ace)
  rownames(S_r) <- sapply(1:nrow(S_r), function(i) sprintf("Dim%d", i))
  colMaps(ace)[[sprintf("%s", ortho_out)]] <- Matrix::t(S_r)
  colMapTypes(ace)[[sprintf("%s", ortho_out)]] <- "reduction"


  V <- reduction.out$V
  colnames(V) <- sapply(1:dim(V)[2], function(i) sprintf("V%d", i))
  rowMaps(ace)[[sprintf("%s_V", ortho_out)]] <- V
  rowMapTypes(ace)[[sprintf("%s_V", ortho_out)]] <- "internal"


  A <- reduction.out$A
  colnames(A) <- sapply(1:dim(A)[2], function(i) sprintf("A%d", i))
  rowMaps(ace)[[sprintf("%s_A", ortho_out)]] <- A
  rowMapTypes(ace)[[sprintf("%s_A", ortho_out)]] <- "internal"


  B <- reduction.out$B
  colnames(B) <- sapply(1:dim(B)[2], function(i) sprintf("B%d", i))
  colMaps(ace)[[sprintf("%s_B", ortho_out)]] <- B
  colMapTypes(ace)[[sprintf("%s_B", ortho_out)]] <- "internal"


  S4Vectors::metadata(ace)[[sprintf("%s_sigma", ortho_out)]] <- reduction.out$sigma

  return(ace)
}


#' @export
orthogonalize.ace.batch.simple <- function(ace,
                                           batch_attr,
                                           reduction_slot = "ACTION",
                                           ortho_out = "ACTION_ortho",
                                           assay_name = NULL) {
  ace <- as(ace, "ACTIONetExperiment")

  if (is.null(assay_name)) {
    if ("default_assay" %in% names(metadata(ace))) {
      message(sprintf("Input assay_name is NULL. Setting assay_name to the metadata(ace)[['default_assay']]"))
      assay_name <- metadata(ace)[["default_assay"]]
    } else {
      message(sprintf("Input assay_name is NULL. Setting assay_name to logcounts"))
      assay_name <- "logcounts"
    }
  }
  .validate_assay(ace, assay_name = assay_name, return_elem = FALSE)


  batch_attr <- ACTIONetExperiment::get.data.or.split(ace, attr = batch_attr, to_return = "data")
  batch_attr <- as.factor(batch_attr)
  design_mat <- stats::model.matrix(~batch_attr)

  ace <- orthogonalize.ace.batch(
    ace,
    design_mat,
    reduction_slot = reduction_slot,
    ortho_out = ortho_out,
    assay_name = assay_name
  )

  return(ace)
}

#' @export
reduce.and.batch.orthogonalize.ace <- function(ace,
                                               design_mat,
                                               reduced_dim = 50,
                                               max_iter = 1000,
                                               assay_name = NULL,
                                               reduction_out = "ACTION",
                                               ortho_out = "ACTION_ortho",
                                               seed = 0,
                                               SVD_algorithm = 0) {
  ace <- as(ace, "ACTIONetExperiment")

  if (is.null(assay_name)) {
    if ("default_assay" %in% names(metadata(ace))) {
      message(sprintf("Input assay_name is NULL. Setting assay_name to the metadata(ace)[['default_assay']]"))
      assay_name <- metadata(ace)[["default_assay"]]
    } else {
      message(sprintf("Input assay_name is NULL. Setting assay_name to logcounts"))
      assay_name <- "logcounts"
    }
  }
  .validate_assay(ace, assay_name = assay_name, return_elem = FALSE)

  if (!is.matrix(design_mat)) {
    err <- sprintf("'design_mat' must be a matrix.\n")
    stop(err)
  }

  ace <- reduce.ace(
    ace = ace,
    reduced_dim = reduced_dim,
    max_iter = max_iter,
    assay_name = assay_name,
    reduction_out = reduction_out,
    seed = seed,
    SVD_algorithm = SVD_algorithm
  )

  ace <- orthogonalize.ace.batch(
    ace = ace,
    design_mat = design_mat,
    reduction_slot = reduction_out,
    ortho_out = ortho_out,
    assay_name = assay_name
  )

  return(ace)
}