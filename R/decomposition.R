#' @export
decomp <- function(X, k, method, W0 = NULL, H0 = NULL, params) {
  if (method == "simplex_regression") {
    W <- as.matrix(W0)
    B <- as.matrix(X)
    H <- run_simplex_regression(A = W, B = B)
    extra <- list()
  } else if (method == "SVD") {
    # if (is.matrix(X)) {
    #   reduction.out <- IRLB_SVD_full(
    #     A = X, dim = k,
    #     iter = params$max_iter, seed = params$seed
    #   )
    # } else {
    #   reduction.out <- IRLB_SVD(
    #     A = X, dim = k,
    #     iter = params$max_iter, seed = params$seed
    #   )
    # }
    # W <- reduction.out$u
    # H <- as.matrix(Diagonal(length(reduction.out$d), reduction.out$d) %*% t(reduction.out$v))
    # extra <- list(d = reduction.out$d)
  } else if (method == "PCA") {
    SVD.out <- decomp(X, k, "SVD", params = params)
    u <- as.matrix(SVD.out$W)
    v <- Matrix::t(as.matrix(Diagonal(length(SVD.out$extra$d), 1 / SVD.out$extra$d) %*% SVD.out$H))
    d <- as.matrix(SVD.out$extra$d)
    if (is.matrix(X)) {
      reduction.out <- SVD2PCA_full(S = X, u = u, d = d, v = v)
    } else {
      reduction.out <- SVD2PCA(S = X, u = u, d = d, v = v)
    }
    W <- reduction.out$x
    H <- reduction.out$rotation
    extra <- list(sdev = reduction.out$sdev)
    extra$d <- extra$sdev * sqrt(nrow(W) - 1)
  } else if (method == "ACTION_reduction") {
    if (is.matrix(X)) {
      reduction.out <- reduce_kernel_full(
        S = X, reduced_dim = k,
        iter = params$max_iter, seed = params$seed, SVD_algorithm = params$SVD_algorithm
      )
    } else {
      reduction.out <- reduce_kernel(
        S = X, reduced_dim = k,
        iter = params$max_iter, seed = params$seed, SVD_algorithm = params$SVD_algorithm
      )
    }
    W <- reduction.out$V
    H <- reduction.out$S_r
    extra <- list(A = reduction.out$A, B = reduction.out$B, reduction.out$sigma)
  } else if (method == "SPA") {
    out <- run_SPA(X, k)
    W <- X[, out$selected_columns]
    H <- run_simplex_regression(W, X)

    extra <- list(selected_columns = out$selected_columns, norms = out$norms)
  } else if (method == "AA") {
    if (is.null(W0)) {
      W0 <- matrix(0, nrow(X), k)
    }
    out <- run_AA(as.matrix(X), W0 = as.matrix(W0))
    W <- out$W
    H <- out$H
    extra <- list(C = out$C)
  } else if (method == "ACTION_decomposition") {
    out <- run_ACTION(S_r = X, k_min = k, k_max = k, thread_no = params$thread_no, max_it = params$max_it, min_delta = 1e-300)
    H <- out$H[[k]]
    C <- out$C[[k]]
    W <- X %*% C
    extra <- list(C = C)
  } else if (method == "ACTION_decomposition_MR") {
    out <- decomp.ACTION_MR(
      ace = NULL,
      S_r = X,
      k_min = params$k_min,
      k_max = params$k_max,
      specificity_th = params$specificity_th,
      min_cells_per_arch = params$min_cells_per_arch,
      unification_th = params$unification_th,
      max_iter = params$max_iter,
      thread_no = params$thread_no,
      reduction_slot = NULL,
      return_raw = TRUE,
      seed = params$seed
    )

    return(out)

  } else if (method == "ACTION_decomposition_MR_with_batch_correction") {
    ACTION.out <- run_ACTION_with_batch_correction(
      S_r = X,
      batch = params$batch,
      k_min = params$k_min,
      k_max = params$k_max,
      thread_no = params$thread_no,
      max_it = params$max_iter,
      min_delta = 1e-300
    )

    # Prune nonspecific and/or unreliable archetypes
    pruning.out <- .run.pruneArchetypes(
      C_trace = ACTION.out$C,
      H_trace = ACTION.out$H,
      ace = NULL,
      specificity_th = params$specificity_th,
      min_cells_per_arch = params$min_cells_per_arch,
      return_raw = TRUE
    )

    # Identiy equivalent classes of archetypes and group them together
    unification.out <- .run.unifyArchetypes(
      ace = NULL,
      S_r = X,
      C_stacked = pruning.out$C_stacked,
      H_stacked = pruning.out$H_stacked,
      reduction_slot = NULL,
      C_stacked_slot = NULL,
      H_stacked_slot = NULL,
      violation_threshold = params$unification_th,
      thread_no = params$thread_no,
      return_raw = TRUE
    )

    H <- unification.out$H_unified
    W <- X %*% unification.out$C_unified
    extra <- list(H = ACTION.out$H, C = ACTION.out$C, H_stacked = pruning.out$H_stacked, C_stacked = pruning.out$C_stacked, H_unified = unification.out$H_unified, C_unified = unification.out$C_unified, assigned_archetype = unification.out$assigned_archetype)

    out <- list(W = W, H = H, extra = extra)
    return(out)
  }

  out <- list(W = W, H = H, extra = extra)
  return(out)
}

decomp.SVD <- function(){

}
if (is.matrix(X)) {
  reduction.out <- IRLB_SVD_full(
    A = X, dim = k,
    iter = params$max_iter, seed = params$seed
  )
} else {
  reduction.out <- IRLB_SVD(
    A = X, dim = k,
    iter = params$max_iter, seed = params$seed
  )
}
W <- reduction.out$u
H <- as.matrix(Diagonal(length(reduction.out$d), reduction.out$d) %*% t(reduction.out$v))
extra <- list(d = reduction.out$d)

#' Run ACTION_decomposition_MR
decomp.ACTION_MR <- function(ace = NULL,
                             S_r = NULL,
                             k_min = 2,
                             k_max = 30,
                             specificity_th = -3,
                             min_cells_per_arch = 2,
                             unification_th = 0,
                             max_iter = 50,
                             thread_no = 0,
                             reduction_slot = "ACTION",
                             return_raw = FALSE,
                             seed = 0) {

  if (return_raw == FALSE && is.null(ace)) {
    err <- sprintf("'ace' cannot be null if 'return_raw=FALSE'")
    stop(err)
  }

  ace <- .validate_ace(ace, allow_null = TRUE)

  if (is.null(S_r)) {
    if (!(reduction_slot %in% names(colMaps(ace)))) {
      err <- sprintf("Attribute '%s' is not in 'colMaps'.\n", reduction_slot)
      stop(err)
    }
    S_r <- Matrix::t(scale(colMaps(ace)[[reduction_slot]]))
  }

  ACTION.out <- run_ACTION(
    S_r = S_r,
    k_min = k_min,
    k_max = k_max,
    thread_no = thread_no,
    max_it = max_iter,
    min_delta = 1e-300
  )

  # Prune nonspecific and/or unreliable archetypes
  pruning.out <- .run.pruneArchetypes(
    C_trace = ACTION.out$C,
    H_trace = ACTION.out$H,
    ace = NULL,
    specificity_th = specificity_th,
    min_cells_per_arch = min_cells_per_arch,
    return_raw = TRUE
  )

  # Identiy equivalent classes of archetypes and group them together
  unification.out <- .run.unifyArchetypes(
    ace = NULL,
    S_r = S_r,
    C_stacked = pruning.out$C_stacked,
    H_stacked = pruning.out$H_stacked,
    reduction_slot = NULL,
    C_stacked_slot = NULL,
    H_stacked_slot = NULL,
    violation_threshold = unification_th,
    thread_no = thread_no,
    return_raw = TRUE
  )

  H <- unification.out$H_unified
  W <- S_r %*% unification.out$C_unified
  extra <- list(H = ACTION.out$H, C = ACTION.out$C, H_stacked = pruning.out$H_stacked, C_stacked = pruning.out$C_stacked, H_unified = unification.out$H_unified, C_unified = unification.out$C_unified, assigned_archetype = unification.out$assigned_archetype)

  out <- list(W = W, H = H, extra = extra)

  if (return_raw == TRUE) {
    return(out)
  } else {
    # Store resulting decompositions
    colMaps(ace)[["H_stacked"]] <- Matrix::t(as(out$extra$H_stacked, "sparseMatrix"))
    colMapTypes(ace)[["H_stacked"]] <- "internal"

    colMaps(ace)[["C_stacked"]] <- as(out$extra$C_stacked, "sparseMatrix")
    colMapTypes(ace)[["C_stacked"]] <- "internal"

    H_unified <- as(Matrix::t(out$extra$H_unified), "sparseMatrix")
    colMaps(ace)[["H_unified"]] <- H_unified
    colMapTypes(ace)[["H_unified"]] <- "internal"

    colMaps(ace)[["C_unified"]] <- as(out$extra$C_unified, "sparseMatrix")
    colMapTypes(ace)[["C_unified"]] <- "internal"

    SummarizedExperiment::colData(ace)[["assigned_archetype"]] <- c(out$extra$assigned_archetype)

    return(ace)
  }
}
