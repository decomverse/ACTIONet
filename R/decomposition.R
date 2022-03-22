.decomp_methods <- c("simplex_regression", "spa", "svd", "pca", "action_reduce", "aa", "action", "action_batchcorr")

#' @export
decomp <- function(X, method, k = NULL, ...) {

  method <- tolower(method)
  method = match.arg(arg = method, choices = .decomp_methods, several.ok = FALSE)

  params = list(...)

  if (method == "simplex_regression") {

    out <- decomp.simplex(
      X = X,
      W0 = params$W0
    )

  } else if (method == "spa") {

    out <- decomp.SPA(
      X = X,
      k = k
    )

  } else if (method == "svd") {

    out <- decomp.SVD(
      X = X,
      k = k,
      max_iter = params$max_iter,
      seed = params$seed
    )

  } else if (method == "pca") {

    out <- decomp.PCA(
      X = X,
      k = k,
      max_iter = params$max_iter,
      seed = params$seed
    )

  } else if (method == "action_reduce") {

    out <- decomp.ACTION_red(
      X = X,
      k = k,
      max_iter = params$max_iter,
      seed = params$seed
    )

  } else if (method == "aa") {

    out <- decomp.AA(
      X = X,
      k = k,
      W0 = params$W0
    )

  } else if (method == "action") {

    out <- decomp.ACTION(
      X = X,
      k = k,
      k_max = params$k_max,
      batch_vec = NULL,
      specificity_th = params$specificity_th,
      min_cells_per_arch = params$min_cells_per_arch,
      unification_th = params$unification_th,
      max_iter = params$max_iter,
      thread_no = params$thread_no,
      algorithm = "default",
      min_delta = 1e-300,
      seed = params$seed
    )

  } else if (method == "action_batchcorr") {

    out <- decomp.ACTION(
      X = X,
      k = k,
      k_max = params$k_max,
      batch_vec = params$batch,,
      specificity_th = params$specificity_th,
      min_cells_per_arch = params$min_cells_per_arch,
      unification_th = params$unification_th,
      max_iter = params$max_iter,
      thread_no = params$thread_no,
      algorithm = "batchcorr",
      min_delta = 1e-300,
      seed = params$seed
    )

  }
  return(out)
}


decomp.simplex <- function(
  X,
  W0,
  return_raw = FALSE
){

  if (!is.matrix(X) || !is.matrix(W0)) {
    err = sprintf("'X' and 'W0' must be of type 'matrix'.\n")
    stop(err)
  }

  W <- as.matrix(W0)
  H <- run_simplex_regression(A = W, B = as.matrix(X))

  if(return_raw == TRUE){
    out <-  H
  } else {
    out <- list(W = W, H = H)
  }
  return(out)
}


decomp.SPA <- function(
  X,
  k = 30,
  return_raw = FALSE
){

  if (!is.matrix(X)) {
    err = sprintf("'X' must be of type 'matrix'.\n")
    stop(err)
  }

  SPA.out <- run_SPA(X, k)

  if(return_raw == TRUE){
    out <-  SPA.out
  } else {
    W <- X[, SPA.out$selected_columns]
    H <- run_simplex_regression(W, X)
    misc <- list(selected_cols = SPA.out$selected_columns, norms = SPA.out$norms)

    out <- list(W = W, H = H, misc = misc)
  }
  return(out)
}


decomp.SVD <- function(
  X,
  k = 30,
  max_iter = 10,
  seed = 0,
  return_raw = FALSE
){

  SVD.out <- .runFastIRLB(
    X = X,
    k = k,
    max_iter = max_iter,
    seed = seed
  )

  if(return_raw == TRUE){
    out <- SVD.out
  } else {
    W <- SVD.out$u
    H <- as.matrix(Matrix::Diagonal(length(SVD.out$d), SVD.out$d) %*% Matrix::t(SVD.out$v))
    misc <- list(d = SVD.out$d)

    out <- list(W = W, H = H, misc = misc)
  }
  return(out)
}


decomp.PCA <- function(
  X,
  k = 30,
  max_iter = 10,
  seed = 0,
  return_raw = FALSE
){

  X <- .validate_matrix(X)

  SVD.out <- .runFastIRLB(
    X = X,
    k = k,
    max_iter = max_iter,
    seed = seed
  )

  if (is.matrix(X)) {
    reduction.out <- SVD2PCA_full(S = X, u = SVD.out$u, d = SVD.out$d, v = SVD.out$v)
  } else {
    reduction.out <- SVD2PCA(S = X, u = SVD.out$u, d = SVD.out$d, v = SVD.out$v)
  }

  if(return_raw == TRUE){
    out <- reduction.out
  } else {
    W <- reduction.out$x
    H <- reduction.out$rotation
    misc <- list(sdev = reduction.out$sdev)
    misc$d <- misc$sdev * sqrt(nrow(W) - 1)

    out <- list(W = W, H = H, misc = misc)
  }
  return(out)
}


decomp.ACTION_red <- function(
  X,
  k = 30,
  max_iter = 10,
  seed = 0,
  return_raw = FALSE
){

  X <- .validate_matrix(X)

  if (is.matrix(X)) {
    reduction.out <- reduce_kernel_full(
      S = X, reduced_dim = k,
      iter = max_iter, seed = seed, SVD_algorithm = 0
    )
  } else {
    reduction.out <- reduce_kernel(
      S = X, reduced_dim = k,
      iter = max_iter, seed = seed, SVD_algorithm = 0
    )
  }

  if(return_raw == TRUE){
    out <- reduction.out
  } else {
    W <- reduction.out$V
    H <- reduction.out$S_r
    misc <- list(A = reduction.out$A, B = reduction.out$B, reduction.out$sigma)

    out <- list(W = W, H = H, misc = misc)
  }
  return(out)
}


decomp.AA <- function(
  X,
  k = NULL,
  W0 = NULL,
  return_raw = FALSE
){

  if (!is.matrix(X)) {
    err = sprintf("'X' must be of type 'matrix'.\n")
    stop(err)
  }

  if (!is.null(W0) && !is.matrix(W0)) {
    err = sprintf("'W0' must be of type 'matrix' if given.\n")
    stop(err)
  }

  if (is.null(W0)) {
    if (is.null(k)) {
      err = sprintf("'k' cannot be 'NULL' if 'W0' is not given.\n")
      stop(err)
    }
    W0 <- matrix(0, NROW(X), k)
  }

  AA.out <- run_AA(X, W0 = W0)

  if(return_raw == TRUE){
    out <-  AA.out
  } else {
    out <- list(
      W = AA.out$W,
      H = AA.out$H,
      misc = list(C = AA.out$C)
    )
  }
  return(out)
}


#' Run ACTION_decomposition
decomp.ACTION <- function(
  X,
  k,
  k_max = NULL,
  batch_vec = NULL,
  specificity_th = -3,
  min_cells_per_arch = 2,
  unification_th = 0,
  max_iter = 50,
  thread_no = 0,
  algorithm = c("default", "batchcorr"),
  min_delta = 1e-300,
  seed = 0,
  return_raw = FALSE
) {

  algorithm <- tolower(algorithm)
  algorithm <- match.arg(algorithm, several.ok = FALSE)

  # S_r <- Matrix::t(scale(X))
  S_r <- Matrix::t(X)

  if(is.null(k_max) || k_max == k) {
    k_max = k
    do_MR = FALSE
  } else if (k_max > k) {
    do_MR = TRUE
  } else {
    err = sprintf("`k_max` must be must be >= `k`/`k_min`.\n")
    stop(err)
  }

  if (algorithm == "batchcorr") {
    batch_vec <- .validate_attr(X, attr = batch_vec, obj_name = "X", attr_name = "batch_vec")
    batch_vec = as.numeric(factor(batch_vec))
    ACTION.out <- run_ACTION_with_batch_correction(
      S_r = S_r,
      batch = batch_vec,
      k_min = k,
      k_max = k_max,
      thread_no = thread_no,
      max_it = max_iter,
      min_delta = min_delta
    )
  } else {
    ACTION.out <- run_ACTION(
      S_r = S_r,
      k_min = k,
      k_max = k_max,
      thread_no = thread_no,
      max_it = max_iter,
      min_delta = min_delta
    )
  }

  if (!do_MR) {
    if(return_raw == TRUE){
      out <- ACTION.out
    } else {
      H <- ACTION.out$H[[k]]
      C <- ACTION.out$C[[k]]
      W <- S_r %*% C
      misc <- list(C = C)
      out <- list(W = W, H = H, misc = misc)
    }
  } else {
    # Prune nonspecific and/or unreliable archetypes
    pruning.out <- .pruneArchetypes(
      C_trace = ACTION.out$C,
      H_trace = ACTION.out$H,
      specificity_th = specificity_th,
      min_cells_per_arch = min_cells_per_arch
    )

    # Identiy equivalent classes of archetypes and group them together
    unification.out <- .unifyArchetypes(
      S_r = S_r,
      C_stacked = pruning.out$C_stacked,
      H_stacked = pruning.out$H_stacked,
      unification_th = unification_th,
      thread_no = thread_no
    )

    if(return_raw == TRUE){
      out <- list(
        ACTION = ACTION.out,
        pruning = pruning.out,
        unification = unification.out
      )
    } else {
      H <- unification.out$H_unified
      W <- S_r %*% unification.out$C_unified
      misc <- list(
        H = ACTION.out$H,
        C = ACTION.out$C,
        H_stacked = pruning.out$H_stacked,
        C_stacked = pruning.out$C_stacked,
        H_unified = unification.out$H_unified,
        C_unified = unification.out$C_unified,
        assigned_archetype = unification.out$assigned_archetype
      )
      out <- list(W = W, H = H, misc = misc)
    }
  }
  return(out)
}


.runFastIRLB <- function(
  X,
  k = 30,
  max_iter = 10,
  seed = 0
){

  X <- .validate_matrix(X)

  if (is.matrix(X)) {
    out <- IRLB_SVD_full(
      A = X,
      dim = k,
      iter = max_iter,
      seed = seed
    )
  } else {
    out <- IRLB_SVD(
      A = X,
      dim = k,
      iter = max_iter,
      seed = seed
    )
  }
  return(out)
}


.pruneArchetypes <- function(
  C_trace,
  H_trace,
  specificity_th = -3,
  min_cells_per_arch = 2
){

  out <- prune_archetypes(
    C_trace = C_trace,
    H_trace = H_trace,
    min_specificity_z_thresh = specificity_th,
    min_cells = min_cells_per_arch
  )

  return(out)
}


.unifyArchetypes <- function(
  S_r,
  C_stacked,
  H_stacked,
  unification_th = 0,
  thread_no = 0
) {

  if (!all(dim(C_stacked) == rev(dim(H_stacked)))) {
    err = sprintf("Dimensions of `C_stacked` and transposed `H_stacked` do not match.\n")
    stop(err)
  }

  out <- unify_archetypes(
    S_r = S_r,
    C_stacked = C_stacked,
    H_stacked = H_stacked,
    violation_threshold = unification_th,
    thread_no = thread_no
  )

  return(out)
}
