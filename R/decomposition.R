.decomp_methods <- c("simplex_regression", "spa", "svd", "pca", "action_reduce", "aa", "action", "action_batchcorr")


decomp.simplex <- function(
  X,
  B,
  return_raw = FALSE
){

  if (!is.matrix(X) || !is.matrix(B)) {
    err = sprintf("'X' and 'W0' must be of type 'matrix'.\n")
    warning(err)
  }

  H <- run_simplex_regression(A = as.matrix(X), B = as.matrix(B))

  if(return_raw == TRUE){
    out <-  H
  } else {
    out <- list(W = X, H = H)
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
    warning(err)
    X = as.matrix(X)
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
  max_iter = 50,
  min_delta = 1e-300,
  return_raw = FALSE
){

  if (!is.matrix(X)) {
    err = sprintf("'X' must be of type 'matrix'.\n")
    warning(err)
    X = as.matrix(X)
  }

  if (!is.null(W0) && !is.matrix(W0)) {
    err = sprintf("'W0' must be of type 'matrix' if given.\n")
    warning(err)
    W0 = as.matrix(W0)
  }

  if (is.null(W0)) {
    if (is.null(k)) {
      err = sprintf("'k' cannot be 'NULL' if 'W0' is not given.\n")
      stop(err)
    }
    W0 <- matrix(0, NROW(X), k)
  }

  AA.out <- run_AA(X, W0 = W0, max_it =  max_iter, min_delta =  min_delta)

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
  max_iter = 50,
  min_delta = 1e-300,
  return_raw = FALSE) {

  if (!is.matrix(X)) {
    err = sprintf("'X' must be of type 'matrix'.\n")
    warning(err)
    X = as.matrix(X)
  }
  
  ACTION.out <- run_ACTION(
    S_r = X,
    k_min = k,
    k_max = k,
    thread_no = 1,
    max_it = max_iter,
    min_delta = min_delta
  )

  if(return_raw == TRUE){
    out <-  ACTION.out
  } else {
    out <- list(
      W = ACTION.outt$W,
      H = ACTION.outt$H,
      misc = list(C = ACTION.out$C)
    )
  }
  return(out)
}


decomp.ACTIONMR <- function(
  X,
  k_min,
  k_max,
  specificity_th = -3,
  min_cells_per_arch = 2,
  unification_backbone_density = 0.5,
  unification_resolution = 1.0,
  unification_min_cluster_size = 3,
  max_iter = 50,
  min_delta = 1e-100,
  thread_no = 1,
  return_raw = FALSE) {

  if (!is.matrix(X)) {
    err = sprintf("'X' must be of type 'matrix'.\n")
    warning(err)
    X = as.matrix(X)
  }

  ACTION.out <- run_ACTION(
    S_r = X,
    k_min = k_min,
    k_max = k_max,
    thread_no = thread_no,
    max_it = max_iter,
    min_delta = min_delta
  )

  # Prune nonspecific and/or unreliable archetypes
  pruning.out <- .pruneArchetypes(
    C_trace = ACTION.out$C,
    H_trace = ACTION.out$H,
    specificity_th = specificity_th,
    min_cells_per_arch = min_cells_per_arch
  )

  # Identiy equivalent classes of archetypes and group them together
  unification.out <- .unifyArchetypes(
    S_r = X,
    C_stacked = pruning.out$C_stacked,
    H_stacked = pruning.out$H_stacked,
    unification_backbone_density = unification_backbone_density,
    unification_resolution = unification_resolution,
    unification_min_cluster_size = unification_min_cluster_size,
    thread_no = thread_no
  )

  H <- unification.out$H_unified
  W <- X %*% unification.out$C_unified
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

  return(out)
}

.runFastIRLB <- function(
  X,
  k = 30,
  max_iter = 10,
  seed = 0
){

  X <- .validate_matrix(X)
  max_iter = max_iter * 100
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
  unification_backbone_density = 0.5,
  unification_resolution = 1.0,
  unification_min_cluster_size = 3,
  thread_no = 0) {

  if (!all(dim(C_stacked) == rev(dim(H_stacked)))) {
    err = sprintf("Dimensions of `C_stacked` and transposed `H_stacked` do not match.\n")
    stop(err)
  }


  out <- unify_archetypes(
    S_r = S_r,
    C_stacked = C_stacked,
    H_stacked = H_stacked,
    backbone_density = unification_backbone_density,
    resolution = unification_resolution,
    min_cluster_size = unification_min_cluster_size,
    thread_no = thread_no
  )

  return(out)
}
