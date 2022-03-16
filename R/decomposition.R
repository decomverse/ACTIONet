#' @export
decomp <- function(X, k, method, W0 = NULL, H0 = NULL, params) {
  if (method == "simplex_regression") {
    W <- as.matrix(W0)
    B <- as.matrix(X)
    H <- run_simplex_regression(A = W, B = B)
    extra <- list()
  } else if (method == "SVD") {

    out <- decomp.SVD(
      X = X,
      reduced_dim = k,
      max_iter = params$max_iter,
      seed = params$seed
    )
    return(out)

  } else if (method == "PCA") {

    out <- decomp.PCA(
      X,
      reduced_dim = k,
      max_iter = params$max_iter,
      seed = params$seed,
    )

    return(out)

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

    out <- decomp.ACTION(
      X = X,
      k = k,
      k_max = NULL,
      batch_vec = NULL,
      specificity_th = NULL,
      min_cells_per_arch = NULL,
      unification_th = NULL,
      max_iter = params$max_iter,
      thread_no = params$thread_no,
      algorithm = "default",
      min_delta = 1e-300,
      seed = params$seed
    )

    return(out)

  } else if (method == "ACTION_decomposition_MR") {

    out <- decomp.ACTION(
      X = X,
      k = params$k_min,
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

    return(out)

  } else if (method == "ACTION_decomposition_MR_with_batch_correction") {

    out <- decomp.ACTION(
      X = X,
      k = params$k_min,
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

    return(out)
  }

  out <- list(W = W, H = H, extra = extra)
  return(out)
}

decomp.SVD <- function(
  X,
  reduced_dim = 50,
  max_iter = 10,
  seed = 0
){

  SVD.out <- .runFastIRLB(
    X = X,
    reduced_dim = reduced_dim,
    max_iter = max_iter,
    seed = seed
  )

  W <- SVD.out$u
  H <- as.matrix(Matrix::Diagonal(length(SVD.out$d), SVD.out$d) %*% Matrix::t(SVD.out$v))
  extra <- list(d = SVD.out$d)

  out <- list(W = W, H = H, extra = extra)

  return(out)
}


decomp.PCA <- function(
  X,
  reduced_dim = 50,
  max_iter = 10,
  seed = 0
){

  X <- .checkSparse(X)

  SVD.out <- .runFastIRLB(
    X = X,
    reduced_dim = reduced_dim,
    max_iter = max_iter,
    seed = seed
  )

  if (is.matrix(X)) {
    reduction.out <- SVD2PCA_full(S = X, u = SVD.out$u, d = SVD.out$d, v = SVD.out$v)
  } else {
    reduction.out <- SVD2PCA(S = X, u = SVD.out$u, d = SVD.out$d, v = SVD.out$v)
  }

  W <- reduction.out$x
  H <- reduction.out$rotation
  extra <- list(sdev = reduction.out$sdev)
  extra$d <- extra$sdev * sqrt(nrow(W) - 1)
  out <- list(W = W, H = H, extra = extra)

  return(out)
}


#' Run ACTION_decomposition_MR
decomp.ACTION_MR <- function(
  ace = NULL,
  S_r = NULL,
  batch_attr = NULL,
  k_min = 2,
  k_max = 30,
  specificity_th = -3,
  min_cells_per_arch = 2,
  unification_th = 0,
  max_iter = 50,
  thread_no = 0,
  reduction_slot = "ACTION",
  algorithm = c("default", "batchcorr"),
  return_raw = FALSE,
  seed = 0
) {

  algorithm <- tolower(algorithm)
  algorithm <- match.arg(algorithm, several.ok = FALSE)

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

  if (algorithm == "batchcorr") {
    batch_vec <- .validate_attr(Matrix::t(S_r), batch_attr)
    batch_vec = as.numeric(factor(batch_vec))
    ACTION.out <- run_ACTION_with_batch_correction(
      S_r = S_r,
      batch = batch_vec,
      k_min = k_min,
      k_max = k_max,
      thread_no = thread_no,
      max_it = max_iter,
      min_delta = 1e-300
    )
  } else {
    ACTION.out <- run_ACTION(
      S_r = S_r,
      k_min = k_min,
      k_max = k_max,
      thread_no = thread_no,
      max_it = max_iter,
      min_delta = 1e-300
    )
  }

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
  seed = 0
) {

  algorithm <- tolower(algorithm)
  algorithm <- match.arg(algorithm, several.ok = FALSE)

  S_r <- Matrix::t(scale(X))

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
    batch_vec <- .validate_attr(X, batch_vec)
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
    ACTION.out$H[k_min:k_max]
    ACTION.out$C[k_min:k_max]

    H <- ACTION.out$H[[k]]
    C <- ACTION.out$C[[k]]
    W <- S_r %*% C
    extra <- list(C = C)
  } else {
    # Prune nonspecific and/or unreliable archetypes
    pruning.out <- .pruneArchetypes(
      C_trace = ACTION.out$C,
      H_trace = ACTION.out$H,
      specificity_th = specificity_th,
      min_cells_per_arch = min_cells_per_arch,
    )

    # Identiy equivalent classes of archetypes and group them together
    unification.out <- .unifyArchetypes(
      S_r = S_r,
      C_stacked = pruning.out$C_stacked,
      H_stacked = pruning.out$H_stacked,
      violation_threshold = unification_th,
      thread_no = thread_no,
    )

    H <- unification.out$H_unified
    W <- S_r %*% unification.out$C_unified
    extra <- list(
      H = ACTION.out$H,
      C = ACTION.out$C,
      H_stacked = pruning.out$H_stacked,
      C_stacked = pruning.out$C_stacked,
      H_unified = unification.out$H_unified,
      C_unified = unification.out$C_unified,
      assigned_archetype = unification.out$assigned_archetype
    )
  }

  out <- list(W = W, H = H, extra = extra)
  return(out)
}


.runFastIRLB <- function(
  X,
  reduced_dim = 50,
  max_iter = 10,
  seed = 0
){

  X <- .checkSparse(X)

  if (is.matrix(X)) {
    out <- IRLB_SVD_full(
      A = X,
      dim = reduced_dim,
      iter = max_iter,
      seed = seed
    )
  } else {
    out <- IRLB_SVD(
      A = X,
      dim = reduced_dim,
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
  violation_threshold = 0,
  thread_no = 0
) {

  if (dim(C_stacked) == rev(dim(H_stacked))) {
    err = sprintf("Dimensions of `C_stacked` and transposed `H_stacked` do not match.\n")
    stop(err)
  }

  out <- unify_archetypes(
    S_r = S_r,
    C_stacked = C_stacked,
    H_stacked = H_stacked,
    violation_threshold = violation_threshold,
    thread_no = thread_no
  )

  return(out)
}
