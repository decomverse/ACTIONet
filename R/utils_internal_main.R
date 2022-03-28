#' Run ACTION_decomposition_MR
.run.ACTION_MR.ace <- function(
  ace,
  S_r = NULL,
  k_min = 2,
  k_max = 30,
  specificity_th = -3,
  min_cells_per_arch = 2,
  unification_th = 0,
  max_iter = 50,
  thread_no = 0,
  unified_suffix = "unified",
  footprint_slot_name = "assigned_archetype",
  reduction_slot = "ACTION",
  seed = 0,
  return_raw = FALSE
) {

  .validate_ace(ace, allow_null = FALSE, return_elem = FALSE)

  if (is.null(S_r)) {
    S_r <- .validate_map(ace, map_slot = reduction_slot, ace_name = "ace")
  } else if (NCOL(ace) != NROW(S_r)) {
    err = sprint("'NCOL(ace)' must match 'NROW(S_r)'.\n")
    stop(err)
  }

  out <- decomp.ACTION(
    X = S_r,
    k = k_min,
    k_max = k_max,
    batch_vec = NULL,
    specificity_th = specificity_th,
    min_cells_per_arch = min_cells_per_arch,
    unification_th = unification_th,
    max_iter = max_iter,
    thread_no = thread_no,
    algorithm = "default",
    min_delta = 1e-300,
    seed = seed,
    return_raw = TRUE
  )

  if(return_raw == TRUE) {
    return(out)
  } else {
    colMaps(ace)[["H_stacked"]] <- Matrix::t(as(out$pruning$H_stacked, "sparseMatrix"))
    colMapTypes(ace)[["H_stacked"]] <- "internal"

    colMaps(ace)[["C_stacked"]] <- as(out$pruning$C_stacked, "sparseMatrix")
    colMapTypes(ace)[["C_stacked"]] <- "internal"

    colMaps(ace)[[sprintf("H_%s", unified_suffix)]] <- as(Matrix::t(out$unification$H_unified), "sparseMatrix")
    colMapTypes(ace)[[sprintf("H_%s", unified_suffix)]] <- "internal"

    colMaps(ace)[[sprintf("C_%s", unified_suffix)]] <- as(out$unification$C_unified, "sparseMatrix")
    colMapTypes(ace)[[sprintf("C_%s", unified_suffix)]] <- "internal"

    colData(ace)[[footprint_slot_name]] <- c(out$unification$assigned_archetype)

    return(ace)
  }
}


.run.layoutNetwork <- function(
  ace,
  compactness_level = 50,
  n_epochs = 1000,
  algorithm = c("tumap", "umap"),
  init_coor_slot = "ACTION",
  net_slot = "ACTIONet",
  thread_no = 0,
  seed = 0,
  return_raw = FALSE
) {

  algorithm = tolower(algorithm)
  algorithm <- match.arg(algorithm, several.ok = FALSE)

  .validate_ace(ace, allow_null = FALSE, return_elem = FALSE)

  G <- .validate_net(
    ace,
    net_slot = net_slot,
    matrix_type = "sparse",
    force_type = TRUE,
  )

  initial_coordinates <- .validate_map(
    ace,
    map_slot = init_coor_slot,
    matrix_type = "dense",
    force_type = TRUE,
  )
  initial_coordinates <- Matrix::t(scale(initial_coordinates))


  vis.out <- layoutNetwork(
    G = G,
    initial_position = initial_coordinates,
    algorithm = algorithm,
    compactness_level = compactness_level,
    n_epochs = n_epochs,
    thread_no = thread_no,
    seed = seed
  )

  if (return_raw == TRUE) {
    return(vis.out)
  } else {
    colMaps(ace)[["ACTIONred"]] <- Matrix::t(initial_coordinates[1:3, ])
    colMapTypes(ace)[["ACTIONred"]] <- "embedding"

    X <- vis.out$coordinates
    colnames(X) <- c("x", "y")
    rownames(X) <- colnames(ace)
    colMaps(ace)$ACTIONet2D <- X
    colMapTypes(ace)[["ACTIONet2D"]] <- "embedding"

    X <- vis.out$coordinates_3D
    colnames(X) <- c("x", "y", "z")
    rownames(X) <- colnames(ace)
    colMaps(ace)$ACTIONet3D <- X
    colMapTypes(ace)[["ACTIONet3D"]] <- "embedding"

    X <- vis.out$colors
    colnames(X) <- c("r", "g", "b")
    rownames(X) <- colnames(ace)
    colMaps(ace)$denovo_color <- X
    colMapTypes(ace)[["denovo_color"]] <- "embedding"

    return(ace)
  }
}


#' Prune nonspecific and/or unreliable archetypes
.run.pruneArchetypes <- function(
  ace,
  C_trace,
  H_trace,
  specificity_th = -3,
  min_cells_per_arch = 2
) {

  .validate_ace(ace, allow_null = FALSE, return_elem = FALSE)

  pruning.out <- .pruneArchetypes(
    C_trace = C_trace,
    H_trace = H_trace,
    specificity_th = specificity_th,
    min_cells_per_arch = min_cells_per_arch
  )

  colMaps(ace)[["H_stacked"]] <- Matrix::t(as(pruning.out$H_stacked, "sparseMatrix"))
  colMapTypes(ace)[["H_stacked"]] <- "internal"

  colMaps(ace)[["C_stacked"]] <- as(pruning.out$C_stacked, "sparseMatrix")
  colMapTypes(ace)[["C_stacked"]] <- "internal"

  return(ace)
}


#' Identiy equivalent classes of archetypes and group them together
.run.unifyArchetypes <- function(
  ace,
  unification_th = 0,
  reduction_slot = "ACTION",
  C_stacked_slot = "C_stacked",
  H_stacked_slot = "H_stacked",
  unified_suffix = "unified",
  footprint_slot_name = "assigned_archetype",
  thread_no = 0,
  return_raw = FALSE
) {

  .validate_ace(ace, allow_null = FALSE, return_elem = FALSE)

  S_r <- .validate_map(
    ace = ace,
    map_slot = reduction_slot,
    matrix_type = "dense",
    force_type = TRUE
  )
  S_r = Matrix::t(S_r)

  C_stacked <- .validate_map(
    ace = ace,
    map_slot = C_stacked_slot,
    matrix_type = "dense",
    force_type = TRUE
  )

  H_stacked <- .validate_map(
    ace = ace,
    map_slot = H_stacked_slot,
    matrix_type = "dense",
    force_type = TRUE
  )
  H_stacked <- Matrix::t(H_stacked)

  unification.out <- .unifyArchetypes(
    S_r = S_r,
    C_stacked = C_stacked,
    H_stacked = H_stacked,
    unification_th = unification_th,
    thread_no = thread_no
  )

  if (return_raw == TRUE) {
    return(unification.out)
  } else {
    Ht_unified <- as(Matrix::t(unification.out$H_unified), "sparseMatrix")
    colMaps(ace)[[sprintf("H_%s", unified_suffix)]] <- Ht_unified
    colMapTypes(ace)[[sprintf("H_%s", unified_suffix)]] <- "internal"

    colMaps(ace)[[sprintf("C_%s", unified_suffix)]] <- as(unification.out$C_unified, "sparseMatrix")
    colMapTypes(ace)[[sprintf("C_%s", unified_suffix)]] <- "internal"

    colData(ace)[[footprint_slot_name]] <- c(unification.out$assigned_archetype)

    return(ace)
  }
}


.smoothPCs <- function(
  ace,
  diffusion_algorithm = "pagerank",
  alpha = 0.9,
  diffusion_it = 5,
  reduction_slot = "ACTION",
  net_slot = "ACTIONet",
  thread_no = 0,
  return_raw = FALSE
) {

  S_r <- .validate_map(
    ace,
    map_slot = reduction_slot,
    matrix_type = "dense",
    force_type = TRUE,
  )

  G <- .validate_net(
    ace,
    net_slot = net_slot,
    matrix_type = "sparse",
    force_type = TRUE,
  )

  vars <- list(
    V = rowMaps(ace)[[sprintf("%s_V", reduction_slot)]],
    A = rowMaps(ace)[[sprintf("%s_A", reduction_slot)]],
    B = colMaps(ace)[[sprintf("%s_B", reduction_slot)]],
    sigma = S4Vectors::metadata(ace)[[sprintf("%s_sigma", reduction_slot)]]
  )

  if (any(sapply(vars, is.null))) {
    nullvars <- paste(names(vars)[which(sapply(vars, is.null))], collapse = ",")
    err <- sprintf("'%s' missing from 'ace'. Perhaps re-run 'reduce.ace()'.\n", nullvars)
    stop(err)
  }

  V <- vars$V
  A <- vars$A
  B <- vars$B
  sigma <- vars$sigma

  U <- as.matrix(S_r %*% Matrix::Diagonal(length(sigma), 1 / sigma))
  SVD.out <- ACTIONet::perturbedSVD(V, sigma, U, -A, B)
  V.smooth <- networkDiffusion(
    obj = G,
    scores = SVD.out$v,
    algorithm = diffusion_algorithm,
    alpha = alpha,
    thread_no = thread_no,
    max_it = diffusion_it
  )

  H <- V.smooth %*% diag(SVD.out$d)

  if (return_raw == TRUE) {
    out <- list(U = U, SVD.out = SVD.out, V.smooth = V.smooth, H = H)
    return(out)
  } else {
    W <- SVD.out$u
    rownames(W) <- rownames(ace)
    smooth_red_name = sprintf("%s_smooth", reduction_slot)
    smooth_U_name = sprintf("%s_U", reduction_slot)
    rowMaps(ace)[[smooth_U_name]] <- W
    colMaps(ace)[[smooth_red_name]] <- H
    rowMapTypes(ace)[[smooth_U_name]] <- colMapTypes(ace)[[smooth_red_name]] <- "internal"
    return(ace)
  }
}
