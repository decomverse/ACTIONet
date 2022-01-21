
.run.layoutNetwork <- function(ace,
                               G = NULL,
                               initial_coordinates = NULL,
                               compactness_level = 50,
                               n_epochs = 1000,
                               layout_alg = c("tumap", "umap"),
                               thread_no = 0,
                               reduction_slot = "ACTION",
                               net_slot = "ACTIONet",
                               seed = 0,
                               return_raw = FALSE) {
  layout_alg <- match.arg(toupper(layout_alg), choices = c("TUMAP", "UMAP"), several.ok = FALSE)

  if (is.null(G)) {
    if (!(net_slot %in% names(colNets(ace)))) {
      err <- sprintf("Attribute '%s' is not in 'colNets'.\n", net_slot)
      stop(err)
    }
    G <- colNets(ace)[[net_slot]]
  }

  if (is.null(initial_coordinates)) {
    if (!(reduction_slot %in% names(colMaps(ace)))) {
      err <- sprintf("Attribute '%s' is not in 'colMaps'.\n", reduction_slot)
      stop(err)
    }
    initial_coordinates <- Matrix::t(scale(colMaps(ace)[[reduction_slot]]))
  }

  vis.out <- layoutNetwork(
    G = G,
    initial_position = initial_coordinates,
    algorithm = layout_alg,
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

.run.archetypeFeatureSpecificity <- function(ace,
                                             S = NULL,
                                             H = NULL,
                                             assay_name = "logcounts",
                                             footprint_slot = "archetype_footprint",
                                             thread_no = 0,
                                             return_raw = FALSE) {
  if (is.null(S)) {
    if (!(assay_name %in% names(assays(ace)))) {
      err <- sprintf("'S' not given and %s is not an assay of the input %s object.\n", assay_name, class(ace))
      stop(err)
    }
    S <- SummarizedExperiment::assays(ace)[[assay_name]]
  }

  if (is.null(H)) {
    if (!(footprint_slot %in% names(colMaps(ace)))) {
      err <- sprintf("'H' not given and %s is not in 'colMaps'.\n", footprint_slot)
      stop(err)
    }
    H <- Matrix::t(colMaps(ace)[[footprint_slot]])
  }

  if (is.matrix(S)) {
    specificity.out <- compute_archetype_feature_specificity_full(S, H, thread_no)
  } else {
    specificity.out <- compute_archetype_feature_specificity(S, H, thread_no)
  }

  specificity.out <- lapply(specificity.out, function(specificity.scores) {
    rownames(specificity.scores) <- rownames(ace)
    colnames(specificity.scores) <- paste("A", 1:ncol(specificity.scores), sep = "")
    return(specificity.scores)
  })

  if (return_raw == TRUE) {
    return(specificity.out)
  } else {
    rowMaps(ace)[["unified_feature_profile"]] <- specificity.out[["archetypes"]]
    rowMapTypes(ace)[["unified_feature_profile"]] <- "internal"

    rowMaps(ace)[["unified_feature_specificity"]] <- specificity.out[["upper_significance"]]
    rowMapTypes(ace)[["unified_feature_specificity"]] <- "reduction"

    return(ace)
  }
}


#' Prune nonspecific and/or unreliable archetypes
.run.pruneArchetypes <- function(C_trace,
                                 H_trace,
                                 ace = NULL,
                                 specificity_th = -3,
                                 min_cells_per_arch = 2,
                                 return_raw = FALSE) {
  if (return_raw == FALSE && is.null(ace)) {
    err <- sprintf("'ace' cannot be null if 'return_raw=FALSE'")
    stop(err)
  }

  if (!is.null(ace)) {
    if (class(ace) != "ACTIONetExperiment") {
      err <- sprintf("'ace' must be 'ACTIONetExperiment'.\n")
      stop(err)
    }
  }

  pruning.out <- prune_archetypes(
    C_trace = C_trace,
    H_trace = H_trace,
    min_specificity_z_thresh = specificity_th,
    min_cells = min_cells_per_arch
  )

  if (return_raw == TRUE) {
    return(pruning.out)
  } else {
    colMaps(ace)[["H_stacked"]] <- Matrix::t(as(pruning.out$H_stacked, "sparseMatrix"))
    colMapTypes(ace)[["H_stacked"]] <- "internal"

    colMaps(ace)[["C_stacked"]] <- as(pruning.out$C_stacked, "sparseMatrix")
    colMapTypes(ace)[["C_stacked"]] <- "internal"

    return(ace)
  }
}


#' Identiy equivalent classes of archetypes and group them together
.run.unifyArchetypes <- function(ace = NULL,
                                 S_r = NULL,
                                 C_stacked = NULL,
                                 H_stacked = NULL,
                                 reduction_slot = "ACTION",
                                 C_stacked_slot = "C_stacked",
                                 H_stacked_slot = "H_stacked",
                                 unified_suffix = "unified",
                                 violation_threshold = 0,
                                 thread_no = 0,
                                 return_raw = FALSE) {
  if (return_raw == FALSE && is.null(ace)) {
    err <- sprintf("'ace' cannot be null if 'return_raw=FALSE'")
    stop(err)
  }

  if (!is.null(ace)) {
    if (class(ace) != "ACTIONetExperiment") {
      err <- sprintf("'ace' must be 'ACTIONetExperiment'.\n")
      stop(err)
    }
  }

  if (is.null(S_r)) {
    if (!(reduction_slot %in% names(colMaps(ace)))) {
      err <- sprintf("Attribute '%s' is not in 'colMaps'.\n", reduction_slot)
      stop(err)
    }
    S_r <- Matrix::t(colMaps(ace)[[reduction_slot]])
  }

  if (is.null(C_stacked)) {
    if (!(C_stacked_slot %in% names(colMaps(ace)))) {
      err <- sprintf("Attribute '%s' is not in 'colMaps'.\n", C_stacked_slot)
      stop(err)
    }
    C_stacked <- as.matrix(colMaps(ace)[[C_stacked_slot]])
  }

  if (is.null(H_stacked)) {
    if (!(H_stacked_slot %in% names(colMaps(ace)))) {
      err <- sprintf("Attribute '%s' is not in 'colMaps'.\n", H_stacked_slot)
      stop(err)
    }
    H_stacked <- Matrix::t(as.matrix(colMaps(ace)[[H_stacked_slot]]))
  }

  unification.out <- unify_archetypes(
    S_r = S_r,
    C_stacked = C_stacked,
    H_stacked = H_stacked,
    violation_threshold = violation_threshold,
    thread_no = thread_no
  )

  if (return_raw == TRUE) {
    return(unification.out)
  } else {
    Ht_unified <- as(Matrix::t(unification.out$H_unified), "sparseMatrix")
    colMaps(ace)[[sprintf("H_%s", unified_suffix)]] <- Ht_unified
    colMapTypes(ace)[[sprintf("H_%s", unified_suffix)]] <- "internal"

    colMaps(ace)[[sprintf("C_%s", unified_suffix)]] <- as(unification.out$C_unified, "sparseMatrix")
    colMapTypes(ace)[[printf("C_%s", unified_suffix)]] <- "internal"

    colData(ace)[["assigned_archetype"]] <- c(unification.out$assigned_archetype)

    return(ace)
  }
}


.smoothPCs <- function(ace = NULL,
                       S_r = NULL,
                       V = NULL,
                       A = NULL,
                       B = NULL,
                       sigma = NULL,
                       G = NULL,
                       reduction_slot = "ACTION",
                       net_slot = "ACTIONet",
                       thread_no = 0,
                       return_raw = FALSE) {
  if (return_raw == FALSE && is.null(ace)) {
    err <- sprintf("'ace' cannot be null if 'return_raw=FALSE'")
    stop(err)
  }

  vars <- list(
    V = V,
    A = A,
    B = B,
    sigma = sigma
  )

  if (is.null(S_r)) {
    if (!(reduction_slot %in% names(colMaps(ace)))) {
      err <- sprintf("Attribute '%s' is not in 'colMaps'.\n", reduction_slot)
      stop(err)
    }
    S_r <- colMaps(ace)[[reduction_slot]]
  }

  if (is.null(G)) {
    if (!(net_slot %in% names(colNets(ace)))) {
      err <- sprintf("Attribute '%s' is not in 'colNets'.\n", net_slot)
      stop(err)
    }
    G <- colNets(ace)[[net_slot]]
  }

  if (any(sapply(vars, is.null))) {
    if (is.null(ace)) {
      err <- sprintf("'ace' cannot be 'NULL' if any of 'V','A','B', or 'sigma' are missing.\n")
      stop(err)
    }

    vars$V <- rowMaps(ace)[[sprintf("%s_V", reduction_slot)]]
    vars$A <- rowMaps(ace)[[sprintf("%s_A", reduction_slot)]]
    vars$B <- colMaps(ace)[[sprintf("%s_B", reduction_slot)]]
    vars$sigma <- S4Vectors::metadata(ace)[[sprintf("%s_sigma", reduction_slot)]]

    if (any(sapply(vars, is.null))) {
      nullvars <- paste(names(vars)[which(sapply(vars, is.null))], collapse = ",")
      err <- sprintf("'%s' missing from 'ace'.\n", nullvars)
      stop(err)
    }
  }


  V <- vars$V
  A <- vars$A
  B <- vars$B
  sigma <- vars$sigma

  # V <- rowMaps[[sprintf("%s_V", reduction_slot)]]
  # A <- rowMaps(ace)[[sprintf("%s_A", reduction_slot)]]
  # B <- colMaps(ace)[[sprintf("%s_B", reduction_slot)]]
  # sigma <- S4Vectors::metadata(ace)[[sprintf("%s_sigma", reduction_slot)]]
  U <- as.matrix(S_r %*% Diagonal(length(sigma), 1 / sigma))
  SVD.out <- ACTIONet::perturbedSVD(V, sigma, U, -A, B)
  V.smooth <- networkDiffusion(G = G, scores = SVD.out$v, algorithm = "pagerank", alpha = 0.9, thread_no = thread_no)

  H <- V.smooth %*% diag(SVD.out$d)


  if (return_raw == TRUE) {
    out <- list(U = U, SVD.out = SVD.out, V.smooth = V.smooth, H = H)
    return(out)
  } else {
    W <- SVD.out$u
    rownames(W) <- rownames(ace)
    rowMaps(ace)[["SVD_U"]] <- W
    colMaps(ace)[["SVD_V_smooth"]] <- H
    rowMapTypes(ace)[["SVD_U"]] <- colMapTypes(ace)[["SVD_V_smooth"]] <- "internal"
    return(ace)
  }
}
