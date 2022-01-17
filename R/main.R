
#' Run main ACTIONet pipeline
#'
#' @param ace `ACTIONetExperiment` containing an appropriate reduction in 'reduction_slot'.
#' @param k_min Minimum depth of decompositions. (default=2)
#' @param k_max Maximum depth of decompositions. (default=30)
#' @param assay_name Name of assay to be used. (default='logcounts')
#' @param reduction_slot Slot in colMaps(ace) containing reduced kernel. (default='ACTION')
#' @param net_slot_out  Name of slot in colMaps(ace) to store ACTIONet adjacency matrix. (default='ACTIONet')
#' @param min_cells_per_arch Minimum number of observations required to construct an archetype. (default=2)
#' @param max_iter_ACTION Maximum number of iterations for ACTION algorithm. (default=50)
#' @param min_specificity_z_thresh Defines the stringency of pruning nonspecific archetypes.
#' The larger the value, the more archetypes will be filtered out. (default=-3)
#' @param distance_metric Distance metric with which to compute cell-to-cell similarity during network construction. Options are 'jsd' (Jensen-Shannon divergence), L2-norm ('l2'), and inner product ('ip'). (default='jsd')
#' @param nn_approach Nearest-neighbor algorithm to use for network construction. Options are k-nearest neighbors ('knn') and k*-nearest neighbors ('k*nn'). (default='k*nn')
#' @param network_density Density factor of ACTIONet graph. (default=1)
#' @param mutual_edges_only Whether to enforce edges to be mutually-nearest-neighbors. (default=TRUE)
#' @param imputation_alpha Network diffusion parameter to smooth PCs (S_r). (default=0.9)
#' @param layout_compactness A value between 0-100, indicating the compactness of ACTIONet layout. (default=50)
#' @param layout_epochs Number of epochs for SGD algorithm. (default=1000)
#' @param layout_algorithm Algorithm for computing plot layout. can be "UMAP", "TUMAP" (faster), or "forced_atlas". (default="TUMAP")
#' @param layout_in_parallel Run layout construction using multiple cores. May result in marginally different outputs across runs due to parallelization-induced randomization. (default=TRUE)
#' @param unification_violation_threshold Archetype unification resolution parameter. (default=0)
#' @param footprint_alpha Archetype smoothing parameter. (default=0.85)
#' @param thread_no Number of parallel threads. (default=0)
#' @param seed Seed for random initialization. (default=0)
#'
#' @return \itemize{
#' \item If full_trace='FALSE'(default): ACTIONetExperiment object.
#' \item If full_trace='TRUE': Named list containing an ACTIONetExperiment object and a log of ACTIONet function calls.
#' }
#'
#' @examples
#' ace <- runACTIONet(ace)
#' @export
runACTIONet <- function(ace,
                        k_min = 2,
                        k_max = 30,
                        assay_name = "logcounts",
                        reduction_slot = "ACTION",
                        net_slot_out = "ACTIONet",
                        min_cells_per_arch = 2,
                        max_iter_ACTION = 50,
                        min_specificity_z_thresh = -3,
                        distance_metric = "jsd",
                        nn_approach = "k*nn",
                        network_density = 1,
                        mutual_edges_only = TRUE,
                        imputation_alpha = 0.9,
                        layout_compactness = 50,
                        layout_epochs = 1000,
                        layout_algorithm = c("tumap", "umap"),
                        layout_in_parallel = TRUE,
                        unification_violation_threshold = 0,
                        footprint_alpha = 0.15,
                        thread_no = 0,
                        seed = 0) {
  if (!(assay_name %in% names(assays(ace)))) {
    err <- sprintf("'%s' is not an assay of the input '%s' object.\n", assay_name, class(ace))
    stop(err)
  }
  if (!(reduction_slot %in% names(colMaps(ace)))) {
    err <- sprintf("Attribute '%s' is not in 'colMaps'/'reducedDims' of '%s' object.\n", assay_name, class(ace))
    stop(err)
  }
  layout_algorithm <- match.arg(toupper(layout_algorithm), choices = c("TUMAP", "UMAP"), several.ok = FALSE)

  ace <- as(ace, "ACTIONetExperiment")
  S <- SummarizedExperiment::assays(ace)[[assay_name]]
  S_r <- Matrix::t(colMaps(ace)[[reduction_slot]])

  # # Set parameters for ACTION_decomposition_MR
  # params <- list()
  # params$k_min <- k_min
  # params$k_max <- k_max
  # params$min_specificity_z_thresh <- min_specificity_z_thresh
  # params$min_cells_per_arch <- min_cells_per_arch
  # params$unification_violation_threshold <- unification_violation_threshold
  # params$max_iter <- max_iter_ACTION
  # params$thread_no <- thread_no
  # params$seed <- seed
  #
  # # Run ACTION_decomposition_MR
  # ACTION.out <- ACTIONet::decomp(X = S_r, params = params, method = "ACTION_decomposition_MR")
  #
  # # Store resulting decompositions
  # colMaps(ace)[["H_stacked"]] <- Matrix::t(as(ACTION.out$extra$H_stacked, "sparseMatrix"))
  # colMapTypes(ace)[["H_stacked"]] <- "internal"
  #
  # colMaps(ace)[["C_stacked"]] <- as(ACTION.out$extra$C_stacked, "sparseMatrix")
  # colMapTypes(ace)[["C_stacked"]] <- "internal"
  #
  # H_unified <- as(Matrix::t(ACTION.out$extra$H_unified), "sparseMatrix")
  # colMaps(ace)[["H_unified"]] <- H_unified
  # colMapTypes(ace)[["H_unified"]] <- "internal"
  #
  # colMaps(ace)[["C_unified"]] <- as(ACTION.out$extra$C_unified, "sparseMatrix")
  # colMapTypes(ace)[["C_unified"]] <- "internal"
  #
  # ace$assigned_archetype <- c(ACTION.out$extra$assigned_archetype)

  ace <- decomp.ACTION_MR(
    ace = ace,
    S_r = S_r,
    k_min = k_min,
    k_max = k_max,
    min_specificity_z_thresh = min_specificity_z_thresh,
    min_cells_per_arch = min_cells_per_arch,
    unification_violation_threshold = unification_violation_threshold,
    max_iter = max_iter_ACTION,
    thread_no = thread_no,
    reduction_slot = reduction_slot,
    return_raw = FALSE,
    seed = seed
  )

  # Build ACTIONet
  set.seed(seed)
  H <- as.matrix(Matrix::t(colMaps(ace)[["H_stacked"]]))
  G <- buildNetwork(
    H = H,
    algorithm = nn_approach,
    distance_metric = distance_metric,
    density = network_density,
    thread_no = thread_no,
    mutual_edges_only = mutual_edges_only
  )
  colNets(ace)[[net_slot_out]] <- G

  # Use graph core of global and induced subgraphs to infer centrality/quality of each cell
  ace$node_centrality <- c(compute_archetype_core_centrality(G, ace$assigned_archetype))

  # Smooth PCs (S_r) for ease of future imputation (same as MAGIC algorithm)
  ace <- .smoothPCs(
    ace = ace,
    S_r = t(S_r),
    V = rowMaps(ace)[[sprintf("%s_V", reduction_slot)]],
    A = rowMaps(ace)[[sprintf("%s_A", reduction_slot)]],
    B = colMaps(ace)[[sprintf("%s_B", reduction_slot)]],
    sigma = S4Vectors::metadata(ace)[[sprintf("%s_sigma", reduction_slot)]],
    G = G,
    reduction_slot = NULL,
    net_slot = NULL,
    return_raw = FALSE
  )


  # V <- rowMaps(ace)[[sprintf("%s_V", reduction_slot)]]
  # A <- rowMaps(ace)[[sprintf("%s_A", reduction_slot)]]
  # B <- colMaps(ace)[[sprintf("%s_B", reduction_slot)]]
  # sigma <- S4Vectors::metadata(ace)[[sprintf("%s_sigma", reduction_slot)]]
  # U <- as.matrix(Matrix::t(S_r) %*% Diagonal(length(sigma), 1 / sigma))
  # SVD.out <- ACTIONet::perturbedSVD(V, sigma, U, -A, B)
  # # V_svd <- SVD.out$v
  # V.smooth <- networkDiffusion(G, SVD.out$v)
  #
  # H <- V.smooth %*% diag(SVD.out$d)
  # W <- SVD.out$u
  # rownames(W) <- rownames(ace)
  # rowMaps(ace)[["SVD_U"]] <- W
  # colMaps(ace)[["SVD_V_smooth"]] <- H
  # rowMapTypes(ace)[["SVD_U"]] <- colMapTypes(ace)[["SVD_V_smooth"]] <- "internal"
  #

  # initial_coordinates <- ACTIONetExperiment:::.tscalet(S_r)
  # colMaps(ace)[["ACTIONred"]] <- Matrix::t(initial_coordinates[1:3, ])
  # colMapTypes(ace)[["ACTIONred"]] <- "embedding"

  # Layout ACTIONet. Now it uses the smoothed S_r
  ace <- .run.layoutNetwork(ace,
    G = G,
    initial_coordinates = ACTIONetExperiment:::.tscalet(S_r),
    compactness_level = layout_compactness,
    n_epochs = layout_epochs,
    layout_alg = layout_algorithm,
    thread_no = ifelse(layout_in_parallel, thread_no, 1),
    reduction_slot = NULL,
    net_slot = NULL,
    seed = seed
  )

  # Smooth archetype footprints
  # Ht_unified <- colMaps(ace)[["H_unified"]]
  archetype_footprint <- networkDiffusion(
    G = G,
    algorithm = "pagerank",
    scores = as.matrix(colMaps(ace)[["H_unified"]]),
    thread_no = thread_no,
    alpha = footprint_alpha
  )
  colMaps(ace)$archetype_footprint <- archetype_footprint

  # Compute gene specificity for each archetype
  ace <- .run.archetypeFeatureSpecificity(
    ace = ace,
    S = S,
    H = Matrix::t(archetype_footprint),
    assay_name = NULL,
    footprint_slot = NULL,
    return_raw = FALSE
  )
  # H <- Matrix::t(archetype_footprint)
  # if (is.matrix(S)) {
  #   specificity.out <- compute_archetype_feature_specificity_full(S, H)
  # } else {
  #   specificity.out <- compute_archetype_feature_specificity(S, H)
  # }
  #
  # specificity.out <- lapply(specificity.out, function(specificity.scores) {
  #   rownames(specificity.scores) <- rownames(ace)
  #   colnames(specificity.scores) <- paste("A", 1:ncol(specificity.scores), sep = "")
  #   return(specificity.scores)
  # })
  #
  # rowMaps(ace)[["unified_feature_profile"]] <- specificity.out[["archetypes"]]
  # rowMapTypes(ace)[["unified_feature_profile"]] <- "internal"
  #
  # rowMaps(ace)[["unified_feature_specificity"]] <- specificity.out[["upper_significance"]]
  # rowMapTypes(ace)[["unified_feature_specificity"]] <- "reduction"

  # ace <- constructBackbone(
  #   ace = ace,
  #   network_density = network_density,
  #   mutual_edges_only = mutual_edges_only,
  #   layout_compactness = layout_compactness,
  #   layout_epochs = layout_epochs / 5,
  #   thread_no = 1,
  #   net_slot = net_slot_out
  # )

  return(ace)
}

#' Run main ACTIONet pipeline
#'
#' @param ace `ACTIONetExperiment` containing an appropriate reduction in 'reduction_slot'.
#' @param k_min Minimum depth of decompositions. (default=2)
#' @param k_max Maximum depth of decompositions. (default=30)
#' @param assay_name Name of assay to be used. (default='logcounts')
#' @param reduction_slot Slot in colMaps(ace) containing reduced kernel. (default='ACTION')
#' @param net_slot_out  Name of slot in colMaps(ace) to store ACTIONet adjacency matrix. (default='ACTIONet')
#' @param min_cells_per_arch Minimum number of observations required to construct an archetype. (default=2)
#' @param max_iter_ACTION Maximum number of iterations for ACTION algorithm. (default=50)
#' @param min_specificity_z_thresh Defines the stringency of pruning nonspecific archetypes.
#' The larger the value, the more archetypes will be filtered out. (default=-3)
#' @param distance_metric Distance metric with which to compute cell-to-cell similarity during network construction. Options are 'jsd' (Jensen-Shannon divergence), L2-norm ('l2'), and inner product ('ip'). (default='jsd')
#' @param nn_approach Nearest-neighbor algorithm to use for network construction. Options are k-nearest neighbors ('knn') and k*-nearest neighbors ('k*nn'). (default='k*nn')
#' @param network_density Density factor of ACTIONet graph. (default=1)
#' @param mutual_edges_only Whether to enforce edges to be mutually-nearest-neighbors. (default=TRUE)
#' @param layout_compactness A value between 0-100, indicating the compactness of ACTIONet layout. (default=50)
#' @param layout_epochs Number of epochs for SGD algorithm. (default=1000)
#' @param layout_algorithm Algorithm for computing plot layout. t-UMAP ("tumap") or UMAP ("umap"). Not case sensitive. (default="tumap")
#' @param layout_in_parallel Run layout construction using multiple cores. May result in marginally different outputs across runs due to parallelization-induced randomization. (default=TRUE)
#' @param unification_violation_threshold Archetype unification resolution parameter. (default=0)
#' @param footprint_alpha Archetype smoothing parameter. (default=0.85)
#' @param thread_no Number of parallel threads. (default=0)
#' @param full_trace Return list of all intermediate output. Intended for debugging. (default='FALSE')
#' @param seed Seed for random initialization. (default=0)
#'
#' @return \itemize{
#' \item If full_trace='FALSE'(default): ACTIONetExperiment object.
#' \item If full_trace='TRUE': Named list containing an ACTIONetExperiment object and a log of ACTIONet function calls.
#' }
#'
#' @examples
#' ace <- reduce.ace(ace)
#' ace <- run.ACTIONet(ace)
#' @export
run.ACTIONet <- function(ace,
                         k_min = 2,
                         k_max = 30,
                         assay_name = "logcounts",
                         reduction_slot = "ACTION",
                         net_slot_out = "ACTIONet",
                         min_cells_per_arch = 2,
                         max_iter_ACTION = 50,
                         min_specificity_z_thresh = -3,
                         distance_metric = "jsd",
                         nn_approach = "k*nn",
                         network_density = 1,
                         mutual_edges_only = TRUE,
                         layout_compactness = 50,
                         layout_epochs = 1000,
                         layout_algorithm = c("tumap", "umap"),
                         layout_in_parallel = TRUE,
                         unification_violation_threshold = 0,
                         footprint_alpha = 0.85,
                         thread_no = 0,
                         imputation_alpha = 0.9,
                         full_trace = FALSE,
                         seed = 0) {
  if (!(assay_name %in% names(assays(ace)))) {
    err <- sprintf("'%s' is not an assay of the input '%s' object.\n", assay_name, class(ace))
    stop(err)
  }
  if (!(reduction_slot %in% names(colMaps(ace)))) {
    err <- sprintf("Attribute '%s' is not in 'colMaps'/'reducedDims' of '%s' object.\n", assay_name, class(ace))
    stop(err)
  }
  layout_algorithm <- match.arg(toupper(layout_algorithm), choices = c("TUMAP", "UMAP"), several.ok = FALSE)

  ace <- as(ace, "ACTIONetExperiment")
  S <- SummarizedExperiment::assays(ace)[[assay_name]]
  S_r <- Matrix::t(colMaps(ace)[[reduction_slot]])

  # Run ACTION
  ACTION.out <- run_ACTION(
    S_r = S_r,
    k_min = k_min,
    k_max = k_max,
    thread_no = thread_no,
    max_it = max_iter_ACTION,
    min_delta = 1e-300
  )

  # Prune nonspecific and/or unreliable archetypes
  pruning.out <- .run.pruneArchetypes(
    C_trace = ACTION.out$C,
    H_trace = ACTION.out$H,
    ace = NULL,
    min_specificity_z_thresh = min_specificity_z_thresh,
    min_cells_per_arch = min_cells_per_arch,
    return_raw = TRUE
  )

  colMaps(ace)[["H_stacked"]] <- Matrix::t(as(pruning.out$H_stacked, "sparseMatrix"))
  colMapTypes(ace)[["H_stacked"]] <- "internal"

  colMaps(ace)[["C_stacked"]] <- as(pruning.out$C_stacked, "sparseMatrix")
  colMapTypes(ace)[["C_stacked"]] <- "internal"


  # Build ACTIONet
  set.seed(seed)
  G <- buildNetwork(
    H = pruning.out$H_stacked,
    algorithm = nn_approach,
    distance_metric = distance_metric,
    density = network_density,
    thread_no = thread_no,
    mutual_edges_only = mutual_edges_only
  )
  colNets(ace)[[net_slot_out]] <- G

  # Use graph core of global and induced subgraphs to infer centrality/quality of each cell
  ace$node_centrality <- c(compute_archetype_core_centrality(G, ace$assigned_archetype))

  # Smooth PCs (S_r) for ease of future imputation (same as MAGIC algorithm)
  S_r_norm <- networkDiffusion(G, scores = Matrix::t(S_r), algorithm = "pagerank", alpha = imputation_alpha, max_it = 5, thread_no = thread_no)
  colMaps(ace)[["ACTIONnorm"]] <- S_r_norm
  colMapTypes(ace)[["ACTIONnorm"]] <- "internal"

  # Layout ACTIONet
  ace <- .run.layoutNetwork(ace,
    G = G,
    initial_coordinates = ACTIONetExperiment:::.tscalet(S_r),
    compactness_level = layout_compactness,
    n_epochs = layout_epochs,
    layout_alg = layout_algorithm,
    thread_no = ifelse(layout_in_parallel, thread_no, 1),
    reduction_slot = NULL,
    net_slot = NULL,
    seed = seed
  )

  # Identiy equivalent classes of archetypes and group them together
  ace <- .run.unifyArchetypes(
    ace = ace,
    S_r = S_r,
    C_stacked = pruning.out$C_stacked,
    H_stacked = pruning.out$H_stacked, ,
    reduction_slot = NULL,
    C_stacked_slot = NULL,
    H_stacked_slot = NULL,
    violation_threshold = unification_violation_threshold,
    thread_no = thread_no,
    return_raw = FALSE
  )

  # Smooth archetype footprints
  archetype_footprint <- networkDiffusion(
    G = G,
    algorithm = "pagerank",
    scores = as.matrix(colMaps(ace)[["H_unified"]]),
    thread_no = thread_no,
    alpha = footprint_alpha
  )
  colMaps(ace)$archetype_footprint <- archetype_footprint

  # Compute gene specificity for each archetype
  ace <- .run.archetypeFeatureSpecificity(
    ace = ace,
    S = S,
    H = Matrix::t(archetype_footprint),
    assay_name = NULL,
    footprint_slot = NULL,
    return_raw = FALSE
  )

  # ace <- constructBackbone(
  #   ace = ace,
  #   network_density = network_density,
  #   mutual_edges_only = mutual_edges_only,
  #   layout_compactness = layout_compactness,
  #   layout_epochs = layout_epochs / 5,
  #   thread_no = 1,
  #   net_slot = net_slot_out
  # )

  if (full_trace == TRUE) {
    trace <- list(
      ACTION.out = ACTION.out,
      pruning.out = pruning.out,
      unification.out = unification.out
    )

    trace$log <- list(
      genes = rownames(ace),
      cells = colnames(ace),
      time = Sys.time()
    )

    out <- list(ace = ace, trace = trace)
    return(out)
  } else {
    return(ace)
  }
}


#' Reconstructs the ACTIONet graph with the new parameters (uses prior decomposition)
#'
#' @param ace ACTIONetExperiment object.
#' @param network_density Density factor of ACTIONet graph. (default=1)
#' @param distance_metric Distance metric with which to compute cell-to-cell similarity during network construction. Options are 'jsd' (Jensen-Shannon divergence), L2-norm ('l2'), and inner product ('ip'). (default='jsd')
#' @param nn_approach Nearest-neighbor algorithm to use for network construction. Options are k-nearest neighbors ('knn') and k*-nearest neighbors ('k*nn'). (default='k*nn')
#' @param mutual_edges_only Whether to enforce edges to be mutually-nearest-neighbors. (default=TRUE)
#' @param layout_compactness A value between 0-100, indicating the compactness of ACTIONet layout (default=50).
#' @param layout_epochs Number of epochs for SGD algorithm. (default=1000)
#' @param layout_algorithm Algorithm for computing plot layout. t-UMAP ("tumap") or UMAP ("umap"). Not case sensitive. (default="tumap")
#' @param layout_in_parallel Run layout construction using multiple cores. May result in marginally different outputs across runs due to parallelization-induced randomization. (default=TRUE)
#' @param thread_no Number of parallel threads. (default=0)
#' @param reduction_slot Slot in colMaps(ace) containing reduced kernel. (default='ACTION')
#' @param new_net_slot Name of slot in colMaps(ace) to store ACTIONet adjacency matrix. (default='ACTIONet')
#' @param seed Seed for random initialization. (default=0)
#'
#' @return ace Updated ace object
#'
#' @examples
#' plot.ACTIONet(ace)
#' ace.updated <- reconstructACTIONet(ace, network_density = 0.1)
#' plot.ACTIONet(ace.updated)
#' @export
reconstructACTIONet <- function(ace,
                                network_density = 1,
                                distance_metric = "jsd",
                                nn_approach = "k*nn",
                                mutual_edges_only = TRUE,
                                layout_compactness = 50,
                                layout_epochs = 1000,
                                layout_algorithm = c("tumap", "umap"),
                                layout_in_parallel = TRUE,
                                thread_no = 0,
                                reduction_slot = "ACTION",
                                new_net_slot = "ACTIONet",
                                H_stacked_slot = "H_stacked",
                                seed = 0) {
  set.seed(seed)
  layout_algorithm <- match.arg(toupper(layout_algorithm), choices = c("TUMAP", "UMAP"), several.ok = FALSE)

  if (!(reduction_slot %in% names(colMaps(ace)))) {
    err <- sprintf("Attribute '%s' is not in 'colMaps'.\n", reduction_slot)
    stop(err)
  }

  # re-Build ACTIONet
  if (!(H_stacked_slot %in% names(colMaps(ace)))) {
    err <- sprintf("Attribute '%s' is not in 'colMaps'.\n", H_stacked_slot)
    stop(err)
  }

  H_stacked <- Matrix::t(as.matrix(colMaps(ace)[[H_stacked_slot]]))

  G <- buildNetwork(
    H = H_stacked,
    algorithm = nn_approach,
    distance_metric = distance_metric,
    density = network_density,
    thread_no = thread_no,
    mutual_edges_only = mutual_edges_only
  )
  colNets(ace)[[new_net_slot]] <- G

  # Layout ACTIONet
  ace <- rerunLayout(
    ace = ace,
    layout_compactness = layout_compactness,
    layout_epochs = layout_epochs,
    layout_algorithm = layout_algorithm,
    thread_no = ifelse(layout_in_parallel, thread_no, 1),
    network_density = network_density,
    mutual_edges_only = mutual_edges_only,
    reduction_slot = reduction_slot,
    net_slot = new_net_slot,
    seed = seed
  )

  return(ace)
}


#' Rerun layout on the ACTIONet graph with new parameters
#'
#' @param ace ACTIONetExperiment object.
#' @param layout_compactness A value between 0-100, indicating the compactness of ACTIONet layout (default=50).
#' @param layout_epochs Number of epochs for SGD algorithm (default=1000).
#' @param layout_algorithm Algorithm for computing plot layout. t-UMAP ("tumap") or UMAP ("umap"). Not case sensitive. (default="tumap")
#' @param network_density Density factor of ACTIONet graph (default=1).
#' @param mutual_edges_only Whether to enforce edges to be mutually-nearest-neighbors (default=TRUE).
#' @param thread_no Number of parallel threads (default=0).
#' @param reduction_slot Slot in colMaps(ace) containing reduced kernel (default='ACTION').
#' @param net_slot Slot in colMaps(ace) containing ACTIONet adjacency matrix (default='ACTIONet').
#' @param seed Seed for random initialization (default=0).
#'
#' @return ace Updated ace object
#'
#' @examples
#' plot.ACTIONet(ace)
#' ace.updated <- rerunLayout(ace, layout_compactness = 20)
#' plot.ACTIONet(ace.updated)
#' @export
rerunLayout <- function(ace,
                        layout_compactness = 50,
                        layout_epochs = 1000,
                        layout_algorithm = c("tumap", "umap"),
                        network_density = 1,
                        mutual_edges_only = TRUE,
                        thread_no = 0,
                        reduction_slot = "ACTION",
                        net_slot = "ACTIONet",
                        seed = 0) {
  layout_alg <- match.arg(toupper(layout_algorithm), choices = c("TUMAP", "UMAP"), several.ok = FALSE)

  ace <- .run.layoutNetwork(
    ace = ace,
    G = NULL,
    initial_coordinates = NULL,
    compactness_level = layout_compactness,
    n_epochs = layout_epochs,
    layout_alg = layout_alg,
    thread_no = thread_no,
    reduction_slot = reduction_slot,
    net_slot = net_slot,
    seed = seed,
    return_raw = FALSE
  )

  # ace <- constructBackbone(
  #   ace = ace,
  #   network_density = network_density,
  #   mutual_edges_only = mutual_edges_only,
  #   layout_compactness = layout_compactness,
  #   layout_epochs = layout_epochs / 5,
  #   thread_no = 1,
  #   net_slot = net_slot_out
  # )

  return(ace)
}


#' @export
rerunArchAggr <- function(ace,
                          assay_name = "logcounts",
                          reduction_slot = "ACTION",
                          C_stacked_slot = "C_stacked",
                          H_stacked_slot = "H_stacked",
                          net_slot = "ACTIONet",
                          unified_suffix = "unified",
                          footprint_alpha = 0.85,
                          network_density = 1,
                          mutual_edges_only = TRUE,
                          layout_compactness = 50,
                          layout_epochs = 100,
                          thread_no = 0,
                          unification_violation_threshold = 0) {
  if (!(assay_name %in% names(assays(ace)))) {
    err <- sprintf("'%s' is not an assay of the input '%s' object.\n", assay_name, class(ace))
    stop(err)
  }

  if (!(reduction_slot %in% names(colMaps(ace)))) {
    err <- sprintf("Attribute '%s' is not in 'colMaps'.\n", reduction_slot)
    stop(err)
  }

  if (!(C_stacked_slot %in% names(colMaps(ace)))) {
    err <- sprintf("Attribute '%s' is not in 'colMaps'.\n", C_stacked_slot)
    stop(err)
  }

  if (!(H_stacked_slot %in% names(colMaps(ace)))) {
    err <- sprintf("Attribute '%s' is not in 'colMaps'.\n", H_stacked_slot)
    stop(err)
  }

  if (!(net_slot %in% names(colNets(ace)))) {
    err <- sprintf("Attribute '%s' is not in 'colNets'.\n", net_slot)
    stop(err)
  }

  S <- SummarizedExperiment::assays(ace)[[assay_name]]
  S_r <- Matrix::t(scale(colMaps(ace)[[reduction_slot]]))
  C_stacked <- as.matrix(colMaps(ace)[[C_stacked_slot]])
  H_stacked <- Matrix::t(as.matrix(colMaps(ace)[[H_stacked_slot]]))
  G <- colNets(ace)[[net_slot]]

  ace <- .run.unifyArchetypes(
    ace = NULL,
    S_r = S_r,
    C_stacked = C_stacked,
    H_stacked = H_stacked,
    reduction_slot = NULL,
    C_stacked_slot = NULL,
    H_stacked_slot = NULL,
    unified_suffix = unified_suffix,
    violation_threshold = unification_violation_threshold,
    thread_no = thread_no,
    return_raw = FALSE
  )

  # unification.out <- unify_archetypes(
  #   S_r = S_r,
  #   C_stacked = C_stacked,
  #   H_stacked = H_stacked,
  #   violation_threshold = unification_violation_threshold,
  #   thread_no = thread_no
  # )

  # colMaps(ace)[[sprintf("H_%s", unified_suffix)]] <- as(Matrix::t(unification.out$H_unified), "sparseMatrix")
  # colMapTypes(ace)[[sprintf("H_%s", unified_suffix)]] <- "internal"
  #
  # colMaps(ace)[[sprintf("C_%s", unified_suffix)]] <- as(unification.out$C_unified, "sparseMatrix")
  # colMapTypes(ace)[[sprintf("C_%s", unified_suffix)]] <- "internal"
  #
  # colData(ace)[["assigned_archetype"]] <- c(unification.out$assigned_archetype)

  # Use graph core of global and induced subgraphs to infer centrality/quality of
  # each cell
  ace$node_centrality <- c(compute_archetype_core_centrality(G, ace$assigned_archetype))

  # Ht_unified <- colMaps(ace)[["H_unified"]]
  archetype_footprint <- networkDiffusion(
    G = G,
    algorithm = "pagerank",
    scores = as.matrix(colMaps(ace)[[sprintf("H_%s", unified_suffix)]]),
    thread_no = thread_no,
    alpha = footprint_alpha
  )
  colMaps(ace)$archetype_footprint <- archetype_footprint

  # Ht_unified <- colMaps(ace)[[sprintf("H_%s", unified_suffix)]]
  # archetype_footprint <- compute_network_diffusion_fast(
  #   G = G,
  #   X0 = Ht_unified,
  #   thread_no = thread_no,
  #   alpha = footprint_alpha
  # )
  # colMaps(ace)$archetype_footprint <- archetype_footprint

  # Compute gene specificity for each archetype
  ace <- .run.archetypeFeatureSpecificity(
    ace = ace,
    S = S,
    H = Matrix::t(archetype_footprint),
    assay_name = NULL,
    footprint_slot = NULL,
    return_raw = FALSE
  )


  # H <- Matrix::t(archetype_footprint)
  #
  # # Compute gene specificity for each archetype
  # if (is.matrix(S)) {
  #   specificity.out <- compute_archetype_feature_specificity_full(S, H)
  # } else {
  #   specificity.out <- compute_archetype_feature_specificity(S, H)
  # }
  #
  # specificity.out <- lapply(specificity.out, function(specificity.scores) {
  #   rownames(specificity.scores) <- rownames(ace)
  #   colnames(specificity.scores) <- paste("A", 1:ncol(specificity.scores), sep = "")
  #   return(specificity.scores)
  # })
  #
  # rowMaps(ace)[[sprintf("%s_feature_profile", unified_suffix)]] <- specificity.out[["archetypes"]]
  # rowMapTypes(ace)[[sprintf("%s_feature_profile", unified_suffix)]] <- "internal"
  #
  # rowMaps(ace)[[sprintf("%s_feature_specificity", unified_suffix)]] <- specificity.out[["upper_significance"]]
  # rowMapTypes(ace)[[sprintf("%s_feature_specificity", unified_suffix)]] <- "reduction"

  # ace <- constructBackbone(
  #   ace = ace,
  #   network_density = network_density,
  #   mutual_edges_only = mutual_edges_only,
  #   layout_compactness = layout_compactness,
  #   layout_epochs = layout_epochs / 5,
  #   thread_no = 1,
  #   net_slot = net_slot_out
  # )

  return(ace)
}

constructBackbone <- function(ace,
                              network_density = 1,
                              mutual_edges_only = TRUE,
                              layout_algorithm = c("tumap", "umap"),
                              layout_compactness = 50,
                              layout_epochs = 100,
                              net_slot = "ACTIONet") {
  # if (!("archetype_footprint" %in% names(colMaps(ace)))) {
  #     G <- colNets(ace)[[net_slot]]
  #     Ht_unified <- colMaps(ace)[["H_unified"]]
  #
  #     archetype_footprint <- compute_network_diffusion_fast(
  #       G = G,
  #       X0 = Ht_unified,
  #       thread_no = thread_no,
  #       alpha = footprint_alpha
  #     )
  #     colMaps(ace)$archetype_footprint <- archetype_footprint
  #   }
  #
  #   W <- exp(scale(ace$archetype_footprint))
  #   W <- as(W, "sparseMatrix")
  #
  #   arch.vis.out <- transform_layout(
  #     W = W,
  #     coor2D = Matrix::t(ace$ACTIONet2D),
  #     coor3D = Matrix::t(ace$ACTIONet3D),
  #     colRGB = Matrix::t(ace$denovo_color),
  #     n_epochs = layout_epochs,
  #     compactness_level = layout_compactness,
  #     thread_no = thread_no
  #   )
  #
  #   arch.G <- computeFullSim(colMaps(ace)$archetype_footprint)
  #   diag(arch.G) <- 0
  #
  #   backbone <- list(
  #     G = arch.G,
  #     coordinates = Matrix::t(arch.vis.out$coordinates),
  #     coordinates_3D = Matrix::t(arch.vis.out$coordinates_3D),
  #     colors = Matrix::t(arch.vis.out$colors)
  #   )
  backbone <- list()
  metadata(ace)$backbone <- backbone

  return(ace)
}