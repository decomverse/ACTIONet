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
#' @param specificity_th Defines the stringency of pruning nonspecific archetypes.
#' The larger the value, the more archetypes will be filtered out. (default=-3)
#' @param network_metric Distance metric with which to compute cell-to-cell similarity during network construction. Options are 'jsd' (Jensen-Shannon divergence), L2-norm ('l2'), and inner product ('ip'). (default='jsd')
#' @param network_algorithm Algorithm to use for network construction. Options are k-nearest neighbors ('knn') and k*-nearest neighbors ('k*nn'). (default='k*nn')
#' @param network_density Density factor of ACTIONet graph. (default=1)
#' @param mutual_edges_only Whether to enforce edges to be mutually-nearest-neighbors. (default=TRUE)
#' @param imputation_alpha Network diffusion parameter to smooth PCs (S_r). (default=0.9)
#' @param layout_compactness A value between 0-100, indicating the compactness of ACTIONet layout. (default=50)
#' @param layout_epochs Number of epochs for SGD algorithm. (default=250)
#' @param layout_algorithm Algorithm for computing plot layout. Options are UMAP ("umap") or t-UMAP ("tumap"). (default="umap")
#' @param layout_parallel Run layout construction using multiple cores. May result in marginally different outputs across runs due to parallelization-induced randomization. (default=TRUE)
#' @param unification_th Archetype unification resolution parameter. (default=0)
#' @param footprint_alpha Archetype smoothing parameter. (default=0.85)
#' @param compute_specificity_parallel Run feature specificity enrichment using multiple cores. Setting this to `TRUE` on large datasets may cause an out of memory crash. (default=FALSE)
#' @param smoothPC Whether to smooth output of 'reduce.ace()'. Performs diffusion smoothing on the components of the reduction. Skipped if 'TRUE' and output of 'reduce.ace()' is missing from object. (default=TRUE)
#' @param thread_no Number of parallel threads. (default=0)
#' @param seed Seed for random initialization. (default=0)
#'
#' @return An ACTIONetExperiment object
#'
#' @examples
#' ace <- runACTIONet(ace)
#' @export
runACTIONet <- function(ace,
                        k_min = 2,
                        k_max = 30,
                        assay_name = NULL,
                        reduction_slot = NULL,
                        reduction_normalization = 1,
                        net_slot_out = "ACTIONet",
                        min_cells_per_arch = 2,
                        max_iter_ACTION = 50,
                        specificity_th = -3,
                        network_metric = "jsd",
                        network_algorithm = "k*nn",
                        network_density = 1,
                        mutual_edges_only = TRUE,
                        imputation_alpha = 0.9,
                        layout_spread = 1.0,
                        layout_min_dist = 1.0,
                        layout_gamma = 1.0,
                        layout_epochs = 250,
                        layout_algorithm = c("umap", "tumap"),
                        layout_parallel = TRUE,
                        unification_backbone_density = 0.5,
                        unification_resolution = 1.0,
                        unification_min_cluster_size = 3,
                        footprint_alpha = 0.85,
                        compute_specificity_parallel = FALSE,
                        smoothPC = TRUE,
                        thread_no = 0,
                        backbone.density = 0.5,
                        seed = 0) {
  if (!is(ace, "ACTIONetExperiment")) {
    ace <- as(ace, "ACTIONetExperiment")
  }

  layout_algorithm <- tolower(layout_algorithm)
  layout_algorithm <- match.arg(layout_algorithm, several.ok = FALSE)

  set.seed(seed)

  if (is.null(assay_name)) {
    if ("default_assay" %in% names(metadata(ace))) {
      assay_name <- metadata(ace)[["default_assay"]]
      message(sprintf("Input assay_name is NULL. Setting assay_name to the metadata(ace)[['default_assay']] => %s", assay_name))
    } else {
      message(sprintf("Input assay_name is NULL. Setting assay_name to logcounts"))
      assay_name <- "logcounts"
    }
  }
  .validate_assay(ace, assay_name = assay_name, return_elem = FALSE)

  if (is.null(reduction_slot)) {
    if ("default_reduction" %in% names(metadata(ace))) {
      reduction_slot <- metadata(ace)[["default_reduction"]]
      message(sprintf("Input reduction_slot is NULL. Setting reduction_slot to the metadata(ace)[['default_reduction']] => %s", reduction_slot))
    } else {
      message(sprintf("Input reduction_slot is NULL. Setting reduction_slot to ACTION"))
      reduction_slot <- "ACTION"
    }
  }
  .validate_map(ace = ace, map_slot = reduction_slot, return_elem = FALSE)

  if (!(sprintf("%s_normalized", reduction_slot) %in% names(colMaps(ace)))) {
    ace <- normalize.reduction(ace, reduction_slot = reduction_slot, reduction_normalization = reduction_normalization)
  }
  reduction_slot <- sprintf("%s_normalized", reduction_slot)

  # Calls "decomp.ACTIONMR" and fills in the appropriate slots in the ace object
  ace <- .run.ACTIONMR.ace(
    ace = ace,
    k_min = k_min,
    k_max = k_max,
    specificity_th = specificity_th,
    min_cells_per_arch = min_cells_per_arch,
    unification_backbone_density = unification_backbone_density,
    unification_resolution = unification_resolution,
    unification_min_cluster_size = unification_min_cluster_size,
    max_iter = max_iter_ACTION,
    thread_no = thread_no,
    unified_suffix = "unified",
    footprint_slot_name = "assigned_archetype",
    reduction_slot = reduction_slot,
  )

  # Build ACTIONet
  G <- buildNetwork(
    H = as.matrix(Matrix::t(colMaps(ace)[["H_stacked"]])),
    algorithm = network_algorithm,
    distance_metric = network_metric,
    density = network_density,
    thread_no = thread_no,
    mutual_edges_only = mutual_edges_only
  )
  colNets(ace)[[net_slot_out]] <- G

  colData(ace)[["node_centrality"]] <- networkCentrality(
    obj = G,
    label_attr = colData(ace)[["assigned_archetype"]],
    algorithm = "local_coreness"
  )

  # Uses "layoutNetwork" to layout the network and fill in the appropriate slots in the ace object
  ace <- .run.layoutNetwork(
    ace = ace,
    algorithm = layout_algorithm,
    n_epochs = layout_epochs,
    spread = layout_spread,
    min_dist = layout_min_dist,
    gamma = layout_gamma,
    initial_coordinates = reduction_slot,
    net_slot = net_slot_out,
    thread_no = ifelse(layout_parallel, thread_no, 1),
    seed = seed
  )


  # # Smooth PCs (S_r) for ease of future imputation (same as MAGIC algorithm)
  # if (smoothPC == TRUE) {
  #   ace <- .smoothPCs(
  #     ace = ace,
  #     diffusion_algorithm = "pagerank",
  #     alpha = imputation_alpha,
  #     diffusion_it = 5,
  #     reduction_slot = reduction_slot,
  #     net_slot = net_slot_out,
  #     thread_no = thread_no
  #   )
  # }


  # Smooth archetype footprints
  archetype_footprint <- networkDiffusion(
    obj = G,
    scores = colMaps(ace)[["H_unified"]],
    algorithm = "pagerank",
    alpha = footprint_alpha,
    thread_no = thread_no,
    max_it = 5,
    res_threshold = 1e-8,
    net_slot = NULL
  )
  colMaps(ace)$archetype_footprint <- archetype_footprint

  # Compute gene specificity for each archetype
  ace <- archetypeFeatureSpecificity(
    ace = ace,
    assay_name = assay_name,
    footprint_slot = "archetype_footprint",
    thread_no = ifelse(compute_specificity_parallel, thread_no, 1),
    return_raw = FALSE
  )

  ace <- constructBackbone(ace, backbone.density = backbone.density)
  return(ace)
}

#' Reconstructs the ACTIONet graph with the new parameters (uses prior decomposition)
#'
#' @param ace ACTIONetExperiment object.
#' @param network_density Density factor of ACTIONet graph. (default=1)
#' @param network_metric Distance metric with which to compute cell-to-cell similarity during network construction. Options are 'jsd' (Jensen-Shannon divergence), L2-norm ('l2'), and inner product ('ip'). (default='jsd')
#' @param algorithm Algorithm to use for network construction. Options are k-nearest neighbors ('knn') and k*-nearest neighbors ('k*nn'). (default='k*nn')
#' @param mutual_edges_only Whether to enforce edges to be mutually-nearest-neighbors. (default=TRUE)
#' @param layout_compactness A value between 0-100, indicating the compactness of ACTIONet layout (default=50).
#' @param layout_epochs Number of epochs for SGD algorithm. (default=250)
#' @param layout_algorithm Algorithm for computing plot layout. Options are UMAP ("umap") or t-UMAP ("tumap"). (default="umap")
#' @param layout_parallel Run layout construction using multiple cores. May result in marginally different outputs across runs due to parallelization-induced randomization. (default=TRUE)
#' @param thread_no Number of parallel threads. (default=0)
#' @param reduction_slot Slot in colMaps(ace) containing reduced kernel. (default='ACTION')
#' @param net_slot_out Name of slot in colMaps(ace) to store ACTIONet adjacency matrix. (default='ACTIONet')
#' @param seed Seed for random initialization. (default=0)
#'
#' @return ace Updated ace object
#'
#' @examples
#' plot.ACTIONet(ace)
#' ace.updated <- rebuildACTIONet(ace, network_density = 0.1)
#' plot.ACTIONet(ace.updated)
#' @export
rebuildACTIONet <- function(ace,
                            network_density = 1,
                            network_metric = "jsd",
                            algorithm = "k*nn",
                            mutual_edges_only = TRUE,
                            layout_spread = 1.0,
                            layout_min_dist = 1.0,
                            layout_gamma = 1.0,
                            layout_epochs = 250,
                            layout_algorithm = c("umap", "tumap"),
                            layout_parallel = TRUE,
                            thread_no = 0,
                            reduction_slot = "ACTION",
                            net_slot_out = "ACTIONet",
                            H_stacked_slot = "H_stacked",
                            reduction_normalization = 1,
                            seed = 0) {
  set.seed(seed)
  layout_algorithm <- tolower(layout_algorithm)
  layout_algorithm <- match.arg(layout_algorithm, several.ok = FALSE)

  .validate_ace(ace = ace, return_elem = FALSE)
  .validate_map(ace = ace, map_slot = reduction_slot, return_elem = FALSE)

  if (!(sprintf("%s_normalized", reduction_slot) %in% names(colMaps(ace)))) {
    ace <- normalize.reduction(ace, reduction_slot = reduction_slot, reduction_normalization = reduction_normalization)
  }
  reduction_slot <- sprintf("%s_normalized", reduction_slot)

  # re-Build ACTIONet
  H_stacked <- .validate_map(
    ace = ace,
    map_slot = H_stacked_slot,
    matrix_type = "dense",
    force_type = TRUE
  )
  H_stacked <- Matrix::t(H_stacked)

  G <- buildNetwork(
    H = H_stacked,
    algorithm = algorithm,
    distance_metric = network_metric,
    density = network_density,
    thread_no = thread_no,
    mutual_edges_only = mutual_edges_only
  )
  colNets(ace)[[net_slot_out]] <- G

  # Layout ACTIONet
  ace <- rerunLayout(
    ace = ace,
    algorithm = layout_algorithm,
    n_epochs = layout_epochs,
    spread = layout_spread,
    min_dist = layout_min_dist,
    gamma = layout_gamma,
    reduction_slot = reduction_slot,
    net_slot = net_slot_out,
    thread_no = ifelse(layout_parallel, thread_no, 1),
    seed = seed
  )

  return(ace)
}


#' Rerun layout on the ACTIONet graph with new parameters
#'
#' @param ace ACTIONetExperiment object.
#' @param layout_compactness A value between 0-100, indicating the compactness of ACTIONet layout (default=50).
#' @param layout_epochs Number of epochs for SGD algorithm (default=250).
#' @param layout_algorithm Algorithm for computing plot layout. Options are UMAP ("umap") or t-UMAP ("tumap"). (default="umap")
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
#' ace.updated <- rerunLayout(ace, compactness = 20)
#' plot.ACTIONet(ace.updated)
#' @export
rerunLayout <- function(ace,
                        algorithm = c("umap", "tumap"),
                        spread = 1.0,
                        min_dist = 1.0,
                        gamma = 1.0,
                        n_epochs = 250,
                        thread_no = 0,
                        reduction_slot = "ACTION",
                        net_slot = "ACTIONet",
                        reduction_normalization = 1,
                        seed = 0) {
  algorithm <- tolower(algorithm)
  algorithm <- match.arg(algorithm, several.ok = FALSE)

  if (!(sprintf("%s_normalized", reduction_slot) %in% names(colMaps(ace)))) {
    ace <- normalize.reduction(ace, reduction_slot = reduction_slot, reduction_normalization = reduction_normalization)
  }
  reduction_slot <- sprintf("%s_normalized", reduction_slot)

  ace <- .run.layoutNetwork(
    ace = ace,
    algorithm = algorithm,
    n_epochs = n_epochs,
    spread = spread,
    min_dist = min_dist,
    gamma = gamma,
    initial_coordinates = reduction_slot,
    net_slot = net_slot,
    thread_no = thread_no,
    seed = seed
  )

  return(ace)
}


#' @export
rerun.archetype.unification <- function(ace,
                                        backbone_density = 0.5,
                                        resolution = 1.0,
                                        min_cluster_size = 3,
                                        footprint_alpha = 0.85,
                                        compute_specificity_parallel = FALSE,
                                        assay_name = "logcounts",
                                        reduction_slot = "ACTION",
                                        C_stacked_slot = "C_stacked",
                                        H_stacked_slot = "H_stacked",
                                        net_slot = "ACTIONet",
                                        unified_suffix = "unified",
                                        backbone.density = 0.5,
                                        reduction_normalization = 1,
                                        thread_no = 0) {
  .validate_ace(ace = ace, return_elem = FALSE)

  G <- .validate_net(
    ace = ace,
    net_slot = net_slot,
    matrix_type = "sparse",
    force_type = FALSE
  )

  if (!(sprintf("%s_normalized", reduction_slot) %in% names(colMaps(ace)))) {
    ace <- normalize.reduction(ace, reduction_slot = reduction_slot, reduction_normalization = reduction_normalization)
  }
  reduction_slot <- sprintf("%s_normalized", reduction_slot)

  ace <- .run.unifyArchetypes(
    ace = ace,
    unification_backbone_density = backbone_density,
    unification_resolution = resolution,
    unification_min_cluster_size = min_cluster_size,
    unification_normalization = 0,
    reduction_slot = reduction_slot,
    C_stacked_slot = C_stacked_slot,
    H_stacked_slot = H_stacked_slot,
    unified_suffix = unified_suffix,
    thread_no = thread_no
  )

  colData(ace)[["node_centrality"]] <- networkCentrality(
    obj = G,
    label_attr = colData(ace)[["assigned_archetype"]],
    algorithm = "local_coreness"
  )

  archetype_footprint <- networkDiffusion(
    obj = G,
    scores = colMaps(ace)[[sprintf("H_%s", unified_suffix)]],
    algorithm = "pagerank",
    alpha = footprint_alpha,
    thread_no = thread_no,
    max_it = 5,
    res_threshold = 1e-8,
    net_slot = NULL
  )
  colMaps(ace)$archetype_footprint <- archetype_footprint

  # Compute gene specificity for each archetype
  ace <- archetypeFeatureSpecificity(
    ace = ace,
    assay_name = assay_name,
    footprint_slot = "archetype_footprint",
    thread_no = ifelse(compute_specificity_parallel, thread_no, 1),
    return_raw = FALSE
  )

  ace <- constructBackbone(ace, backbone.density = backbone_density)

  return(ace)
}

constructBackbone <- function(ace,
                              backbone.density = 0.5) {
  footprint <- ace$archetype_footprint
  footprint <- normalize_mat(footprint, 1)
  arch.graph <- buildNetwork(footprint, density = backbone.density)
  arch.coors <- as.matrix(Matrix::t(colMaps(ace)$C_unified) %*% ace$ACTIONet2D)
  arch.coors_3D <- as.matrix(Matrix::t(colMaps(ace)$C_unified) %*% ace$ACTIONet3D)
  arch.colors <- as.matrix(Matrix::t(colMaps(ace)$C_unified) %*% ace$denovo_color)
  arch.colors[arch.colors > 1] <- 1
  backbone <- list(graph = arch.graph, coordinates = arch.coors, coordinates_3D = arch.coors_3D, colors = arch.colors)

  metadata(ace)$backbone <- backbone

  return(ace)
}


normalize.reduction <- function(ace, reduction_slot = "ACTION", reduction_normalization = 1) {
  X <- Matrix::t(colMaps(ace)[[sprintf("%s", reduction_slot)]])

  X_norm <- ACTIONet::normalize_mat(X, reduction_normalization)

  colMaps(ace)[[sprintf("%s_normalized", reduction_slot)]] <- Matrix::t(X_norm)
  metadata(ace)[["reduction_normalization"]] <- reduction_normalization

  return(ace)
}