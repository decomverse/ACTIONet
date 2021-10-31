run.ACTIONet.simple <- function(ace,
                         batch = NULL,
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
                         layout_algorithm = 0,
                         layout_in_parallel = TRUE,
                         unification_violation_threshold = 0,
                         footprint_alpha = 0.85,
                         thread_no = 0,
                         full_trace = FALSE,
                         seed = 0) {
  if (!(assay_name %in% names(assays(ace)))) {
    err <- sprintf("Attribute %s is not an assay of the input ace\n", assay_name)
    stop(err)
  }

  ace <- as(ace, "ACTIONetExperiment")

  S <- SummarizedExperiment::assays(ace)[[assay_name]]
  S_r <- Matrix::t(colMaps(ace)[[reduction_slot]])

  # Set parameters for ACTION_decomposition_MR
  params <- list()
  params$batch <- batch
  params$k_min <- k_min
  params$k_max <- k_max
  params$min_specificity_z_thresh <- min_specificity_z_thresh
  params$min_cells_per_arch <- min_cells_per_arch
  params$unification_violation_threshold <- unification_violation_threshold
  params$max_iter <- max_iter_ACTION
  params$thread_no <- thread_no
  params$seed <- seed

  # Run ACTION_decomposition_MR
  ACTION.out <- ACTIONet::decomp(X = S_r, params = params, method = "ACTION_decomposition_MR")

  # Store resulting decompositions
  colMaps(ace)[["H_stacked"]] <- Matrix::t(as(ACTION.out$extra$H_stacked, "sparseMatrix"))
  colMapTypes(ace)[["H_stacked"]] <- "internal"

  colMaps(ace)[["C_stacked"]] <- as(ACTION.out$extra$C_stacked, "sparseMatrix")
  colMapTypes(ace)[["C_stacked"]] <- "internal"

  H_unified <- as(Matrix::t(ACTION.out$extra$H_unified), "sparseMatrix")
  colMaps(ace)[["H_unified"]] <- H_unified
  colMapTypes(ace)[["H_unified"]] <- "internal"

  colMaps(ace)[["C_unified"]] <- as(ACTION.out$extra$C_unified, "sparseMatrix")
  colMapTypes(ace)[["C_unified"]] <- "internal"

  ace$assigned_archetype <- c(ACTION.out$extra$assigned_archetype)

  # Build ACTIONet
  set.seed(seed)
  H <- as.matrix(Matrix::t(colMaps(ace)[["H_stacked"]]))
  G <- build_ACTIONet(
    H_stacked = H,
    density = network_density,
    thread_no = thread_no,
    mutual_edges_only = mutual_edges_only,
    distance_metric = distance_metric,
    nn_approach = nn_approach
  )
  colNets(ace)[[net_slot_out]] <- G


  # Layout ACTIONet
  initial_coordinates <- ACTIONetExperiment:::.tscalet(S_r)
  colMaps(ace)[["ACTIONred"]] <- Matrix::t(initial_coordinates[1:3, ])
  colMapTypes(ace)[["ACTIONred"]] <- "embedding"

  ace <- .run.layout_ACTIONet(ace,
    G = G,
    S_r = initial_coordinates,
    compactness_level = layout_compactness,
    n_epochs = layout_epochs,
    layout_alg = layout_algorithm,
    thread_no = ifelse(layout_in_parallel, thread_no, 1),
    reduction_slot = NULL,
    net_slot = NULL,
    seed = seed
  )

  # Use graph core of global and induced subgraphs to infer centrality/quality of each cell
  ace$node_centrality <- c(compute_archetype_core_centrality(G, ace$assigned_archetype))

  # Smooth archetype footprints
  Ht_unified <- colMaps(ace)[["H_unified"]]
  archetype_footprint <- compute_network_diffusion_fast(
    G = G,
    X0 = Ht_unified,
    thread_no = thread_no,
    alpha = footprint_alpha
  )
  colMaps(ace)$archetype_footprint <- archetype_footprint


  # Compute gene specificity for each archetype
  H <- Matrix::t(archetype_footprint)
  if (is.matrix(S)) {
    specificity.out <- compute_archetype_feature_specificity_full(S, H)
  } else {
    specificity.out <- compute_archetype_feature_specificity(S, H)
  }

  specificity.out <- lapply(specificity.out, function(specificity.scores) {
    rownames(specificity.scores) <- rownames(ace)
    colnames(specificity.scores) <- paste("A", 1:ncol(specificity.scores), sep = "")
    return(specificity.scores)
  })

  rowMaps(ace)[["unified_feature_profile"]] <- specificity.out[["archetypes"]]
  rowMapTypes(ace)[["unified_feature_profile"]] <- "internal"

  rowMaps(ace)[["unified_feature_specificity"]] <- specificity.out[["upper_significance"]]
  rowMapTypes(ace)[["unified_feature_specificity"]] <- "reduction"

  ace <- construct.backbone(
    ace = ace,
    network_density = network_density,
    mutual_edges_only = mutual_edges_only,
    layout_compactness = layout_compactness,
    layout_epochs = layout_epochs / 5,
    thread_no = 1,
    ACTIONet_slot = net_slot_out
  )

  return(ace)
}


#' Filter columns and rows of `ACTIONetExperiment` or `SummarizedExperiment`-like object.
#' @export
filter.ace.test <- function(
  ace,
  assay_name = "counts",
  min_cells_per_feat = NULL,
  min_feats_per_cell = NULL,
  min_umis_per_cell = NULL,
  max_umis_per_cell = NULL,
  max_mito_fraction = NULL,
  species = c("mmusculus", "hsapiens"),
  features_use = NULL,
  return_fil_ace = TRUE
) {

    init_dim = dim(ace)
    init_dnames = dimnames(ace)

    X = SummarizedExperiment::assays(ace)[[assay_name]]
    X = as(X, "dgCMatrix")
    dimnames(X) = list(1:NROW(X), 1:NCOL(X))

    i = 0
    repeat {
        prev_dim = dim(X)
        rows_mask = rep(TRUE, NROW(ace))
        cols_mask = rep(TRUE, NCOL(ace))
        if (!is.null(min_umis_per_cell)) {
            umi_mask = Matrix::colSums(X) >= min_umis_per_cell
            cols_mask = cols_mask & umi_mask
        }

        if (!is.null(max_umis_per_cell)) {
            umi_mask = Matrix::colSums(X) <= max_umis_per_cell
            cols_mask = cols_mask & umi_mask
        }

        if (!is.null(min_feats_per_cell)) {
            feature_mask = Matrix::colSums(as(X > 0, "dgCMatrix")) >= min_feats_per_cell
            cols_mask = cols_mask & feature_mask
        }

        if (!is.null(min_cells_per_feat)) {
            if ((min_cells_per_feat < 1) & (min_cells_per_feat > 0)) {
                min_fc = min_cells_per_feat * prev_dim[2]
            } else {
                min_fc = min_cells_per_feat
            }
            cell_count_mask = ACTIONetExperiment:::fastRowSums(as(X > 0, "dgCMatrix")) >= min_fc
            rows_mask = rows_mask & cell_count_mask
        }

        X <- X[rows_mask, cols_mask]
        invisible(gc())
        i = i + 1
        if (all(dim(X) == prev_dim)) {
            break
        }
    }
    ace = ace[as.numeric(rownames(X)), as.numeric(colnames(X))]
    invisible(gc())

    if (!is.null(max_mito_fraction)){
      mt_frac = get_mtRNA_stats(
        ace,
        by = NULL,
        groups_use = NULL,
        assay = assay_name,
        species = species,
        metric = "pct",
        features_use = features_use
      )

      ace = ace[, mt_frac <= max_mito_fraction]
      invisible(gc())
    }

    if (return_fil_ace){
      return(ace)
    } else {
      fil_cols_mask = !(init_dnames[[2]] %in% colnames(ace))
      fil_rows_mask = !(init_dnames[[1]] %in% rownames(ace))

      fil_cols_list = data.frame(
        name = init_dnames[[2]][fil_cols_mask],
        idx = which(fil_cols_mask)
      )

      fil_rows_list = data.frame(
        name = init_dnames[[1]][fil_rows_mask],
        idx = which(fil_rows_mask)
      )

      fil_list = list(
        cols_filtered = fil_cols_list,
        rows_filtered = fil_rows_list
      )

      return(fil_list)
    }
}
