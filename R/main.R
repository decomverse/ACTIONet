#' A wrapper function to call all main functions of the ACTIONet
#'
#' @param ace Reduced `ACTIONetExperiment (ace)` object (output of reduce.ace() function).
#' @param k_max Maximum depth of decompositions (default=30).
#' @param min_specificity_z_threshold Defines the stringency of pruning nonspecific archetypes.
#' The larger the value, the more archetypes will be filtered out (default=-1).
#' @param network_density Density factor of ACTIONet graph (default=1).
#' @param mutual_edges_only Whether to enforce edges to be mutually-nearest-neighbors (default=TRUE).
#' @param layout_compactness A value between 0-100, indicating the compactness of ACTIONet layout (default=50).
#' @param layout_epochs Number of epochs for SGD algorithm (default=1000).
#' @param thread_no Number of parallel threads (default=0).
#' @param reduction_slot Slot in the colMaps(ace) that holds reduced kernel (default='S_r').
#' @param assay_name Name of assay to be used (default='logcounts').
#'
#' @return A named list: \itemize{
#' \item ace: ACTIONetExperiment object (derived from SummarizedExperiment)
#' \item ACTIONet.trace: Log of ACTIONet function calls
#'}
#'
#' @examples
#' ace = reduce(ace)
#' ACTIONet.out = run.ACTIONet(ace)
#' ace = ACTIONet.out$ace # main output
#' trace = ACTIONet.out$trace # for backup
#' @export
run.ACTIONet <- function(ace, k_max = 30, min_cells_per_arch = 2, min_specificity_z_threshold = -3, 
    network_density = 1, mutual_edges_only = TRUE, layout_compactness = 50, layout_epochs = 1000, 
    unification_alpha = 0.99, unification_outlier_threshold = 2, unification_sim_threshold = 0, 
    layout_in_parallel = TRUE, thread_no = 0, assay_name = "logcounts", reduction_slot = "ACTION", 
    footprint_alpha = 0.85, max_iter_ACTION = 50, full_trace = FALSE) {
    if (!(assay_name %in% names(assays(ace)))) {
        err = sprintf("Attribute %s is not an assay of the input ace\n", assay_name)
        stop(err)
    }
    
    ace = as(ace, "ACTIONetExperiment")
    
    
    S = assays(ace)[[assay_name]]
    S_r = Matrix::t(ACTIONet::colMaps(ace)[[reduction_slot]])
    
    # Run ACTION
    ACTION.out = run_ACTION(S_r, k_min = 2, k_max = k_max, thread_no = thread_no, 
        max_it = max_iter_ACTION, min_delta = 1e-300)
    
    # Prune nonspecific and/or unreliable archetypes
    pruning.out = prune_archetypes(ACTION.out$C, ACTION.out$H, min_specificity_z_threshold = min_specificity_z_threshold, 
        min_cells = min_cells_per_arch)
    
    colMaps(ace)[["H_stacked"]] = Matrix::t(as(pruning.out$H_stacked, "sparseMatrix"))
    colMapTypes(ace)[["H_stacked"]] = "internal"
    
    colMaps(ace)[["C_stacked"]] = as(pruning.out$C_stacked, "sparseMatrix")
    colMapTypes(ace)[["C_stacked"]] = "internal"
    
    
    # Build ACTIONet
    set.seed(0)
    G = build_ACTIONet(H_stacked = pruning.out$H_stacked, density = network_density, 
        thread_no = thread_no, mutual_edges_only = mutual_edges_only)
    colNets(ace)$ACTIONet = G
    
    
    # Layout ACTIONet
    initial.coordinates = t(scale(t(S_r)))
    colMaps(ace)[["ACTIONred"]] = Matrix::t(initial.coordinates[1:3, ])
    colMapTypes(ace)[["ACTIONred"]] = "embedding"
    
    if (layout_in_parallel == FALSE) {
        vis.out = layout_ACTIONet(G, S_r = initial.coordinates, compactness_level = layout_compactness, 
            n_epochs = layout_epochs, thread_no = 1)
    } else {
        # WARNING! This makes the results none reproducible
        vis.out = layout_ACTIONet(G, S_r = initial.coordinates, compactness_level = layout_compactness, 
            n_epochs = layout_epochs, thread_no = thread_no)
    }
    
    X = vis.out$coordinates
    colnames(X) = c("x", "y")
    rownames(X) = colnames(ace)
    colMaps(ace)$ACTIONet2D = X
    colMapTypes(ace)[["ACTIONet2D"]] = "embedding"
    
    X = vis.out$coordinates_3D
    colnames(X) = c("x", "y", "z")
    rownames(X) = colnames(ace)
    colMaps(ace)$ACTIONet3D = X
    colMapTypes(ace)[["ACTIONet3D"]] = "embedding"
    
    X = vis.out$colors
    colnames(X) = c("r", "g", "b")
    rownames(X) = colnames(ace)
    colMaps(ace)$denovo_color = X
    colMapTypes(ace)[["denovo_color"]] = "embedding"
    
    
    # Identiy equivalent classes of archetypes and group them together
    unification.out = unify_archetypes(G = G, S_r = S_r, C_stacked = pruning.out$C_stacked, 
        alpha = unification_alpha, outlier_threshold = unification_outlier_threshold, 
        sim_threshold = unification_sim_threshold, thread_no)
    metadata(ace)$selected_archetypes = unification.out$selected_archetypes
    metadata(ace)$selected_archetypes_ontology = unification.out$ontology
    metadata(ace)$selected_archetypes_ontology_annotations = unification.out$ontology_node_attributes
    
    Ht_unified = as(Matrix::t(unification.out$H_unified), "sparseMatrix")
    colMaps(ace)[["H_unified"]] = Ht_unified
    colMapTypes(ace)[["H_unified"]] = "internal"
    
    colMaps(ace)[["C_unified"]] = as(unification.out$C_unified, "sparseMatrix")
    colMapTypes(ace)[["C_unified"]] = "internal"
    
    ace$assigned_archetype = unification.out$assigned_archetype
    
    # Use graph core of global and induced subgraphs to infer centrality/quality of
    # each cell
    ace$node_centrality = compute_archetype_core_centrality(G, ace$assigned_archetype)
    
    # Smooth archetype footprints
    Ht_unified = colMaps(ace)[["H_unified"]]
    archetype_footprint = compute_network_diffusion(G, Ht_unified, alpha = footprint_alpha, 
        thread_no = thread_no)
    colMaps(ace)$archetype_footprint = archetype_footprint
    
    H = Matrix::t(archetype_footprint)
    
    # Compute gene specificity for each archetype
    if (is.matrix(S)) {
        specificity.out = compute_archetype_feature_specificity_full(S, H)
    } else {
        specificity.out = compute_archetype_feature_specificity(S, H)
    }
    
    specificity.out = lapply(specificity.out, function(specificity.scores) {
        rownames(specificity.scores) = rownames(ace)
        colnames(specificity.scores) = paste("A", 1:ncol(specificity.scores), sep = "")
        return(specificity.scores)
    })
    rowMaps(ace)[["unified_feature_profile"]] = specificity.out[["archetypes"]]
    rowMapTypes(ace)[["unified_feature_profile"]] = "internal"
    
    rowMaps(ace)[["unified_feature_specificity"]] = specificity.out[["upper_significance"]]
    rowMapTypes(ace)[["unified_feature_specificity"]] = "reduction"
    
    ace = construct.backbone(ace, network_density = network_density, mutual_edges_only = mutual_edges_only, 
        layout_compactness = layout_compactness, layout_epochs = layout_epochs/5, 
        thread_no = 1)
    
    if (full_trace == T) {
        # Prepare output
        trace = list(ACTION.out = ACTION.out, pruning.out = pruning.out, vis.out = vis.out, 
            unification.out = unification.out)
        trace$log = list(genes = rownames(ace), cells = colnames(ace), time = Sys.time())
        
        metadata(ace)$trace = trace
    }
    
    return(ace)
}


#' Reconstructs the ACTIONet graph with the new parameters (uses prior decomposition)
#'
#' @param ace ACTIONetExperiment object containing the results
#' @param network_density Density factor of ACTIONet graph (default=1).
#' @param mutual_edges_only Whether to enforce edges to be mutually-nearest-neighbors (default=TRUE).
#' @param compactness_level A value between 0-100, indicating the compactness of ACTIONet layout (default=50)
#' @param n_epochs Number of epochs for SGD algorithm (default=500).
#' @param thread_no Number of parallel threads (default=0)
#' @param reduction_slot Slot in the colMaps(ace) that holds reduced kernel (default='S_r')
#'
#' @return ace Updated ace object
#'
#' @examples
#' plot.ACTIONet(ace)
#' ace.updated = reconstruct.ACTIONet(ace, network_density = 0.1)
#' plot.ACTIONet(ace.updated)
#' @export
reconstruct.ACTIONet <- function(ace, network_density = 1, mutual_edges_only = TRUE, 
    layout_compactness = 50, layout_epochs = 1000, thread_no = 0, layout_in_parallel = FALSE, 
    reduction_slot = "ACTION") {
    set.seed(0)
    
    # re-Build ACTIONet
    H_stacked = Matrix::t(as.matrix(colMaps(ace)[["H_stacked"]]))
    
    G = build_ACTIONet(H_stacked = H_stacked, density = network_density, thread_no = thread_no, 
        mutual_edges_only = mutual_edges_only)
    colNets(ace)$ACTIONet = G
    
    
    # Layout ACTIONet
    initial.coordinates = Matrix::t(scale(ACTIONet::colMaps(ace)[[reduction_slot]]))
    if (layout_in_parallel == FALSE) {
        vis.out = layout_ACTIONet(G, S_r = initial.coordinates, compactness_level = layout_compactness, 
            n_epochs = layout_epochs, thread_no = 1)
    } else {
        # WARNING! This makes the results none reproducible
        vis.out = layout_ACTIONet(G, S_r = initial.coordinates, compactness_level = layout_compactness, 
            n_epochs = layout_epochs, thread_no = thread_no)
    }
    
    X = vis.out$coordinates
    colnames(X) = c("x", "y")
    rownames(X) = colnames(ace)
    colMaps(ace)$ACTIONet2D = X
    colMapTypes(ace)[["ACTIONet2D"]] = "embedding"
    
    X = vis.out$coordinates_3D
    colnames(X) = c("x", "y", "z")
    rownames(X) = colnames(ace)
    colMaps(ace)$ACTIONet3D = X
    colMapTypes(ace)[["ACTIONet3D"]] = "embedding"
    
    X = vis.out$colors
    colnames(X) = c("r", "g", "b")
    rownames(X) = colnames(ace)
    colMaps(ace)$denovo_color = X
    colMapTypes(ace)[["denovo_color"]] = "embedding"
    
    ace = construct.backbone(ace, network_density = network_density, mutual_edges_only = mutual_edges_only, 
        layout_compactness = layout_compactness, layout_epochs = layout_epochs/5, 
        thread_no = 1)
    
    return(ace)
}



#' Rerun layout on the ACTIONet graph with new parameters
#'
#' @param ace ACTIONetExperiment object containing the results
#' @param layout_compactness A value between 0-100, indicating the compactness of ACTIONet layout (default=50)
#' @param layout_epochs Number of epochs for SGD algorithm (default=500).
#' @param thread_no Number of parallel threads (default=8)
#' @param reduction_slot Slot in the colMaps(ace) that holds reduced kernel (default='S_r')
#'
#' @return ace Updated ace object
#'
#' @examples
#' plot.ACTIONet(ace)
#' ace.updated = rerun.layout(ace, layout_compactness = 20)
#' plot.ACTIONet(ace.updated)
#' @export
rerun.layout <- function(ace, layout_compactness = 50, layout_epochs = 1000, thread_no = 0, 
    network_density = 1, mutual_edges_only = T, reduction_slot = "ACTIONet3D", net_slot = "ACTIONet") {
    G = colNets(ace)[[net_slot]]
    
    # re-Layout ACTIONet
    S_r = Matrix::t(ACTIONet::colMaps(ace)[[reduction_slot]])
    
    initial.coordinates = t(scale(t(S_r)))
    vis.out = layout_ACTIONet(G, S_r = initial.coordinates, compactness_level = layout_compactness, 
        n_epochs = layout_epochs, thread_no = thread_no)
    
    X = vis.out$coordinates
    colnames(X) = c("x", "y")
    rownames(X) = colnames(ace)
    colMaps(ace)$ACTIONet2D = X
    colMapTypes(ace)[["ACTIONet2D"]] = "embedding"
    
    X = vis.out$coordinates_3D
    colnames(X) = c("x", "y", "z")
    rownames(X) = colnames(ace)
    colMaps(ace)$ACTIONet3D = X
    colMapTypes(ace)[["ACTIONet3D"]] = "embedding"
    
    X = vis.out$colors
    colnames(X) = c("r", "g", "b")
    rownames(X) = colnames(ace)
    colMaps(ace)$denovo_color = X
    colMapTypes(ace)[["denovo_color"]] = "embedding"
    
    ace = construct.backbone(ace, network_density = network_density, mutual_edges_only = mutual_edges_only, 
        layout_compactness = layout_compactness, layout_epochs = layout_epochs/5, 
        thread_no = 1)
    
    return(ace)
}

#' @export
rerun.archetype.aggregation <- function(ace, assay_name = "logcounts", reduction_slot = "ACTION", 
    unified_suffix = "unified", footprint_alpha = 0.85, network_density = 1, mutual_edges_only = TRUE, 
    layout_compactness = 50, layout_epochs = 100, thread_no = 0, unification_alpha = 0.99, 
    unification_outlier_threshold = 2, unification_sim_threshold = 0) {
    
    S = assays(ace)[[assay_name]]
    S_r = Matrix::t(ACTIONet::colMaps(ace)[[reduction_slot]])
    C_stacked = as.matrix(colMaps(ace)[["C_stacked"]])
    H_stacked = Matrix::t(as.matrix(colMaps(ace)[["H_stacked"]]))
    G = colNets(ace)[["ACTIONet"]]
    
    unification.out = unify_archetypes(G = G, S_r = S_r, C_stacked = C_stacked, alpha = unification_alpha, 
        outlier_threshold = unification_outlier_threshold, sim_threshold = unification_sim_threshold, 
        thread_no)
    metadata(ace)$selected_archetypes = unification.out$selected_archetypes
    metadata(ace)$selected_archetypes_ontology = unification.out$ontology
    metadata(ace)$selected_archetypes_ontology_annotations = unification.out$ontology_node_attributes
    
    
    colMaps(ace)[[sprintf("H_%s", unified_suffix)]] = as(Matrix::t(unification.out$H_unified), 
        "sparseMatrix")
    colMapTypes(ace)[[sprintf("H_%s", unified_suffix)]] = "internal"
    
    colMaps(ace)[[sprintf("C_%s", unified_suffix)]] = as(unification.out$C_unified, 
        "sparseMatrix")
    colMapTypes(ace)[[sprintf("C_%s", unified_suffix)]] = "internal"
    
    colData(ace)[["assigned_archetype"]] = unification.out$assigned_archetype
    
    # Use graph core of global and induced subgraphs to infer centrality/quality of
    # each cell
    ace$node_centrality = compute_archetype_core_centrality(G, ace$assigned_archetype)
    
    Ht_unified = colMaps(ace)[[sprintf("H_%s", unified_suffix)]]
    archetype_footprint = compute_network_diffusion(G, Ht_unified, alpha = footprint_alpha, 
        thread_no = thread_no)
    colMaps(ace)$archetype_footprint = archetype_footprint
    
    
    H = Matrix::t(archetype_footprint)
    
    # Compute gene specificity for each archetype
    if (is.matrix(S)) {
        specificity.out = compute_archetype_feature_specificity_full(S, H)
    } else {
        specificity.out = compute_archetype_feature_specificity(S, H)
    }
    
    specificity.out = lapply(specificity.out, function(specificity.scores) {
        rownames(specificity.scores) = rownames(ace)
        colnames(specificity.scores) = paste("A", 1:ncol(specificity.scores), sep = "")
        return(specificity.scores)
    })
    
    rowMaps(ace)[[sprintf("%s_feature_profile", unified_suffix)]] = specificity.out[["archetypes"]]
    rowMapTypes(ace)[[sprintf("%s_feature_profile", unified_suffix)]] = "internal"
    
    rowMaps(ace)[[sprintf("%s_feature_specificity", unified_suffix)]] = specificity.out[["upper_significance"]]
    rowMapTypes(ace)[[sprintf("%s_feature_specificity", unified_suffix)]] = "reduction"
    
    ace = construct.backbone(ace, network_density = network_density, mutual_edges_only = mutual_edges_only, 
        layout_compactness = layout_compactness, layout_epochs = layout_epochs/5, 
        thread_no = 1)
    
    return(ace)
}

construct.backbone <- function(ace, network_density = 1, mutual_edges_only = TRUE, 
    layout_compactness = 50, layout_epochs = 100, thread_no = 1, footprint_alpha = 0.85, 
    ACTIONet_slot = "ACTIONet") {
    if (!("archetype_footprint" %in% names(colMaps(ace)))) {
        G = colNets(ace)[[ACTIONet_slot]]
        Ht_unified = colMaps(ace)[["H_unified"]]
        archetype_footprint = compute_network_diffusion(G, Ht_unified, alpha = footprint_alpha, 
            thread_no = thread_no)
        colMaps(ace)$archetype_footprint = archetype_footprint
    }
    
    
    W = exp(scale(ace$archetype_footprint))
    W = as(W, "sparseMatrix")
    arch.vis.out = transform_layout(W, coor2D = Matrix::t(ace$ACTIONet2D), coor3D = Matrix::t(ace$ACTIONet3D), 
        colRGB = Matrix::t(ace$denovo_color), n_epochs = layout_epochs, compactness_level = layout_compactness, 
        thread_no = thread_no)
    # arch.G = build_ACTIONet(H_stacked = colMaps(ace)$archetype_footprint, density =
    # network_density, thread_no = thread_no, mutual_edges_only = mutual_edges_only)
    arch.G = ACTIONet::computeFullSim(colMaps(ace)$archetype_footprint)
    diag(arch.G) = 0
    
    backbone = list(G = arch.G, coordinates = Matrix::t(arch.vis.out$coordinates), 
        coordinates_3D = Matrix::t(arch.vis.out$coordinates_3D), colors = Matrix::t(arch.vis.out$colors))
    
    metadata(ace)$backbone = backbone
    
    return(ace)
}
