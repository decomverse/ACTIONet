#' A wrapper function to call all main functions of the ACTIONet
#'
#' @param ace Reduced `ACTIONetExperiment (ACE)` object (output of reduce.ace() function).
#' @param k_max Maximum depth of decompositions (default=30).
#' @param min_specificity_z_threshold Defines the stringency of pruning nonspecific archetypes. 
#' The larger the value, the more archetypes will be filtered out (default=-1).
#' @param network_density Density factor of ACTIONet graph (default=1).
#' @param mutual_edges_only Whether to enforce edges to be mutually-nearest-neighbors (default=TRUE).
#' @param compactness_level A value between 0-100, indicating the compactness of ACTIONet layout (default=50)
#' @param n_epochs Number of epochs for SGD algorithm (default=500).
#' @param thread_no Number of parallel threads (default=0)
#' @param reduction.slot Slot in the colFactors(ace) that holds reduced kernel (default="S_r")
#' @param data.slot Corresponding slot in the `ace` object the normalized counts (default="logcounts")
#' @param renormalize.logcounts.slot Name of the new assay with updated logcounts adjusted using archetypes
#' If it is NULL, values of logcounts(ace) would be directly used without renormalization for computing speicificity scores
#' 
#' @return A named list: \itemize{
#' \item ace: ACTIONetExperiment object (derived from RangedSummarizedExperiment)
#' \item ACTIONet.trace: Log of ACTIONet function calls 
#'}
#' 
#' @examples
#' ace = reduce(ace)
#' ACTIONet.out = run.ACTIONet(ace)
#' ace = ACTIONet.out$ace # main output
#' trace = ACTIONet.out$trace # for backup
run.ACTIONet <- function(ace, k_max = 30, min.cells.per.arch = 2, min_specificity_z_threshold = 0, network_density = 1, mutual_edges_only = TRUE, layout_compactness = 50, layout_epochs = 500, layout.in.parallel = FALSE, thread_no = 0, data.slot = "logcounts", reduction.slot = "ACTION", unification.resolution = 1, AA_delta = 1e-300, full.trace = T) {		
    if (!(data.slot %in% names(assays(ace)))) {
        R.utils::printf("Attribute %s is not an assay of the input ace\n", data.slot)
        return()
    }
    
    ace = as(ace, "ACTIONetExperiment") 

    
	S = assays(ace)[[data.slot]]
    S_r = t(ACTIONet::colFactors(ace)[[reduction.slot]])
    
    # Run ACTION
	ACTION.out = run_ACTION(S_r, k_min = 2, k_max = k_max, thread_no = thread_no, max_it = 50, min_delta = AA_delta)
    
    # Prune nonspecific and/or unreliable archetypes
    pruning.out = prune_archetypes(ACTION.out$C, ACTION.out$H, min_specificity_z_threshold = min_specificity_z_threshold, min_cells = min.cells.per.arch)

	colFactors(ace)[["H_stacked"]] = as(pruning.out$H_stacked, 'sparseMatrix')
	colFactors(ace)[["C_stacked"]] = as(Matrix::t(pruning.out$C_stacked), 'sparseMatrix')

    
    # Build ACTIONet
    set.seed(0)
    G = build_ACTIONet(H_stacked = pruning.out$H_stacked, density = network_density, thread_no=thread_no, mutual_edges_only = mutual_edges_only)
	colNets(ace)$ACTIONet = G
	
	
    # Layout ACTIONet
	initial.coordinates = t(scale(t(S_r)))
	if(layout.in.parallel == FALSE) {		
		vis.out = layout_ACTIONet(G, S_r = initial.coordinates, compactness_level = layout_compactness, n_epochs = layout_epochs, thread_no = 1)
    } else { # WARNING! This makes the results none reproducible
		vis.out = layout_ACTIONet(G, S_r = initial.coordinates, compactness_level = layout_compactness, n_epochs = layout_epochs, thread_no = thread_no)
	}
	
    colFactors(ace)$ACTIONet2D = vis.out$coordinates
    colFactors(ace)$ACTIONet3D = vis.out$coordinates_3D
    colFactors(ace)$denovo_color = vis.out$colors


	# Identiy equivalent classes of archetypes and group them together
	unification.out = unify_archetypes(S_r, pruning.out$C_stacked, pruning.out$H_stacked, min_overlap = 0, resolution = unification.resolution)

	colFactors(ace)[["H_unified"]] = as(unification.out$H_unified, 'sparseMatrix')
	colFactors(ace)[["C_unified"]] = as(Matrix::t(unification.out$C_unified), 'sparseMatrix');
	ace$assigned_archetype = unification.out$assigned_archetype

	# Use graph core of global and induced subgraphs to infer centrality/quality of each cell
	ace$node_centrality = compute_archetype_core_centrality(G, ace$assigned_archetype)
    

	# Compute gene specificity for each archetype	
	if(is.matrix(S)) {
		specificity.out = compute_archetype_feature_specificity_full(S, unification.out$H_unified)
	} else {
		specificity.out = compute_archetype_feature_specificity(S, unification.out$H_unified)
	}
	
	specificity.out = lapply(specificity.out, function(specificity.scores) {
		rownames(specificity.scores) = rownames(ace)
		colnames(specificity.scores) = paste("A", 1:ncol(specificity.scores), sep = "")
		return(specificity.scores)
	})
	rowFactors(ace)[["H_unified_profile"]] = specificity.out[["archetypes"]]
	rowFactors(ace)[["H_unified_upper_significance"]] = specificity.out[["upper_significance"]]
	rowFactors(ace)[["H_unified_lower_significance"]] = specificity.out[["lower_significance"]]


	if(full.trace == T) {
		# Prepare output
		trace = list(ACTION.out = ACTION.out, pruning.out = pruning.out, vis.out = vis.out, unification.out = unification.out)
		trace$log = list(genes = rownames(ace), cells = colnames(ace), time = Sys.time())
			 
		out = list(ace = ace, trace = trace)
		
		return(out)
	} else {
		return(ace)
	}
}


#' Reconstructs the ACTIONet graph with the new parameters (uses prior decomposition)
#'
#' @param ace ACTIONetExperiment object containing the results
#' @param network_density Density factor of ACTIONet graph (default=1).
#' @param mutual_edges_only Whether to enforce edges to be mutually-nearest-neighbors (default=TRUE).
#' @param compactness_level A value between 0-100, indicating the compactness of ACTIONet layout (default=50)
#' @param n_epochs Number of epochs for SGD algorithm (default=500).
#' @param thread_no Number of parallel threads (default=0)
#' @param reduction.slot Slot in the colFactors(ace) that holds reduced kernel (default="S_r")
#' 
#' @return ace Updated ace object
#' 
#' @examples
#' plot.ACTIONet(ace)
#' ace.updated = reconstruct.ACTIONet(ace, network_density = 0.1)
#' plot.ACTIONet(ace.updated)
reconstruct.ACTIONet <- function(ace, network_density = 1, mutual_edges_only = TRUE, layout_compactness = 50, layout_epochs = 500, thread_no = 0, layout.in.parallel = FALSE, reduction.slot = "ACTION2D") {
    set.seed(0)
	
    # re-Build ACTIONet
	H_stacked = as.matrix(colFactors(ace)[["H_stacked"]])
	
    G = build_ACTIONet(H_stacked = H_stacked, density = network_density, thread_no=thread_no, mutual_edges_only = mutual_edges_only)
	colNets(ace)$ACTIONet = G
	
	
    # Layout ACTIONet
	initial.coordinates = t(scale(ACTIONet::colFactors(ace)[[reduction.slot]]))
	if(layout.in.parallel == FALSE) {		
		vis.out = layout_ACTIONet(G, S_r = initial.coordinates, compactness_level = layout_compactness, n_epochs = layout_epochs, thread_no = 1)
    } else { # WARNING! This makes the results none reproducible
		vis.out = layout_ACTIONet(G, S_r = initial.coordinates, compactness_level = layout_compactness, n_epochs = layout_epochs, thread_no = thread_no)
	}
	
    colFactors(ace)$ACTIONet2D = vis.out$coordinates
    colFactors(ace)$ACTIONet3D = vis.out$coordinates_3D
    colFactors(ace)$denovo_color = vis.out$colors

	return(ace)
}



#' Rerun layout on the ACTIONet graph with new parameters
#'
#' @param ace ACTIONetExperiment object containing the results
#' @param compactness_level A value between 0-100, indicating the compactness of ACTIONet layout (default=50)
#' @param n_epochs Number of epochs for SGD algorithm (default=500).
#' @param thread_no Number of parallel threads (default=8)
#' @param reduction.slot Slot in the colFactors(ace) that holds reduced kernel (default="S_r")
#' 
#' @return ace Updated ace object
#' 
#' @examples
#' plot.ACTIONet(ace)
#' ace.updated = rerun.layout(ace, layout_compactness = 20)
#' plot.ACTIONet(ace.updated)
rerun.layout <- function(ace, layout_compactness = 50, layout_epochs = 500, thread_no = 1, reduction.slot = "ACTION2D") {
    G = colNets(ace)[["ACTIONet"]]
    	
    # re-Layout ACTIONet
    S_r = t(ACTIONet::colFactors(ace)[[reduction.slot]])
    
	initial.coordinates = t(scale(t(S_r)))
	vis.out = layout_ACTIONet(G, S_r = initial.coordinates, compactness_level = layout_compactness, n_epochs = layout_epochs, thread_no = thread_no)
	    
    colFactors(ace)$ACTIONet2D = vis.out$coordinates
    colFactors(ace)$ACTIONet3D = vis.out$coordinates_3D
    colFactors(ace)$denovo_color = vis.out$colors

	return(ace)
}


rerun.archetype.aggregation <- function(ace, resolution = 1, data.slot = "logcounts", reduction.slot = "ACTION", unified_suffix = "unified") {
	S = assays(ace)[[data.slot]]
    S_r = t(ACTIONet::colFactors(ace)[[reduction.slot]])
	C_stacked = Matrix::t(as.matrix(colFactors(ace)[["C_stacked"]]))
	H_stacked = as.matrix(colFactors(ace)[["H_stacked"]])
    G = colNets(ace)[["ACTIONet"]]
	
	unification.out = unify_archetypes(S_r, C_stacked, H_stacked, min_overlap = 0, resolution = resolution)

	R.utils::printf("resolution = %d -> %d states\n", resolution, length(unique(unification.out$assigned_archetype)))
	
	colFactors(ace)[[sprintf("H_%s", unified_suffix)]] = as(unification.out$H_unified, 'sparseMatrix')
	colFactors(ace)[[sprintf("C_%s", unified_suffix)]] = as(Matrix::t(unification.out$C_unified), 'sparseMatrix');
	colData(ace)[[sprintf("%s_archetypes", unified_suffix)]]  = unification.out$assigned_archetype

	# Use graph core of global and induced subgraphs to infer centrality/quality of each cell
	ace$node_centrality = compute_archetype_core_centrality(G, ace$assigned_archetype)
    

	# Compute gene specificity for each archetype	
	if(is.matrix(S)) {
		specificity.out = compute_archetype_feature_specificity_full(S, unification.out$H_unified)
	} else {
		specificity.out = compute_archetype_feature_specificity(S, unification.out$H_unified)
	}
	
	specificity.out = lapply(specificity.out, function(specificity.scores) {
		rownames(specificity.scores) = rownames(ace)
		colnames(specificity.scores) = paste("A", 1:ncol(specificity.scores), sep = "")
		return(specificity.scores)
	})
	
	rowFactors(ace)[[sprintf("H_%s_profile", unified_suffix)]] = specificity.out[["archetypes"]]
	rowFactors(ace)[[sprintf("H_%s_upper_significance", unified_suffix)]] = specificity.out[["upper_significance"]]
	rowFactors(ace)[[sprintf("H_%s_lower_significance", unified_suffix)]] = specificity.out[["lower_significance"]]
	
	return(ace)
}

regroup.archetypes <- function(ace, unification.resolution = 1, data.slot = "logcounts", reduction.slot = "ACTION") {
	S = assays(ace)[[data.slot]]
	S_r = Matrix::t(colFactors(ace)[[reduction.slot]])
	
	H_stacked = as.matrix(colFactors(ace)[["H_stacked"]])
	C_stacked = as.matrix(Matrix::t(colFactors(ace)[["C_stacked"]]))
	G = colNets(ace)$ACTIONet
	
	# Identiy equivalent classes of archetypes and group them together
	unification.out = unify_archetypes(S_r = S_r, C_stacked = C_stacked, H_stacked = H_stacked, min_overlap = 0, resolution = unification.resolution)

	colFactors(ace)[["H_unified"]] = as(unification.out$H_unified, 'sparseMatrix')
	colFactors(ace)[["C_unified"]] = as(Matrix::t(unification.out$C_unified), 'sparseMatrix');
	ace$assigned_archetype = unification.out$assigned_archetype

	# Use graph core of global and induced subgraphs to infer centrality/quality of each cell
	ace$node_centrality = compute_archetype_core_centrality(G, ace$assigned_archetype)
    

	# Compute gene specificity for each archetype	
	if(is.matrix(S)) {
		specificity.out = compute_archetype_feature_specificity_full(S, unification.out$H_unified)
	} else {
		specificity.out = compute_archetype_feature_specificity(S, unification.out$H_unified)
	}
	
	specificity.out = lapply(specificity.out, function(specificity.scores) {
		rownames(specificity.scores) = rownames(ace)
		colnames(specificity.scores) = paste("A", 1:ncol(specificity.scores), sep = "")
		return(specificity.scores)
	})
	rowFactors(ace)[["H_unified_profile"]] = specificity.out[["archetypes"]]
	rowFactors(ace)[["H_unified_upper_significance"]] = specificity.out[["upper_significance"]]
	rowFactors(ace)[["H_unified_lower_significance"]] = specificity.out[["lower_significance"]]

	return(ace)
}
