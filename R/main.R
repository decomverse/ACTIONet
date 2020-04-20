#' A wrapper function to call all main functions of the ACTIONet
#'
#' @param sce Reduced `SingleCellExperiment (SCE)` object (output of reduce.ace() function).
#' @param k_max Maximum depth of decompositions (default=30).
#' @param min_specificity_z_threshold Defines the stringency of pruning nonspecific archetypes. 
#' The larger the value, the more archetypes will be filtered out (default=-1).
#' @param network_density Density factor of ACTIONet graph (default=1).
#' @param mutual_edges_only Whether to enforce edges to be mutually-nearest-neighbors (default=TRUE).
#' @param compactness_level A value between 0-100, indicating the compactness of ACTIONet layout (default=50)
#' @param n_epochs Number of epochs for SGD algorithm (default=500).
#' @param thread_no Number of parallel threads (default=4)
#' @param reduction.slot Slot in the ReducedDims(ace) that holds reduced kernel (default="S_r")
#' @param data.slot Corresponding slot in the `ace` object the normalized counts (default="logcounts")
#' @param renormalize.logcounts.slot Name of the new assay with updated logcounts adjusted using archetypes
#' If it is NULL, values of logcounts(sce) would be directly used without renormalization for computing speicificity scores
#' 
#' @return A named list: \itemize{
#' \item ace: ACTIONetExperiment object (derived from SingleCellExperiment)
#' \item ACTIONet.trace: Log of ACTIONet function calls 
#'}
#' 
#' @examples
#' sce = reduce(sce)
#' ACTIONet.out = run.ACTIONet(sce)
#' ace = ACTIONet.out$ace # main output
#' trace = ACTIONet.out$trace # for backup
run.ACTIONet <- function(sce, k_max = 30, AA_delta = 0.01, min_specificity_z_threshold = -1, network_density = 1, mutual_edges_only = TRUE, layout_compactness = 50, layout_epochs = 500, thread_no = 4, data.slot = "logcounts", reduction.slot = "S_r", renormalize.logcounts.slot = "logcounts_adjusted") {
    if (!(data.slot %in% names(assays(sce)))) {
        R.utils::printf("Attribute %s is not an assay of the input ace\n", data.slot)
        return()
    }
    
    ace = as(sce, "ACTIONetExperiment")
    
    # Run ACTION
	ACTION.out = run_ACTION(t(SingleCellExperiment::reducedDims(ace)[[reduction.slot]]), k_min = 2, k_max = k_max, thread_no = thread_no, AA_delta)
    
    # Prune nonspecific and/or unreliable archetypes
    pruning.out = prune_archetypes(ACTION.out$C, ACTION.out$H, min_specificity_z_threshold = min_specificity_z_threshold)
    
    # Build ACTIONet
    set.seed(0)
    G = build_ACTIONet(H_stacked = pruning.out$H_stacked, density = network_density, thread_no=thread_no, mutual_edges_only = mutual_edges_only)
	colNets(ace)$ACTIONet = G
	
	
    # Layout ACTIONet
	initial.coordinates = t(scale(reducedDims(ace)[[reduction.slot]]))
    vis.out = layout_ACTIONet(G, S_r = initial.coordinates, compactness_level = layout_compactness, n_epochs = layout_epochs)
    
    reducedDims(ace)$ACTIONet2D = vis.out$coordinates
    reducedDims(ace)$ACTIONet3D = vis.out$coordinates_3D
    ace$denovo_color = rgb(vis.out$colors)


	# Identiy equivalent classes of archetypes and group them together
	unification.out = unify_archetypes(G, S_r, pruning.out$C_stacked, pruning.out$H_stacked)
	colFactors(ace)[["archetype_footprint"]] = unification.out$H_unified
	colFactors(ace)[["archetype_cell_contributions"]] = t(unification.out$C_unified)
	ace$archetype_assignment = unification.out$sample_assignments

	# Use graph core of global and induced subgraphs to infer centrality/quality of each cell
	ace$node_centrality = compute_archetype_core_centrality(G, ace$archetype_assignment)

	
    
    # Re-normalize input (~gene expression) matrix and compute feature (~gene) specificity scores
    if(is.null(renormalize.logcounts.slot)) {
		renormalize.logcounts.slot = data.slot
	} else {
		S = assays(ace)[[data.slot]]
		norm.out = ACTIONet::renormalize_input_matrix(S, unification.out$sample_assignments)
		assays(ace)[[renormalize.logcounts.slot]] = norm.out$S_norm
	}
	
	# Compute gene specificity for each archetype	
	S.norm = assays(ace)[[renormalize.logcounts.slot]]
	specificity.out = compute_archetype_feature_specificity(S.norm, unification.out$H_unified)
	specificity.out = lapply(specificity.out, function(specificity.scores) {
		rownames(specificity.scores) = rownames(ace)
		colnames(specificity.scores) = paste("A", 1:ncol(specificity.scores), sep = "")
		return(specificity.scores)
	})
	rowFactors(ace)[["archetype_gene_profile"]] = specificity.out[["average_profile"]]
	rowFactors(ace)[["archetype_gene_specificity"]] = specificity.out[["upper_significance"]]
	
	
	# Prepare output
	trace = list(ACTION.out = ACTION.out, pruning.out = pruning.out, vis.out = vis.out, unification.out = unification.out)
    trace$log = list(genes = rownames(ace), cells = colnames(ace), time = Sys.time())
      
    out = list(ace = ace, trace = trace)
    
    return(out)
}


#' Report the top-rank features from a given enrichment table
#'
#' @param top.features Number of features to return
#' @param reorder.columns Whether to optimally re-order columns of the enrichment table
#' 
#' @return Sorted table with the selected top-ranked
#' 
#' @examples
#' feature.enrichment.table = as.matrix(rowFactors(ace)[["archetype_gene_specificity"]])
#' enrichment.table.top = select.top.k.features(feature.enrichment.table, 3)
select.top.k.features <- function(feature.enrichment.table, top.features = 3, normalize = F, reorder.columns = T) {
	
	W0 = (feature.enrichment.table)
	if(normalize == T)
		W0 = doubleNorm(W0)
	
	IDX = matrix(0, nrow = top.features, ncol = ncol(W0))
	VV = matrix(0, nrow = top.features, ncol = ncol(W0))
	W = (W0)
	for(i in 1:nrow(IDX)) {
	  W.m = as(MWM_hungarian(W), 'dgTMatrix')
	  IDX[i, W.m@j+1] = W.m@i + 1
	  VV[i, W.m@j+1] = W.m@x
	  
	  W[IDX[i, W.m@j+1], ] = 0
	}

	if(reorder.columns == T) {
		feature.enrichment.table.aggregated = apply(IDX, 2, function(perm) as.numeric(Matrix::colMeans(W0[perm, ])))
		CC = cor(feature.enrichment.table.aggregated)
		D = as.dist(1 - CC)
		cols = seriation::get_order(seriation::seriate(D, "OLO"))	
		rows = as.numeric(IDX[, cols])
	} else {
		cols = 1:ncol(W0)
		rows = unique(as.numeric(IDX))
	}
	W = feature.enrichment.table[rows, cols]
	
	return(W)
}


