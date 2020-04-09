#' A wrapper function to call all main functions of the ACTIONet
#'
#' @param sce Reduced `SingleCellExperiment` object (output of reduce.sce() function).
#' @param k_max Maximum depth of decompositions (default=30).
#' @param min_specificity_z_threshold Defines the stringency of pruning nonspecific archetypes. 
#' The larger the value, the more archetypes will be filtered out (default=-1).
#' @param network_density Density factor of ACTIONet graph (default=1).
#' @param mutual_edges_only Whether to enforce edges to be mutually-nearest-neighbors (default=TRUE).
#' @param compactness_level A value between 0-100, indicating the compactness of ACTIONet layout (default=50)
#' @param n_epochs Number of epochs for SGD algorithm (default=500).
#' @param thread_no Number of parallel threads (default=4)
#' @param reduction_slot Slot in the ReducedDims(sce) that holds reduced kernel (default="S_r")
#' @param sce.data.attr Corresponding slot in the `sce` object the normalized counts (default="logcounts")
#' 
#' @return A named list: \itemize{
#' \item ACTION.out: Output of run_ACTION() with C and H traces.
#' \item pruning.out: Output of prune_archetypes() with C_stacked and H_stacked inside.
#' \item ACTIONet: ACTIONet graph, constructed by build_ACTIONet(), in the igraph format
#' \item vis.out: Output of the layout_ACTIONet() function. It includes coordinates and colors for ACTIONet nodes.
#' \item unification.out: Output of the unify_archetypes() function. Includes C_unified, H_unified, and sample_assignments.
#'}
#' 
#' @examples
#' sce = reduce(sce)
#' ACTIONet.out = run_ACTIONet(sce)
run.ACTIONet <- function(sce, k_max = 30, min_specificity_z_threshold = -1, network_density = 1, mutual_edges_only = TRUE, layout_compactness = 50, n_epochs = 500, thread_no = 4, sce.data.attr = "logcounts", reduction_slot = "S_r") {
    if (!(sce.data.attr %in% names(SummarizedExperiment::assays(sce)))) {
        R.utils::printf("Attribute %s is not an assay of the input sce\n", sce.data.attr)
        return()
    }
    
    # Run ACTION
	ACTION.out = run_ACTION(t(SingleCellExperiment::reducedDims(sce)[[reduction_slot]]), k_min = 2, k_max = k_max, thread_no = thread_no)
    
    # Prune nonspecific and/or unreliable archetypes
    pruning.out = prune_archetypes(ACTION.out$C, ACTION.out$H, min_specificity_z_threshold = min_specificity_z_threshold)

    
    # Build ACTIONet
    set.seed(0)
    G = build_ACTIONet(H_stacked = pruning.out$H_stacked, density = network_density, thread_no=thread_no, mutual_edges_only = mutual_edges_only)

    # Layout ACTIONet
	initial.coordinates = t(scale(reducedDims(sce)[[reduction_slot]]))
    vis.out = layout_ACTIONet(G, S_r = initial.coordinates, compactness_level = layout_compactness, n_epochs = n_epochs)
    
    # Construct igraph object
    ACTIONet = igraph::graph_from_adjacency_matrix(G, mode = "undirected", weighted = TRUE)
    
    coor = vis.out$coordinates
    coor3D = vis.out$coordinates_3D
    V(ACTIONet)$x = coor[, 1]
    V(ACTIONet)$y = coor[, 2]
    V(ACTIONet)$x3D = coor3D[, 1]
    V(ACTIONet)$y3D = coor3D[, 2]
    V(ACTIONet)$z3D = coor3D[, 3]
    V(ACTIONet)$color = rgb(vis.out$colors)

	# Identiy equivalent classes of archetypes and group them together
	unification.out = unify_archetypes(G, S_r, pruning.out$C_stacked, pruning.out$H_stacked)

    # compute_node_connectivity(ACTIONet, unification.out$sample_assignments)  
    
    
    # Re-normalize input (~gene expression) matrix and compute feature (~gene) specificity scores
    S = assays(sce)[[sce.data.attr]]
    norm.out = ACTIONet::renormalize_input_matrix(S, unification.out$sample_assignments)
	S.norm = norm.out$S_norm

	specificity.out = compute_feature_specificity(S.norm, H)
	specificity.out = sapply(specificity.out, function(specificity.scores) {
		rownames(specificity.scores) = rownames(sce)
		colnames(specificity.scores) = paste("A", 1:ncol(specificity.scores))
		return(specificity.scores)
	})
	specificity.scores = SingleCellExperiment(assays = specificity.out)
		
	
	trace = list(ACTION.out = ACTION.out, pruning.out = pruning.out, vis.out = vis.out, unification.out = unification.out)
    ACTIONet.out = list(S.norm = S.norm, ACTIONet = ACTIONet, specificity.scores = specificity.scores, trace = trace)    
      
    
    if( ('cell.hashtag' %in% names(colData(sce))) ) {
		cells = sce$cell.hashtag
	} else if ( !is.null(colnames(sce)) ) {
		cells = colnames(sce)
	} else {
		cells = as.character(1:ncol(sce))
	}

    ACTIONet.out$log = list(genes = rownames(sce), cells = cells, time = Sys.time())
    
    return(ACTIONet.out)
}


#' Report the top-rank marker genes for each archetype
#'
#' @param ACTIONet.out Main ACTIONet output object
#' @param top.genes Number of genes to return
#' 
#' @return A data.frame with each column representing the top-ranked genes for an archetype.
#' 
#' @examples
#' ACTIONet.out.pruned = get.archetype.markers(ACTIONet.out, 10)
get.archetype.markers <- function(ACTIONet.out, top.genes = 20) {
	cell.state.DE = SummarizedExperiment::assays(ACTIONet.out$specificity.scores)[["upper_significance"]]
	Top.genes = apply(cell.state.DE, 2, function(x) rownames(cell.state.DE )[order(x, decreasing = T)[1:top.genes]])
	
	df = as.data.frame(Top.genes)
	
	return(df)
}

#' Report the top-ranked marker genes for each arbitrary annotation
#'
#' @param ACTIONet.out Output list to be pruned
#' @param top.genes Number of genes to return
#' 
#' @return A data.frame with each column representing the top-ranked genes for an archetype.
#' 
#' @examples
#' ACTIONet.out.pruned = get.archetype.markers(ACTIONet.out, 10)
get.annotation.markers <- function(ACTIONet.out, annotation.name, top.genes = 20) {
	idx = which((names(ACTIONet.out$annotations) == annotation.name) | (sapply(ACTIONet.out$annotations, function(X) X$annotation.name == annotation.name)))
	if(length(idx) == 0) {
		R.utils::printf('Annotation %s not found\n', annotation.name)
		return()
	}
	
	R.utils::printf('Annotation found: name = %s, tag = %s\n', names(ACTIONet.out$annotations)[[idx]], ACTIONet.out$annotations[[idx]]$annotation.name)
	
	if(is.null(ACTIONet.out$annotations[[idx]]$DE.profile)) {
		print("Please run compute.annotations.feature.specificity() first")
		return()		
	}
	DE.profile = as.matrix(SummarizedExperiment::assays(ACTIONet.out$annotations[[idx]]$DE.profile)[["significance"]])

	Top.genes = apply(DE.profile, 2, function(x) rownames(DE.profile)[order(x, decreasing = T)[1:top.genes]])
	
	df = as.data.frame(Top.genes)

	return(df)
}

# impute.genes.using.archetype <- function(ACTIONet.out, genes, prune = FALSE) {
#     require(igraph)
#     
#     genes = intersect(unique(genes), rownames(ACTIONet.out$unification.out$cellstates.core))
# 
# 	
# 	Z = as.matrix(ACTIONet.out$unification.out$cellstates.core[genes, ])
# 	H = ACTIONet.out$unification.out$H.core
# 	
# 	imputed.gene.expression = t(Z %*% H)
# 	colnames(imputed.gene.expression) = genes
# 	
#     return(imputed.gene.expression)
# }
