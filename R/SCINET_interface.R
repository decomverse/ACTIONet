run.SCINET.archetype <- function(ACTIONet.out, G=NULL, core = T, min.edge.weight = 2, spec.sample_no = 1000, thread_no = 8, compute.topo.specificity = T) {
	require(SCINET)
	
	print("Preprocessing the baseline interactome");
	if( is.null(G) ) {
		if(!exists('PCNet') ) {
			data("PCNet")
		}
		Adj = PCNet		
	} else if(is.matrix(G) | is.sparseMatrix(G)) {
		Adj = as(G, 'sparseMatrix')
		Adj@x = rep(1, length(Adj@x))
	} else if(is.igraph(G)){		
		Adj = as(get.adjacency(G), 'sparseMatrix')
	}
      
	if(core == T) {
		if (("unification.out" %in% names(ACTIONet.out))) {
			print("Using unification.out$DE.core (merged archetypes)")
			DE.profile = as.matrix(log1p(SummarizedExperiment::assays(ACTIONet.out$unification.out$DE.core)[["significance"]]))
		} else {
			print("unification.out is not in ACTIONet.out. Running unify.cell.states() first ...")
			ACTIONet.out$unification.out = unify.cell.states(ACTIONet.out, sce, reduction_slot = reduction_slot, sce.data.attr = sce.data.attr)    
			DE.profile = as.matrix(log1p(SummarizedExperiment::assays(ACTIONet.out$unification.out$DE.core)[["significance"]]))
		}
	} else {
		if (("archetype.differential.signature" %in% names(ACTIONet.out))) {
			print("Using archetype.differential.signature (all archetypes)")
			DE.profile = as.matrix(log1p(ACTIONet.out$archetype.differential.signature))
		} else {
			print("archetype.differential.signature is not in ACTIONet.out. Running compute.archetype.feature.specificity() first ...")
			ACTIONet.out$archetype.differential.signature = compute.archetype.feature.specificity(ACTIONet.out, sce, mode = specificity.mode, sce.data.attr = sce.data.attr)
			DE.profile = as.matrix(log1p(ACTIONet.out$archetype.differential.signature))
		}
	}                  
	
	
	common.genes = intersect(rownames(DE.profile), rownames(PCNet))
	if(length(common.genes) == 0) {
		print("No common genes found. Check rownames (or vertex names) for the input graph");
		return(ACTIONet.out)
	}
	A = DE.profile[common.genes, ]
	G = Adj[common.genes, common.genes]
	
	
	print("Constructing networks");
	gene.activity.scores = SCINET::compute_gene_activities_full(A = A, thread_no = thread_no)
	cellstate.nets = SCINET::construct_cell_networks(net = G, gene_activities = gene.activity.scores, thread_no = thread_no)
	cellstate.nets.list = as.list(cellstate.nets)
	
	print("Prunning nodes/edges ... ");
	cellstate.nets.list.igraph = lapply(cellstate.nets.list, function(G.Adj) {
		G.Adj@x[G.Adj@x < min.edge.weight] = 0
		filter.mask = Matrix::colSums(G.Adj) == 0
		G = igraph::graph_from_adjacency_matrix(G.Adj[!filter.mask, !filter.mask], mode = "undirected", weighted = T)
		V(G)$name = common.genes[!filter.mask]
		if(compute.topo.specificity == TRUE) {
			z.scores = topo.spec(G, spec.sample_no)
			V(G)$specificity = 1 / (1 + exp(-z.scores))
		}
		
		# G = igraph::graph_from_adjacency_matrix(G.Adj, mode = "undirected", weighted = T)
		# V(G)$name = common.genes
		# G = delete_edges(G, E(G)[E(G)$weight < min.edge.weight])
		# if(compute.topo.specificity == TRUE) {
		# 	z.scores = topo.spec(G, spec.sample_no)
		# 	V(G)$specificity = 1 / (1 + exp(-z.scores))
		# }
		# G = delete_vertices(G, V(G)[strength(G) == 0])				
		
		return(G)
	})

	if(is.null(colnames(DE.profile))) {		
		names(cellstate.nets.list.igraph) = 1:ncol(DE.profile)
	} else {
		names(cellstate.nets.list.igraph) = colnames(DE.profile)
	}
	
	return(cellstate.nets.list.igraph)
}

run.SCINET.annotation <- function(ACTIONet.out, annotation_name, G=NULL, min.edge.weight = 2, spec.sample_no = 1000, thread_no = 8, compute.topo.specificity = T) {
	require(SCINET)
	
	print("Preprocessing the baseline interactome");
	if( is.null(G) ) {
		if(!exists('PCNet') ) {
			data("PCNet")
		}
		Adj = PCNet		
	} else if(is.matrix(G) | is.sparseMatrix(G)) {
		Adj = as(G, 'sparseMatrix')
		Adj@x = rep(1, length(Adj@x))
	} else if(is.igraph(G)){		
		Adj = as(get.adjacency(G), 'sparseMatrix')
	}


	cl.idx = which(names(ACTIONet.out$annotations) == annotation_name)
	if(length(cl.idx) == 0) {
		R.utils::printf('Error in run.annotation.SCINET: annotation_name "%s" not found\n', annotation_name)
		return(ACTIONet.out)
	}		
	
	if( is.null(ACTIONet.out$annotations[[cl.idx]]$DE.profile) ) {
		print("DE.profile of the annotation is missing. Computing it from scratch (please wait) ... ");		
		ACTIONet.out = compute.annotations.feature.specificity(ACTIONet.out, sce, annotation_name)
	}
	DE.profile  = log1p(as.matrix(SummarizedExperiment::assays(ACTIONet.out$annotations[[cl.idx]]$DE.profile)[["significance"]]))
	
	common.genes = intersect(rownames(DE.profile), rownames(PCNet))
	if(length(common.genes) == 0) {
		print("No common genes found. Check rownames (or vertex names) for the input graph");
		return(ACTIONet.out)
	}
	A = DE.profile[common.genes, ]
	G = Adj[common.genes, common.genes]
	
	
	print("Constructing networks");
	gene.activity.scores = SCINET::compute_gene_activities_full(A = A, thread_no = thread_no)
	cellstate.nets = SCINET::construct_cell_networks(net = G, gene_activities = gene.activity.scores, thread_no = thread_no)
	cellstate.nets.list = as.list(cellstate.nets)
	
	print("Prunning nodes/edges ... ");
	cellstate.nets.list.igraph = lapply(cellstate.nets.list, function(G.Adj) {
		G.Adj@x[G.Adj@x < min.edge.weight] = 0
		filter.mask = Matrix::colSums(G.Adj) == 0
		G = igraph::graph_from_adjacency_matrix(G.Adj[!filter.mask, !filter.mask], mode = "undirected", weighted = T)
		V(G)$name = common.genes[!filter.mask]
		if(compute.topo.specificity == TRUE) {
			z.scores = topo.spec(G, spec.sample_no)
			V(G)$specificity = 1 / (1 + exp(-z.scores))
		}
		
		# G = igraph::graph_from_adjacency_matrix(G.Adj, mode = "undirected", weighted = T)
		# V(G)$name = common.genes
		# G = delete_edges(G, E(G)[E(G)$weight < min.edge.weight])
		# if(compute.topo.specificity == TRUE) {
		# 	z.scores = topo.spec(G, spec.sample_no)
		# 	V(G)$specificity = 1 / (1 + exp(-z.scores))
		# }
		# G = delete_vertices(G, V(G)[strength(G) == 0])
						
		return(G)
	})

	if(is.null(colnames(DE.profile))) {		
		names(cellstate.nets.list.igraph) = 1:ncol(DE.profile)
	} else {
		names(cellstate.nets.list.igraph) = colnames(DE.profile)
	}
	return(cellstate.nets.list.igraph)
}


run.SCINET.gene.scores <- function(gene.scores.mat, G=NULL, min.edge.weight = 2, spec.sample_no = 1000, thread_no = 8, compute.topo.specificity = T) {
	require(SCINET)
	
	print("Preprocessing the baseline interactome");
	if( is.null(G) ) {
		if(!exists('PCNet') ) {
			data("PCNet")
		}
		Adj = PCNet		
	} else if(is.matrix(G) | is.sparseMatrix(G)) {
		Adj = as(G, 'sparseMatrix')
		Adj@x = rep(1, length(Adj@x))
	} else if(is.igraph(G)){		
		Adj = as(get.adjacency(G), 'sparseMatrix')
	}


	common.genes = intersect(rownames(gene.scores.mat), rownames(Adj))
	if(length(common.genes) == 0) {
		print("No common genes found. Check rownames (or vertex names) for the input graph");
		return()
	}
	A = gene.scores.mat[common.genes, ]
	G = Adj[common.genes, common.genes]
	
	
	
	print("Constructing networks");
	gene.activity.scores = SCINET::RIN_transform(A = A, thread_no = thread_no)
	cellstate.nets = SCINET::construct_cell_networks(net = G, gene_activities = gene.activity.scores, thread_no = thread_no)
	cellstate.nets.list = as.list(cellstate.nets)
	
	print("Prunning nodes/edges ... ");
	cellstate.nets.list.igraph = lapply(cellstate.nets.list, function(G.Adj) {
		G.Adj@x[G.Adj@x < min.edge.weight] = 0
		filter.mask = Matrix::colSums(G.Adj) == 0
		G = igraph::graph_from_adjacency_matrix(G.Adj[!filter.mask, !filter.mask], mode = "undirected", weighted = T)
		V(G)$name = common.genes[!filter.mask]
		if(compute.topo.specificity == TRUE) {
			z.scores = topo.spec(G, spec.sample_no)
			V(G)$specificity.z = z.scores
			V(G)$specificity = 1 / (1 + exp(-z.scores))
		}
		# G = igraph::graph_from_adjacency_matrix(G.Adj, mode = "undirected", weighted = T)
		# V(G)$name = common.genes
		# G = delete_edges(G, E(G)[E(G)$weight < min.edge.weight])
		# if(compute.topo.specificity == TRUE) {
		# 	z.scores = topo.spec(G, spec.sample_no)
		# 	V(G)$specificity = 1 / (1 + exp(-z.scores))
		# }
		# G = delete_vertices(G, V(G)[strength(G) == 0])				
		return(G)
	})
	
	if(is.null(colnames(gene.scores.mat))) {		
		names(cellstate.nets.list.igraph) = 1:ncol(gene.scores.mat)
	} else {
		names(cellstate.nets.list.igraph) = colnames(gene.scores.mat)
	}

	return(cellstate.nets.list.igraph)
}
