preprocessing(pp).*
	reduce_kernel
	SVD2ACTION
	SVD2PCA
	renormalize_input_matrix
	
decomposition
	run_ACTION
	run_AA
	run_SPA
	run_simplex_regression
	
	prune_archetypes
	unify_archetypes -- gives clusters
	
graph
	build_ACTIONet
	
graphtools.*
	embedding/visualization
		UMAP/LargeVis: layout_ACTIONet -- Add LargeVis interface, potentially update denovo colors
	diffusion/imputation:
		PR: compute_network_diffusion
		L1-reg PR: compute_sparse_network_diffusion
	clustering (Leiden+arch)
	coreness: compute_archetype_core_centrality
		

statistics:
	compute_archetype_feature_specificity
	compute_cluster_feature_specificity
	assess_enrichment
