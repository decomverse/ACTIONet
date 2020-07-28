compute_RNA_archetype_to_RNA_archetype_alignment <- function(reference_ace, query_ace, reference_slot_name = "unified", query_slot_name = "unified", deflate = FALSE, reduced_dim = 50, specificity_filter_threshold = 1) {
	reference_profile = rowMaps(reference_ace)[[sprintf("%s_feature_specificity", reference_slot_name)]]
	query_profile = rowMaps(query_ace)[[sprintf("%s_feature_specificity", query_slot_name)]]

	g1 = rownames(reference_ace)[apply(reference_profile, 1, max) > specificity_filter_threshold]
	g2 = rownames(query_ace)[apply(query_profile, 1, max) > specificity_filter_threshold]
	common.genes = intersect(g1, g2)

	reference_profile = reference_profile[common.genes, ]
	query_profile = query_profile[common.genes, ]

	alignment = compute_pairwise_alignment(reference_profile = reference_profile, query_profile = query_profile, reduced_dim = reduced_dim, deflate = deflate)
	
	return(alignment)
}


compute_RNA_cell_to_RNA_cell_alignment <- function(reference_ace, query_ace, archetype_alignment = NULL, specificity_filter_threshold = 1, alignment_threshold = 0.1, footprint_threshold = 0.1, reference_slot_name = "unified", query_slot_name = "unified", deflate = FALSE, reduced_dim = 50) {	
	if(is.null(archetype_alignment)) {
		archetype_alignment = compute_RNA_archetype_to_RNA_archetype_alignment(reference_ace = reference_ace, query_ace = query_ace, reference_slot_name = reference_slot_name, query_slot_name = query_slot_name, deflate = deflate, reduced_dim = reduced_dim, specificity_filter_threshold = specificity_filter_threshold)
	}
	
	cell_to_cell_alignment = compute_cell_alignments_from_archetype_alignments(reference_ace = reference_ace, query_ace = query_ace, alignment = archetype_alignment, alignment_threshold = alignment_threshold, reference_slot_name = reference_slot_name, query_slot_name = query_slot_name, alignment_threshold = alignment_threshold, footprint_threshold = footprint_threshold)

	return(cell_to_cell_alignment)		
	
	
}

compute_bulkRNA_to_RNA_archetype_alignment <- function(bulk, query_ace, bulk_assay_slot = "logcounts", query_slot_name = "unified", deflate = FALSE, reduced_dim = 50, specificity_filter_threshold = 1) {
	reference_profile = assays(bulk)[[bulk_assay_slot]]
	query_profile = rowMaps(query_ace)[[sprintf("%s_feature_specificity", query_slot_name)]]

	filtered_query_genes = rownames(query_ace)[apply(query_profile, 1, max) > specificity_filter_threshold]

	common.genes = intersect(rownames(bulk), filtered_query_genes)
	reference_profile = reference_profile[common.genes, ]
	query_profile = query_profile[common.genes, ]

	alignment = compute_pairwise_alignment(reference_profile = reference_profile, query_profile = query_profile, reduced_dim = reduced_dim, deflate = deflate)
	
	return(alignment)	
	
}
