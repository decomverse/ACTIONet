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
	
	cell_to_cell_alignment = compute_cell_alignments_from_archetype_alignments(reference_ace = reference_ace, query_ace = query_ace, alignment = archetype_alignment, reference_slot_name = reference_slot_name, query_slot_name = query_slot_name, alignment_threshold = alignment_threshold, footprint_threshold = footprint_threshold)

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

bidirectional_RNA_archetype_to_RNA_archetype_alignment <- function(reference_ace, query_ace, reference_labels = NULL, query_labels = NULL, reference_name = "Reference") {
	gradPal = circlize::colorRamp2(seq(0, 0.9, by=0.1), rev(pals::brewer.spectral(10)))
	
	alignment.forward = compute_RNA_archetype_to_RNA_archetype_alignment(reference_ace, query_ace)
	alignment.backward = Matrix::t(compute_RNA_archetype_to_RNA_archetype_alignment(query_ace, reference_ace))
	alignment = pmin(alignment.forward, alignment.backward)


	rownames(alignment) = rownames(alignment.backward) = rownames(alignment.forward) = paste("A", 1:nrow(alignment), "-", reference_labels, sep = "")
	colnames(alignment) = colnames(alignment.backward) = colnames(alignment.forward) = paste("A", 1:ncol(alignment), "-", query_labels, sep = "")
	
	
	W = alignment
	W[W < 0] = 0
	Wm = as(MWM_hungarian(W), 'dgTMatrix')

	ii = Wm@i+1
	jj = Wm@j+1
	perm = order(Wm@x, decreasing = T)
	subW = alignment[ii, jj]
	subW.forward = alignment.forward[ii, jj]
	subW.backward = alignment.backward[ii, jj]


	library(seriation)
	# CC = (cor(subW) + cor(subW.forward) + cor(subW.backward))/3
	CC = cor(subW)
	perm = get_order(seriate(as.dist(1-cor(CC)), "OLO"))
	subW = subW[perm, perm]
	subW.backward = subW.backward[perm, perm]
	subW.forward = subW.forward[perm, perm]

	gradPal = circlize::colorRamp2(seq(0, 0.9, by=0.1), rev(pals::brewer.spectral(10)))

	ht = Heatmap(subW, cluster_columns = F, cluster_rows = F, rect_gp = gpar(col = "black"), name = "mutual", col = gradPal, row_names_side = "left", row_title = reference_name, column_title = "Mutual", column_names_gp = gpar(fontsize = 14, fontface="bold"), row_names_gp = gpar(fontsize = 14, fontface="bold"), column_title_gp = gpar(fontsize = 18, fontface="bold"), row_title_gp = gpar(fontsize = 18, fontface="bold")) + Heatmap(subW.forward, cluster_columns = F, cluster_rows = F, rect_gp = gpar(col = "black"), name = "forward", col = gradPal, row_names_side = "left", row_title = reference_name, column_title = "Forward", column_names_gp = gpar(fontsize = 14, fontface="bold"), row_names_gp = gpar(fontsize = 14, fontface="bold"), column_title_gp = gpar(fontsize = 18, fontface="bold"), row_title_gp = gpar(fontsize = 18, fontface="bold")) + Heatmap(subW.backward, cluster_columns = F, cluster_rows = F, rect_gp = gpar(col = "black"), name = "backward", col = gradPal, row_names_side = "left", row_title = reference_name, column_title = "Backward", column_names_gp = gpar(fontsize = 14, fontface="bold"), row_names_gp = gpar(fontsize = 14, fontface="bold"), column_title_gp = gpar(fontsize = 18, fontface="bold"), row_title_gp = gpar(fontsize = 18, fontface="bold"), row_names_max_width = unit(100, "cm"), column_names_max_height = unit(100, "cm")) 
}

