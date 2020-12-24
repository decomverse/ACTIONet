compute_pairwise_alignment <- function(reference_profile, query_profile, reduced_dim = 50, 
    log_transform = FALSE, deflate = T, compute.Z = FALSE, sample.no = 1000, seed = 0) {
    set.seed(seed)
    
    if (log_transform == T) {
        if (max(query_profile) > 100) {
            query_profile = log1p(query_profile)
        }
        
        if (max(reference_profile) > 100) {
            reference_profile = log1p(reference_profile)
        }
    }
    
    query_profile = scale(query_profile, center = F)
    reference_profile = scale(reference_profile, center = F)
    
    # Deflate?
    if (deflate == TRUE) {
        query_profile.red = reduce_kernel_full(query_profile, reduced_dim = reduced_dim, 
            iter = 1000)
        reference_profile.red = reduce_kernel_full(reference_profile, reduced_dim = reduced_dim, 
            iter = 1000)
        
        
        
        query_profile.deflated = query_profile + query_profile.red$A %*% t(query_profile.red$B)
        reference_profile.deflated = reference_profile + reference_profile.red$A %*% 
            t(reference_profile.red$B)
        
        query_profile_centered = t(scale(t(query_profile.deflated), scale = FALSE))
        reference_profile_centered = t(scale(t(reference_profile.deflated), scale = FALSE))
    } else {
        query_profile_centered = t(scale(t(query_profile), scale = FALSE))
        reference_profile_centered = t(scale(t(reference_profile), scale = FALSE))
    }
    
    reduced_dim = min(min(min(dim(query_profile)), min(dim(reference_profile))) - 
        1, reduced_dim)
    
    # pc.out1 = prcomp(t(query_profile_centered), retx = TRUE, rank = reduced_dim)
    # pc.out2 = prcomp(t(reference_profile_centered), retx = TRUE, rank =
    # reduced_dim) V.query = pc.out1$rotation[, 1:reduced_dim] V.reference =
    # pc.out2$rotation[, 1:reduced_dim]
    
    query_profile_centered[is.na(query_profile_centered)] = 0
    reference_profile_centered[is.na(reference_profile_centered)] = 0
    
    query_profile.red = ACTIONet::IRLB_SVD_full(t(query_profile_centered), dim = reduced_dim, 
        iters = 100)
    reference_profile.red = ACTIONet::IRLB_SVD_full(t(reference_profile_centered), 
        dim = reduced_dim, iters = 100)
    
    V.query = query_profile.red$v[, 1:reduced_dim]
    V.reference = reference_profile.red$v[, 1:reduced_dim]
    
    V.alignment = V.query %*% (t(V.query) %*% V.reference)
    S_r.query = t(V.alignment) %*% query_profile_centered
    S_r.reference = t(V.reference) %*% reference_profile_centered
    
    X = t(run_simplex_regression(S_r.query, S_r.reference))
    
    return(X)
}


compute_merged_ace_from_cell_alignment <- function(reference_ace, query_ace, cell_to_cell_alignment, 
    reference_cell_labels = NULL, query_cell_labels = NULL, n_epochs = 500) {
    
    if (is.null(colnames(reference_ace))) 
        colnames(reference_ace) = paste("refCell", 1:ncol(reference_ace))
    
    if (is.null(colnames(query_ace))) 
        colnames(query_ace) = paste("refCell", 1:ncol(query_ace))
    
    proj.out = transform_layout(cell_to_cell_alignment, t(reference_ace$ACTIONet2D), 
        t(reference_ace$ACTIONet3D), t(reference_ace$denovo_color), n_epochs = n_epochs)
    
    combined.ace = ACTIONetExperiment(rowData = DataFrame(Id = c(rownames(reference_ace), 
        rownames(query_ace))), colData = DataFrame(Id = c(colnames(reference_ace), 
        colnames(query_ace)), dataset = c(rep("reference", ncol(reference_ace)), 
        rep("query", ncol(query_ace)))))
    
    if (!is.null(reference_cell_labels) & !is.null(query_cell_labels)) {
        combined.ace$Labels = c(as.character(reference_cell_labels), as.character(query_cell_labels))
    }
    
    
    X = rbind(reference_ace$ACTIONet2D, Matrix::t(proj.out$coordinates))
    colnames(X) = c("x", "y")
    rownames(X) = colnames(combined.ace)
    colMaps(combined.ace)$ACTIONet2D = X
    colMapTypes(combined.ace)[["ACTIONet2D"]] = "embedding"
    X = rbind(reference_ace$ACTIONet3D, Matrix::t(proj.out$coordinates_3D))
    colnames(X) = c("x", "y", "z")
    rownames(X) = colnames(combined.ace)
    colMaps(combined.ace)$ACTIONet3D = X
    colMapTypes(combined.ace)[["ACTIONet3D"]] = "embedding"
    X = rbind(reference_ace$denovo_color, Matrix::t(proj.out$colors))
    colnames(X) = c("r", "g", "b")
    rownames(X) = colnames(combined.ace)
    colMaps(combined.ace)$denovo_color = X
    colMapTypes(combined.ace)[["denovo_color"]] = "embedding"
    
    return(combined.ace)
}


plot_pairwise_alignment <- function(alignment, reference_labels = NULL, query_labels = NULL, 
    normalize = FALSE) {
    W = alignment
    W[W < 0] = 0
    Wm = as(MWM_hungarian(W), "dgTMatrix")
    
    ii = Wm@i + 1
    jj = Wm@j + 1
    perm = order(Wm@x, decreasing = T)
    subW = alignment[ii, jj]
    rownames(subW) = reference_labels[ii]
    colnames(subW) = query_labels[jj]
    
    library(seriation)
    perm = get_order(seriate(as.dist(1 - cor(subW)), "OLO"))
    subW = subW[perm, perm]
    
    gradPal = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu"))))(100)
    
    if (length(alignment) < 10000) 
        border.color = "black" else border.color = NA
    
    X = subW
    if (normalize == TRUE) 
        X = doubleNorm(X)
    
    ht = Heatmap(X, cluster_columns = F, cluster_rows = F, rect_gp = gpar(col = border.color), 
        name = "Alignment score", col = gradPal, row_names_side = "left", row_title = "Reference", 
        column_title = "Query", column_names_gp = gpar(fontsize = 14, fontface = "bold"), 
        row_names_gp = gpar(fontsize = 14, fontface = "bold"), column_title_gp = gpar(fontsize = 18, 
            fontface = "bold"), row_title_gp = gpar(fontsize = 18, fontface = "bold"), 
        row_names_max_width = unit(100, "cm"), column_names_max_height = unit(100, 
            "cm"))
    
    ht
}


compute_cell_alignments_from_archetype_alignments <- function(reference_ace, query_ace, 
    alignment, alignment_threshold = 0.1, footprint_threshold = 0.1, reference_slot_name = "unified", 
    query_slot_name = "unified") {
    
    W = alignment
    W[W < alignment_threshold] = 0
    Wm = as(MWM_hungarian(W), "dgTMatrix")
    
    ii = Wm@i + 1
    jj = Wm@j + 1
    w = Wm@x
    
    reference_footprint = colMaps(reference_ace)[[sprintf("H_%s", reference_slot_name)]]
    reference_footprint = reference_footprint[, ii]
    
    query_footprint = colMaps(query_ace)[[sprintf("H_%s", query_slot_name)]]
    query_footprint = query_footprint[, jj]
    
    
    cell_to_cell_alignments = lapply(1:length(ii), function(j) {
        u = reference_footprint[, j]
        u = u/max(u)
        v = query_footprint[, j]
        v = v/max(v)
        match.out = MWM_rank1(u = u, v, footprint_threshold, footprint_threshold)
        
        match.ii = match.out[1, ]
        match.jj = match.out[2, ]
        match.w = (w[j] * u[match.out[1, ]] * v[match.out[2, ]] * w[j])^(1/3)
        # df = data.frame(RNA.ace$Celltype[match.ii], ATAC.ace$Celltype[match.jj],
        # match.w)
        
        M = sparseMatrix(i = match.ii, j = match.jj, x = match.w, dims = c(length(u), 
            length(v)))
        return(M)
    })
    
    cell_to_cell_alignment = Reduce("+", cell_to_cell_alignments)
    
    cell_to_cell_alignment@x[cell_to_cell_alignment@x > 1] = 1
    
    return(cell_to_cell_alignment)
    
}


annotate_cells_from_alignment <- function(ace, alignment, unification.slot = "H_unified") {
    
    cell.scores.mat = colMaps(ace)[[unification.slot]]
    cell.enrichment.mat = cell.scores.mat %*% Matrix::t(alignment)
    
    newLabels = colnames(cell.enrichment.mat)[apply(cell.enrichment.mat, 1, which.max)]
    newLabels.confidence = apply(cell.enrichment.mat, 1, max)
    
    newLabels.corrected = correct.cell.annotations(ace, newLabels, LFR.threshold = 1)
    
    newLabels.annot = annotate.archetypes.using.labels(ace, newLabels.corrected)
    
    M = as(MWM_hungarian(newLabels.annot$Enrichment), "dgTMatrix")
    
    newLabels.CPal = colorspace::lighten(rgb(metadata(ace)$backbone$colors[M@i + 
        1, ]), 0.25)
    names(newLabels.CPal) = colnames(newLabels.annot$Enrichment)[M@j + 1]
    
    
    
    
    out = list(Labels = newLabels.corrected, Labels.confidence = newLabels.confidence, 
        Enrichment = cell.enrichment.mat, Pal = newLabels.CPal)
    
    
    return(out)
}
