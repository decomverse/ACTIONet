#' Filter columns and rows of `ACTIONetExperiment` or `SummarizedExperiment` object.

filter.ace <- function(ace, assay.name = "counts", min_cells_per_feat = NULL, min_feats_per_cell = NULL, min_umis_per_cell = NULL, max_umis_per_cell = NULL, 
    return_fil_ace = TRUE) {
    org_dim = dim(ace)
    ace.fil = ace
    
    i = 0
    repeat {
        prev_dim = dim(ace.fil)
        
        rows_mask = rep(TRUE, nrow(ace.fil))
        cols_mask = rep(TRUE, ncol(ace.fil))
        if (!is.null(min_umis_per_cell)) {
            umi_mask = Matrix::colSums(assays(ace.fil)[[assay.name]]) >= min_umis_per_cell
            cols_mask = cols_mask & umi_mask
        }
        
        if (!is.null(max_umis_per_cell)) {
            umi_mask = Matrix::colSums(assays(ace.fil)[[assay.name]]) <= max_umis_per_cell
            cols_mask = cols_mask & umi_mask
        }
        
        if (!is.null(min_feats_per_cell)) {
            feature_mask = Matrix::colSums(assays(ace.fil)[[assay.name]] > 0) >= min_feats_per_cell
            cols_mask = cols_mask & feature_mask
        }
        
        if (!is.null(min_cells_per_feat)) {
            if ((min_cells_per_feat < 1) & (min_cells_per_feat > 0)) {
                min_cells_per_feat = min_cells_per_feat * org_dim[2]
            }
            
            cell_count_mask = Matrix::rowSums(assays(ace.fil)[[assay.name]] > 0) >= min_cells_per_feat
            rows_mask = rows_mask & cell_count_mask
        }
        ace.fil <- ace.fil[rows_mask, cols_mask]
        
        i = i + 1
        if (all(dim(ace.fil) == prev_dim)) {
            break
        }
    }
    invisible(gc())
    
    if (return_fil_ace) 
        return(ace.fil) else {
        fil_cols_mask = !(colnames(ace) %in% colnames(ace.fil))
        fil_rows_mask = !(rownames(ace) %in% rownames(ace.fil))
        fil_cols_list = data.frame(name = colnames(ace)[fil_cols_mask], idx = which(fil_cols_mask))
        fil_rows_list = data.frame(name = rownames(ace)[fil_rows_mask], idx = which(fil_rows_mask))
        
        fil_list = list(cols_filtered = fil_cols_list, rows_filtered = fil_rows_list)
        return(fil_list)
    }
}

filter.ace.by.attr <- function(ace, by, assay_name = "counts", min_cells_per_feat = NULL, min_feats_per_cell = NULL, min_umis_per_cell = NULL, 
    max_umis_per_cell = NULL) {
    
    # if( length(by) == 1) { IDX = split(1:dim(ace)[2], droplevels(colData(ace))[[by]]) } else { IDX = split(1:dim(ace)[2], by) }
    IDX = get_ace_split_IDX(ace, by)
    
    if (any(duplicated(rownames(ace)))) {
        msg = sprintf("Adding suffix to duplicate rownames.")
        warning(msg)
        rownames(ace) = make.names(rownames(ace), unique = T)
    }
    if (any(duplicated(colnames(ace)))) {
        msg = sprintf("Adding suffix to duplicate colnames.")
        warning(msg)
        colnames(ace) = make.names(colnames(ace), unique = T)
    }
    
    fil_names <- lapply(IDX, function(idx) {
        fil_list <- filter.ace(ace[, idx], assay.name = assay_name, min_cells_per_feat = min_cells_per_feat, min_umis_per_cell = min_umis_per_cell, 
            max_umis_per_cell = max_umis_per_cell, min_feats_per_cell = min_feats_per_cell, return_fil_ace = F)
        return(fil_list)
    })
    
    fil_col = lapply(fil_names, function(i) i[["cols_filtered"]]$name) %>% Reduce(union, .)
    fil_row = lapply(fil_names, function(i) i[["rows_filtered"]]$name) %>% Reduce(intersect, .)
    keep_row = which(!(rownames(ace) %in% fil_row))
    keep_col = which(!(colnames(ace) %in% fil_col))
    
    ace.fil = ace[keep_row, keep_col]
    colData(ace.fil) <- droplevels(colData(ace.fil))
    rowData(ace.fil) <- droplevels(rowData(ace.fil))
    invisible(gc())
    return(ace.fil)
}
