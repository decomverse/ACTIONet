#' Filter columns and rows of `ACTIONetExperiment` or `SummarizedExperiment`-like object.
#' @export
filter.ace <- function(
  ace,
  assay_name = "counts",
  min_cells_per_feat = NULL,
  min_feats_per_cell = NULL,
  min_umis_per_cell = NULL,
  max_umis_per_cell = NULL,
  max_mito_fraction = NULL,
  species = c("hsapiens", "mmusculus"),
  features_use = NULL,
  return_fil_ace = TRUE
) {

    init_dim = dim(ace)
    init_dnames = dimnames(ace)

    X = SummarizedExperiment::assays(ace)[[assay_name]]
    X = as(X, "dgCMatrix")
    dimnames(X) = list(1:NROW(X), 1:NCOL(X))

    i = 0
    repeat {
        prev_dim = dim(X)
        rows_mask = rep(TRUE, NROW(X))
        cols_mask = rep(TRUE, NCOL(X))
        if (!is.null(min_umis_per_cell)) {
            umi_mask = Matrix::colSums(X) >= min_umis_per_cell
            cols_mask = cols_mask & umi_mask
        }

        if (!is.null(max_umis_per_cell)) {
            umi_mask = Matrix::colSums(X) <= max_umis_per_cell
            cols_mask = cols_mask & umi_mask
        }

        if (!is.null(min_feats_per_cell)) {
            feature_mask = Matrix::colSums(as(X > 0, "dgCMatrix")) >= min_feats_per_cell
            cols_mask = cols_mask & feature_mask
        }

        if (!is.null(min_cells_per_feat)) {
            if ((min_cells_per_feat < 1) & (min_cells_per_feat > 0)) {
                min_fc = min_cells_per_feat * prev_dim[2]
            } else {
                min_fc = min_cells_per_feat
            }
            cell_count_mask = ACTIONetExperiment:::fastRowSums(as(X > 0, "dgCMatrix")) >= min_fc
            rows_mask = rows_mask & cell_count_mask
        }

        X <- X[rows_mask, cols_mask]
        invisible(gc())
        i = i + 1
        if (all(dim(X) == prev_dim)) {
            break
        }
    }
    ace = ace[as.numeric(rownames(X)), as.numeric(colnames(X))]
    invisible(gc())

    if (!is.null(max_mito_fraction)){
      mt_frac = get_mtRNA_stats(
        ace,
        by = NULL,
        groups_use = NULL,
        assay = assay_name,
        species = match.arg(species),
        metric = "pct",
        features_use = features_use
      )

      ace = ace[, mt_frac <= max_mito_fraction]
      invisible(gc())
    }

    if (return_fil_ace){
      return(ace)
    } else {
      fil_cols_mask = !(init_dnames[[2]] %in% colnames(ace))
      fil_rows_mask = !(init_dnames[[1]] %in% rownames(ace))

      fil_cols_list = data.frame(
        name = init_dnames[[2]][fil_cols_mask],
        idx = which(fil_cols_mask)
      )

      fil_rows_list = data.frame(
        name = init_dnames[[1]][fil_rows_mask],
        idx = which(fil_rows_mask)
      )

      fil_list = list(
        cols_filtered = fil_cols_list,
        rows_filtered = fil_rows_list
      )

      return(fil_list)
    }
}


#' Filter columns and rows of `ACTIONetExperiment` or `SummarizedExperiment` object by column attribute.
#' @export
filter.ace.by.attr <- function(
  ace,
  by,
  assay_name = "counts",
  min_cells_per_feat = NULL,
  min_feats_per_cell = NULL,
  min_umis_per_cell = NULL,
  max_umis_per_cell = NULL
) {

    IDX = .get_attr_or_split_idx(ace, by)

    if (any(duplicated(rownames(ace)))) {
        msg = sprintf("Adding suffix to duplicate rownames.\n")
        warning(msg)
        rownames(ace) = make.unique(rownames(ace))
    }
    if (any(duplicated(colnames(ace)))) {
        msg = sprintf("Adding suffix to duplicate colnames.\n")
        warning(msg)
        colnames(ace) = make.unique(colnames(ace))
    }

    fil_names <- lapply(IDX, function(idx) {
        fil_list <- filter.ace(
          ace = ace[, idx],
          assay_name = assay_name,
          min_cells_per_feat = min_cells_per_feat,
          min_umis_per_cell = min_umis_per_cell,
          max_umis_per_cell = max_umis_per_cell,
          min_feats_per_cell = min_feats_per_cell,
          return_fil_ace = FALSE
        )

        return(fil_list)
    })

    fil_col = lapply(fil_names, function(i) i[["cols_filtered"]]$name)
    fil_col = Reduce(union, fil_col)


    fil_row = lapply(fil_names, function(i) i[["rows_filtered"]]$name)
    fil_row = Reduce(union, fil_row)

    keep_row = which(!(rownames(ace) %in% fil_row))
    keep_col = which(!(colnames(ace) %in% fil_col))

    ace = ace[keep_row, keep_col]
    colData(ace) <- droplevels(colData(ace))
    rowData(ace) <- droplevels(rowData(ace))

    invisible(gc())
    return(ace)
}


# #' Filter columns and rows of `ACTIONetExperiment` or `SummarizedExperiment`-like object.
# #' @export
# filter.ace <- function(
#   ace,
#   assay_name = "counts",
#   min_cells_per_feat = NULL,
#   min_feats_per_cell = NULL,
#   min_umis_per_cell = NULL,
#   max_umis_per_cell = NULL,
#   max_mito_fraction = NULL,
#   species = "mmusculus",
#   features_use = NULL,
#   return_fil_ace = TRUE
# ) {
#
#     init_dim = dim(ace)
#     init_rn = rownames(ace)
#     init_cn = colnames(ace)
#
#     if (!is.null(max_mito_fraction)){
#       mt_frac = get_mtRNA_stats(
#         ace,
#         by = NULL,
#         groups_use = NULL,
#         assay = assay_name,
#         species = species,
#         metric = "pct",
#         features_use = features_use
#       )
#
#       ace = ace[, mt_frac <= max_mito_fraction]
#     }
#
#     i = 0
#     repeat {
#         prev_dim = dim(ace)
#         rows_mask = rep(TRUE, NROW(ace))
#         cols_mask = rep(TRUE, NCOL(ace))
#         if (!is.null(min_umis_per_cell)) {
#             umi_mask = Matrix::colSums(SummarizedExperiment::assays(ace)[[assay_name]]) >= min_umis_per_cell
#             cols_mask = cols_mask & umi_mask
#         }
#
#         if (!is.null(max_umis_per_cell)) {
#             umi_mask = Matrix::colSums(SummarizedExperiment::assays(ace)[[assay_name]]) <= max_umis_per_cell
#             cols_mask = cols_mask & umi_mask
#         }
#
#         if (!is.null(min_feats_per_cell)) {
#             feature_mask = Matrix::colSums(SummarizedExperiment::assays(ace)[[assay_name]] > 0) >= min_feats_per_cell
#             cols_mask = cols_mask & feature_mask
#         }
#
#         if (!is.null(min_cells_per_feat)) {
#             if ((min_cells_per_feat < 1) & (min_cells_per_feat > 0)) {
#                 min_fc = min_cells_per_feat * init_dim[2]
#             } else {
#                 min_fc = min_cells_per_feat
#             }
#             cell_count_mask = ACTIONetExperiment:::fastRowSums(SummarizedExperiment::assays(ace)[[assay_name]] > 0) >= min_fc
#             rows_mask = rows_mask & cell_count_mask
#         }
#
#         ace <- ace[rows_mask, cols_mask]
#         invisible(gc())
#         i = i + 1
#         if (all(dim(ace) == prev_dim)) {
#             break
#         }
#     }
#     invisible(gc())
#
#     if (return_fil_ace){
#       return(ace)
#     } else {
#       fil_cols_mask = !(init_cn %in% colnames(ace))
#       fil_rows_mask = !(init_rn %in% rownames(ace))
#
#       fil_cols_list = data.frame(
#         name = init_cn[fil_cols_mask],
#         idx = which(fil_cols_mask)
#       )
#
#       fil_rows_list = data.frame(
#         name = init_rn[fil_rows_mask],
#         idx = which(fil_rows_mask)
#       )
#
#       fil_list = list(
#         cols_filtered = fil_cols_list,
#         rows_filtered = fil_rows_list
#       )
#
#       return(fil_list)
#     }
# }
