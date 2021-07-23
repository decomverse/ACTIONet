.check_and_load_package <- function(pkg_names) {
    for (pk in pkg_names) {
        if (!require(pk, character.only = T)) {
            err <- sprintf("Package '%s' is not installed.\n", pk)
            stop(err)
        }
    }
}

.get_attr_or_split_idx <- function(ace, attr, groups_use = NULL, return_vec = FALSE,
                                   d = 2) {
    if (length(attr) == 1) {
        split_vec <- switch(d,
            SummarizedExperiment::rowData(ace)[[attr]],
            SummarizedExperiment::colData(ace)[[attr]]
        )
    } else {
        if (length(attr) != dim(ace)[d]) {
            err <- sprintf("'attr' length does not match %s of ace.\n", ifelse(d ==
                1, "NROW", "NCOL"))
            stop(err)
        }
        split_vec <- attr
    }

    if (is.null(split_vec)) {
        stop(sprintf("Invalid split conditions.\n"))
    } else {
        split_vec <- as.character(split_vec)
    }

    idx <- 1:dim(ace)[d]

    if (!is.null(groups_use)) {
        sub_idx <- which(split_vec %in% groups_use)
        split_vec <- split_vec[sub_idx]
        if (is.null(split_vec)) {
              stop(sprintf("Invalid split conditions.\n"))
          }
        IDX_out <- split(idx[sub_idx], split_vec)
        return(IDX_out)
    } else {
        IDX_out <- split(idx, split_vec)
    }

    if (return_vec) {
        return(split_vec)
    } else {
        return(IDX_out)
    }
}

.preprocess_annotation_features <- function(ace, features_use = NULL) {
    if (is.null(features_use)) {
        features_use <- rownames(ace)
    } else {
        features_use <- .get_attr_or_split_idx(ace,
            attr = features_use, return_vec = TRUE,
            d = 1
        )
    }
    return(features_use)
}

.preprocess_annotation_labels <- function(labels, ace = NULL) {
    if (is.null(labels)) {
        return(NULL)
    }

    if (is.character(labels)) {
        if (length(labels) == 1) {
              labels <- .get_attr_or_split_idx(ace, attr = labels, return_vec = TRUE)
          } else {
              labels <- factor(labels)
          }
    }

    if ((length(labels) > 1) & is.logical(labels)) {
        labels <- factor(as.numeric(labels), levels = c(0, 1), labels = c("No", "Yes"))
    }

    if (is.factor(labels)) {
        v <- as.numeric(labels)
        names(v) <- levels(labels)[v]
        labels <- v
    }
    if (is.matrix(labels)) {
        L <- as.numeric(labels)
        names(L) <- names(labels)
        labels <- L
    }

    if (is.null(names(labels)) | length(unique(names(labels))) > 100) {
        names(labels) <- as.character(labels)
    }

    return(labels)
}

.tscalet <- function(A, center = TRUE, scale = TRUE) {
    A <- Matrix::t(scale(Matrix::t(A), center = center, scale = scale))
    return(A)
}

.preprocess_design_matrix_and_var_names <- function(design_mat, variable_name = NULL) {
    if (is.null(variable_name)) {
        variable_name <- colnames(design_mat)[ncol(design_mat)]
    }
    vn_idx <- which(variable_name == colnames(design_mat))[1]
    colnames(design_mat) <- make.names(colnames(design_mat), unique = TRUE, allow_ = FALSE)
    variable_name <- colnames(design_mat)[vn_idx]

    out <- list(design_mat = design_mat, variable_name = variable_name)
    return(out)
}