
.preprocess_annotation_features <- function(ace, features_use = NULL) {
    if (is.null(features_use)) {
        features_use <- rownames(ace)
    } else {
        features_use <- ACTIONetExperiment:::.get_attr_or_split_idx(ace,
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
              labels <- ACTIONetExperiment:::.get_attr_or_split_idx(ace, attr = labels, return_vec = TRUE)
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
