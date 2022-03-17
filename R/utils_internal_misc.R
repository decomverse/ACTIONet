
.get_feature_vec <- function(ace, features_use = NULL) {
    if (is.null(features_use)) {
        features_use <- rownames(ace)
    } else {
        features_use <- ACTIONetExperiment::get.data.or.split(
          ace = ace,
          attr = features_use,
          to_return = "data",
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
        labels <- factor(ACTIONetExperiment::get.data.or.split(ace, attr = labels, to_return = "data"))
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

.validate_ace <- function(ace, allow_null = FALSE){

  if (is.null(ace) && !allow_null) {
    err <- sprintf("`ace` cannot be `NULL`.\n")
    stop(err)
  }

  if (class(ace) != "ACTIONetExperiment") {
    err <- sprintf("`ace` must be `ACTIONetExperiment`.\n")
    stop(err)
  }

  return(ace)
}


.validate_map <- function(M, map_slot = NULL, row = FALSE){
  if( is.null(M) ){
    err = sprintf("`M` must be given.\n")
    stop(err)
  }

  if (is(M, "ACTIONetExperiment")) {
    if (row == TRUE) {

      if (!(map_slot %in% names(rowMaps(M)))) {
        err <- sprintf("`%s` is not an attribute of `rowMaps`.\n", map_slot)
        stop(err)
      }
      M <- rowMaps(M)[[map_slot]]

    } else {

      if (!(map_slot %in% names(colMaps(M)))) {
        err <- sprintf("`%s` is not an attribute of `colMaps`.\n", map_slot)
        stop(err)
      }
      M <- colMaps(M)[[map_slot]]

    }
  } else {
    if (!is.matrix(M) && !ACTIONetExperiment:::is.sparseMatrix(M)) {
     err = sprintf("`M` must be `matrix` or `sparseMatrix`.\n")
     stop(err)
   }
   if (ACTIONetExperiment:::is.sparseMatrix(M) && !is(M, "dgCMatrix")) {
     M = as(M, "dgCMatrix")
   }
  }

  return(M)
}


.validate_net <- function(G, net_slot = NULL, row = FALSE){
  if( is.null(G) ){
    err = sprintf("`G` must be given.\n")
    stop(err)
  }

  if (is(G, "ACTIONetExperiment")) {
    if (row == TRUE) {

      if (!(net_slot %in% names(rowNets(G)))) {
        err <- sprintf("`%s` is not an attribute of `rowNets`.\n", net_slot)
        stop(err)
      }
      G <- rowNets(G)[[net_slot]]

    } else {

      if (!(net_slot %in% names(colNets(G)))) {
        err <- sprintf("`%s` is not an attribute of `colNets`.\n", net_slot)
        stop(err)
      }
      G <- colNets(G)[[net_slot]]

    }
  } else {
    if (!is.matrix(G) && !ACTIONetExperiment:::is.sparseMatrix(G)) {
      err = sprintf("`G` must be `matrix` or `sparseMatrix`.\n")
      stop(err)
    }
    if(!is(G, "dgCMatrix")) {
      G = as(G, "dgCMatrix")
    }
  }

  return(G)
}

.validate_attr <- function(data, attr, return_type = "data"){

  if( is.null(data) ){
    err = sprintf("`data` cannot be `NULL`.\n")
    stop(err)
  }

  if( is.null(attr) ){
    err = sprintf("`attr` cannot be `NULL`.\n")
    stop(err)
  }

  if (is(data, "ACTIONetExperiment")) {

    attr = ACTIONetExperiment::get.data.or.split(data, attr = attr, to_return = return_type)

  } else {

    attr = c(attr)
    if ( !(length(attr) == NCOL(data)) ){
      err = sprintf("`length(attr)` must equal `NCOL(data)`.\n")
      stop(err)
    }

  }

  return(attr)
}

.validate_matrix <- function(X){

  if(ACTIONetExperiment:::is.sparseMatrix(X)) {
    if(!is(X, "dgCMatrix")) {
      X = as(X, "dgCMatrix")
    }
  } else if (!is.matrix(X)) {
    err = sprintf("`X` must be `matrix` or `sparseMatrix`.\n")
    stop(err)
  }

  return(X)
}
