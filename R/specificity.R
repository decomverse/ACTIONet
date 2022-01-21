#' Computes feature specificity scores for each cluster.
#' @export
clusterFeatureSpecificity <- function(
  ace = NULL,
  S = NULL,
  cluster_attr = NULL,
  output_prefix = NULL,
  assay_name = "logcounts",
  return_raw = FALSE
) {

  if( is.null(ace) && is.null(S) ){
    err = sprintf("Either 'ace' or 'S' must be given.\n")
    stop(err)
  }

  if( is.null(cluster_attr) ){
    err = sprintf("'cluster_attr' cannot be 'NULL'.\n")
    stop(err)
  }

  if( is.null(output_prefix) && return_raw == FALSE){
    err = sprintf("'output_prefix' cannot be 'NULL' if 'return_raw=FALSE'.\n")
    stop(err)
  }

  if (is.null(S)) {
    if ( !( assay_name %in% names(SummarizedExperiment::assays(ace)) ) ) {
      err <- sprintf("'%s' is not an assay of 'ace'.\n", assay_name)
      stop(err)
    }
    S <- SummarizedExperiment::assays(ace)[[assay_name]]
  } else {
    S =  as(S, "dgCMatrix")
  }

  if(!is.null(ace)){
      sa = ACTIONetExperiment::get.data.or.split(ace, attr = cluster_attr, to_return = "levels")
      clusters = sa[["index"]]
      keys = sa[["keys"]]
  } else {
    if ( length(cluster_attr) != NCOL(S) ){
      err = sprintf("'length(label_attr)' must equal 'NROW(G)'.\n")
      stop(err)
    }
    lf = factor(cluster_attr)
    clusters = as.numeric(lf)
    keys = levels(lf)
  }

  # Compute gene specificity for each cluster
  if (is.matrix(S)) {
    specificity.out <- compute_cluster_feature_specificity_full(S, clusters)
  } else {
    specificity.out <- compute_cluster_feature_specificity(S, clusters)
  }

  specificity.out <- lapply(specificity.out, function(scores) {
    colnames(scores) <- keys
    return(scores)
  })

  if (return_raw == TRUE || is.null(ace)){
    return(specificity.out)
  } else {

    specificity.out <- lapply(specificity.out, function(scores) {
      rownames(scores) <- rownames(ace)
      return(scores)
    })

    X <- specificity.out[["upper_significance"]]

    rowMaps(ace)[[sprintf("%s_feature_specificity", output_prefix)]] <- X
    rowMapTypes(ace)[[sprintf("%s_feature_specificity", output_prefix)]] <- "reduction"

    return(ace)
  }
}


#' Computes feature specificity scores for each archetype.
#' @export
archetypeFeatureSpecificity <- function(
  ace = NULL,
  S = NULL,
  H = NULL,
  assay_name = "logcounts",
  footprint_slot = "archetype_footprint",
  thread_no = 0,
  return_raw = FALSE
) {

  if( is.null(ace) && (is.null(S) || is.null(H)) ){
    err = sprintf("'S' and 'H' cannot be 'NULL' if 'ace=NULL'.\n")
    stop(err)
  }

  if (is.null(S)) {
    if (!(assay_name %in% names(assays(ace)))) {
      err <- sprintf("'S' not given and %s is not an assay of the input %s object.\n", assay_name, class(ace))
      stop(err)
    }
    S <- SummarizedExperiment::assays(ace)[[assay_name]]
  }

  if (is.null(H)) {
    if (!(footprint_slot %in% names(colMaps(ace)))) {
      err <- sprintf("'H' not given and %s is not in 'colMaps'.\n", footprint_slot)
      stop(err)
    }
    H <- Matrix::t(colMaps(ace)[[footprint_slot]])
  }

  if (is.matrix(S)) {
    specificity.out <- compute_archetype_feature_specificity_full(S, H, thread_no)
  } else {
    specificity.out <- compute_archetype_feature_specificity(S, H, thread_no)
  }

  specificity.out <- lapply(specificity.out, function(scores) {
    colnames(scores) <- paste("A", 1:ncol(scores), sep = "")
    return(scores)
  })

  if (return_raw == TRUE || is.null(ace)) {
    return(specificity.out)
  } else {

    specificity.out <- lapply(specificity.out, function(scores) {
      rownames(scores) <- rownames(ace)
      return(scores)
    })

    rowMaps(ace)[["unified_feature_profile"]] <- specificity.out[["archetypes"]]
    rowMapTypes(ace)[["unified_feature_profile"]] <- "internal"

    rowMaps(ace)[["unified_feature_specificity"]] <- specificity.out[["upper_significance"]]
    rowMapTypes(ace)[["unified_feature_specificity"]] <- "reduction"

    return(ace)
  }
}
