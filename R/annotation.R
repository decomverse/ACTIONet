.preprocess_annotation_markers <- function(markers, feature_set) {
  if (is.matrix(markers) || ACTIONetExperiment:::is.sparseMatrix(markers)) {
    common_features <- sort(unique(intersect(feature_set, rownames(markers))))
    row_idx <- match(common_features, rownames(markers))
    X <- markers[row_idx, ]
    rownames(X) <- common_features
  } else if (is.list(markers)) {
    if (is.null(names(markers))) {
      names(markers) <- sapply(1:length(markers), function(i) {
        msg <- sprintf("Annotation %d", i)
        print(msg)
      })
    }
    X <- do.call(cbind, lapply(markers, function(gs) {
      genes <- unlist(strsplit(gs, "[+]|[-]"))
      if (length(grep("[+|-]$", gs)) != 0) {
        x <- as.numeric(grepl("[+]", gs)) - as.numeric(grepl("[-]", gs))
        genes <- as.character(sapply(gs, function(gene) substr(gene, 1, stringr::str_length(gene) - 1)))
      } else {
        x <- rep(1, length(genes))
        genes <- gs
      }
      common.genes <- sort(unique(intersect(genes, feature_set)))
      idx1 <- match(common.genes, genes)
      idx2 <- match(common.genes, feature_set)

      v <- sparseMatrix(i = idx2, j = rep(1, length(idx2)), x = x[idx1], dims = c(length(feature_set), 1))
      return(v)
    }))
    colnames(X) <- names(markers)
    rownames(X) <- feature_set
  } else if (is.data.frame(markers) || (length(which(is(markers) == "DFrame")) != 0)) { # marker, cell type, [weight]
    if (ncol(markers) == 2) {
      markers$weight <- 1
    }
    UL <- sort(unique(markers[, 2]))
    X <- do.call(cbind, lapply(UL, function(nn) {
      idx <- which(markers[, 2] == nn)
      genes <- markers[idx, 1]
      x <- markers[idx, 3]

      common.genes <- sort(unique(intersect(genes, feature_set)))
      idx1 <- match(common.genes, genes)
      idx2 <- match(common.genes, feature_set)

      v <- sparseMatrix(i = idx2, j = rep(1, length(idx2)), x = x[idx1], dims = c(length(feature_set), 1))
      return(v)
    }))
    colnames(X) <- UL
    rownames(X) <- feature_set
  }

  return(X)
}


#' Annotate archetypes using prior cell annotations
#' (It uses t-test on the archetype footprint matrix (H))
#'
#' @param ace ACTIONet output object
#' @param labels Annotation of interest (clusters, celltypes, etc.) to test enrichment
#'
#' @return A named list: \itemize{
#' \item Label: Inferred archetype labels
#' \item Confidence: Confidence of inferred labels
#' \item Enrichment: Full enrichment matrix
#' }
#'
#' @examples
#' arch.annot <- annotate.archetypes.using.labels(ace, sce$celltypes)
#' @export
annotate.archetypes.using.labels <- function(ace,
                                             labels,
                                             archetype.slot = "H_unified", algorithm = "ttest") {
  Labels <- .preprocess_annotation_labels(labels, ace)

  if (is.matrix(ace) | ACTIONetExperiment:::is.sparseMatrix(ace)) {
    profile <- as.matrix(ace)
  } else {
    profile <- Matrix::t(colMaps(ace)[[archetype.slot]])
  }
  Annot <- names(Labels)[match(sort(unique(Labels)), Labels)]

  if (algorithm == "wilcox") {
    wilcox.out <- presto::wilcoxauc(profile, Annot[Labels])
    Enrichment <- do.call(cbind, split(-log10(wilcox.out$pval) * sign(wilcox.out$auc - 0.5), wilcox.out$group))
  } else {
    # Using t-statistics
    Enrichment <- sapply(Annot, function(label) {
      mask <- names(Labels) == label
      class.profile <- profile[, mask]
      null.profile <- profile[, !mask]

      N.class <- sum(mask)
      N.null <- sum(!mask)

      if ((N.class < 3) | (N.null < 3)) {
        return(rep(0, nrow(profile)))
      }

      mu.class <- ACTIONetExperiment:::fastRowMeans(class.profile)
      mu.null <- ACTIONetExperiment:::fastRowMeans(null.profile)

      sigma_sq.class <- apply(class.profile, 1, var)
      sigma_sq.null <- apply(null.profile, 1, var)


      delta.mean <- mu.class - mu.null
      t.stat <- delta.mean / sqrt((sigma_sq.class / N.class) + (sigma_sq.null / N.null))
      return(t.stat)
    })
  }

  Enrichment[is.na(Enrichment)] <- 0
  archetypeLabels <- Annot[apply(Enrichment, 1, which.max)]
  Labels.confidence <- apply(Enrichment, 1, max)

  rownames(Enrichment) <- paste("A", 1:nrow(Enrichment), "-", archetypeLabels,
    sep = ""
  )

  out <- list(
    Label = archetypeLabels,
    Confidence = Labels.confidence,
    Enrichment = Enrichment
  )

  return(out)
}


#' Annotate clusters using known marker genes
#' (It uses permutation test on cluster specificity scores)
#'
#' @param ace ACTIONet output object
#' @param markers A list of lists (each a set of markers for a given cell type)
#' @param rand_sample_no Number of random permutations (default=1000)
#'
#' @return A named list: \itemize{
#' \item Label: Inferred archetype labels
#' \item Confidence: Confidence of inferred labels
#' \item Enrichment: Full enrichment matrix
#' }
#'
#' @examples
#' data("curatedMarkers_human") # pre-packaged in ACTIONet
#' markers <- curatedMarkers_human$Blood$PBMC$Monaco2019.12celltypes$marker.genes
#' arch.annot <- annotate.archetypes.using.markers(ace, markers = markers)
#' @export
annotate.archetypes.using.markers <- function(ace,
                                              markers,
                                              features_use = NULL,
                                              significance_slot = "unified_feature_specificity") {
  features_use <- .get_feature_vec(ace, features_use)
  marker_mat <- .preprocess_annotation_markers(markers, features_use)

  marker_stats <- Matrix::t(assess.geneset.enrichment.from.archetypes(ace, marker_mat)$logPvals)
  colnames(marker_stats) <- colnames(marker_mat)

  marker_stats[!is.finite(marker_stats)] <- 0
  annots <- colnames(marker_mat)[apply(marker_stats, 1, which.max)]
  conf <- apply(marker_stats, 1, max)

  out <- list(
    Label = annots,
    Confidence = conf,
    Enrichment = marker_stats
  )

  return(out)
}

#' Infer cell annotations from imputed gene expression for all cells.
#'
#' @param ace ACTIONetExperiment object
#' @param markers A named list of marker genes.
#' @param features_use A vector of features of length NROW(ace) or the name of a column of rowData(ace) containing the genes given in 'markers'.
#' @param alpha_val Random-walk parameter for gene imputation.
#' @param thread_no Number of parallel threads used for gene imputation.
#' @param net_slot Name of slot in colNets(ace) containing the network to use for gene expression imputation (default="ACTIONet").
#' @param assay_name Name of assay for which to impute gene expression (default="logcounts").
#' @return A named list: \itemize{
#' \item Label: Inferred cell type labels
#' \item Confidence: Confidence of inferred labels
#' \item Enrichment: Cell type score matrix.
#' }
#'
#' @examples
#' data("curatedMarkers_human") # pre-packaged in ACTIONet
#' markers <- curatedMarkers_human$Blood$PBMC$Monaco2019.12celltypes$marker.genes
#' annots <- annotate.cells.using.markers(ace, markers = markers)
#' plot.ACTIONet(ace, annots$Label, annots$Confidence)
#' @export
annotate.cells.using.markers <- function(ace,
                                         markers,
                                         features_use = NULL,
                                         alpha_val = 0.9,
                                         thread_no = 0,
                                         net_slot = "ACTIONet",
                                         assay_name = "logcounts",
                                         max_iter = 5) {
  features_use <- .get_feature_vec(ace, features_use)
  marker_mat <- .preprocess_annotation_markers(markers, features_use)

  # marker_mat = as(sapply(marker_set, function(gs) as.numeric(features_use %in% gs) ), "sparseMatrix")
  G <- colNets(ace)[[net_slot]]
  S <- as(SummarizedExperiment::assays(ace)[[assay_name]], "sparseMatrix")

  marker_stats <- compute_marker_aggregate_stats(
    G = G,
    S = S,
    marker_mat = marker_mat,
    alpha = alpha_val,
    max_it = max_iter,
    thread_no = thread_no
  )
  colnames(marker_stats) <- colnames(marker_mat)

  marker_stats[!is.finite(marker_stats)] <- 0
  annots <- colnames(marker_mat)[apply(marker_stats, 1, which.max)]
  conf <- apply(marker_stats, 1, max)

  out <- list(
    Label = annots,
    Confidence = conf,
    Enrichment = marker_stats
  )

  return(out)
}


#' Annotates cells by interpolating marker-based archetype annotations
#'
#' @param ace ACTIONet output object
#' @param markers A list of lists (each a set of markers for a given cell type)
#' @param rand_sample_no Number of random permutations (default=1000)
#'
#' @return A named list: \itemize{
#' \item Label: Inferred archetype labels
#' \item Confidence: Confidence of inferred labels
#' \item Enrichment: Full enrichment matrix
#' }
#'
#' @examples
#'
#' data("curatedMarkers_human") # pre-packaged in ACTIONet
#' marker_set <- curatedMarkers_human$Blood$PBMC$Monaco2019.12celltypes$marker.genes
#' cell.annotations <- annotate.cells.from.archetypes.using.markers(ace, markers)
#' labels <- cell.annotations$Labels
#' @export
annotate.cells.from.archetypes.using.markers <- function(ace,
                                                         markers,
                                                         unified_suffix = "unified") {
  marker_set <- markers
  significance_slot <- sprintf("%s_feature_specificity", unified_suffix)
  arch.annot <- annotate.archetypes.using.markers(
    ace = ace,
    markers = marker_set,
    significance_slot = significance_slot
  )

  enrichment.mat <- arch.annot$Enrichment

  H.slot <- sprintf("H_%s", unified_suffix)
  cell.enrichment.mat <- map.cell.scores.from.archetype.enrichment(
    ace = ace,
    enrichment_mat = enrichment.mat,
    normalize = TRUE,
    H.slot = H.slot
  )
  cell.annotations <- colnames(cell.enrichment.mat)[apply(
    cell.enrichment.mat, 1,
    which.max
  )]

  Labels <- colnames(cell.enrichment.mat)[apply(cell.enrichment.mat, 1, which.max)]
  Labels.confidence <- apply(cell.enrichment.mat, 1, max)

  res <- list(
    Label = Labels,
    Confidence = Labels.confidence,
    Enrichment = cell.enrichment.mat
  )

  return(res)
}

#' Interpolates cell scores from archetype enrichment matrix
#'
#' @param ace ACTIONet output object
#' @param enrichment_mat Enrichment matrix with rows corresponding to archetypes and columns to an arbitrary annotation
#' @param normalize If TRUE, enrichment matrix will be first doubly-normalized
#'
#' @return Enrichment map of size cell x annotation
#'
#' @examples
#'
#' data("curatedMarkers_human") # pre-packaged in ACTIONet
#' marker_set <- curatedMarkers_human$Blood$PBMC$Monaco2019.12celltypes$marker.genes
#' arch.annot <- annotate.archetypes.using.markers(ace, markers = markers)
#' enrichment.mat <- arch.annot$enrichment
#' cell.enrichment.mat <- map.cell.scores.from.archetype.enrichment(ace, enrichment.mat)
#' cell.assignments <- colnames(cell.enrichment.mat)[apply(cell.enrichment.mat, 1, which.max)]
#' @export
map.cell.scores.from.archetype.enrichment <- function(ace,
                                                      enrichment_mat,
                                                      normalize = FALSE,
                                                      H.slot = "H_unified") {
  cell.scores.mat <- colMaps(ace)[[H.slot]]

  if (nrow(enrichment_mat) != ncol(cell.scores.mat)) {
    print("Flipping enrichment matrix")
    enrichment_mat <- Matrix::t(enrichment_mat)
  }

  if (normalize == TRUE) {
    enrichment.scaled <- doubleNorm(enrichment_mat)
  } else {
    enrichment.scaled <- enrichment_mat
    enrichment.scaled[enrichment.scaled < 0] <- 0
    if (max(enrichment.scaled) > 50) {
      enrichment.scaled <- log1p(enrichment.scaled)
    }
  }

  cell.enrichment.mat <- cell.scores.mat %*% enrichment.scaled
  colnames(cell.enrichment.mat) <- colnames(enrichment_mat)
  rownames(cell.enrichment.mat) <- colnames(ace)

  return(cell.enrichment.mat)
}


#' @export
annotateCells <- function(ace, markers, algorithm = "parametric", alpha = 0.85, network_normalization_method = "pagerank_sym", post_correction = F, thread_no = 0, features_use = NULL, assay_name = "logcounts", net_slot = "ACTIONet", specificity_slot = "unified_feature_specificity", H_slot = "H_unified") {
  if (!(net_slot %in% names(colNets(ace)))) {
    warning(sprintf("net_slot does not exist in colNets(ace)."))
    return()
  } else {
    G <- colNets(ace)[[net_slot]]
  }

  features_use <- .get_feature_vec(ace, features_use)
  marker_mat <- .preprocess_annotation_markers(markers, features_use)

  S <- assays(ace)[[assay_name]]
  S <- as(S, "sparseMatrix")

  network_normalization_code <- 0
  if (network_normalization_method == "pagerank_sym") {
    network_normalization_code <- 2
  }
  if (algorithm == "parametric") {
    out <- aggregate_genesets(G, S, marker_mat, network_normalization_code, alpha, thread_no)

    marker_stats <- out[["stats_norm_smoothed"]]
    colnames(marker_stats) <- colnames(marker_mat)
    marker_stats[!is.finite(marker_stats)] <- 0

    marker_stats_raw <- out[["stats_norm"]]
    marker_stats_raw[marker_stats_raw < 0] <- 0
    G.norm <- Matrix::t(normalize_spmat(G, 1))

    logPvals <- assess_label_enrichment(G.norm, marker_stats_raw)
    annots <- colnames(marker_mat)[apply(logPvals, 1, which.max)]
    conf <- apply(logPvals, 1, max)
  } else if (algorithm == "archetypes") {
    arch_enrichment_mat <- Matrix::t(assess.geneset.enrichment.from.archetypes(ace, marker_mat, specificity.slot = specificity_slot)$logPvals)
    arch_enrichment_mat[!is.finite(arch_enrichment_mat)] <- 0

    marker_stats <- as.matrix(map.cell.scores.from.archetype.enrichment(
      ace = ace,
      enrichment_mat = arch_enrichment_mat,
      normalize = TRUE,
      H.slot = H_slot
    ))
    marker_stats <- networkDiffusion(obj = G, scores = marker_stats, algorithm = network_normalization_method, alpha = alpha, thread_no = thread_no)

    annots <- colnames(marker_mat)[apply(marker_stats, 1, which.max)]
    conf <- apply(marker_stats, 1, max)
  } else {
    warning(sprintf("Algorithm %s not found. Reverting back to aggregate_genesets_weighted_enrichment_permutation", algorithm))
    return
  }

  if (post_correction == TRUE) {
    annots_corrected <- correct.cell.labels(ace, annots)
  } else {
    annots_corrected <- annots
  }

  out <- list(Label = annots_corrected, Confidence = conf, Enrichment = marker_stats)

  return(out)
}

#' @export
scoreCells <- function(ace, markers, algorithm = "gmm2", pre_imputation_algorithm = "none", gene_scaling_method = 0,
                       pre_alpha = 0.15, post_alpha = 0.9, network_normalization_method = "pagerank_sym", diffusion_it = 5, thread_no = 0, features_use = NULL, TFIDF_prenorm = 1, assay_name = "logcounts", net_slot = "ACTIONet", specificity_slot = "unified_feature_specificity", H_slot = "H_unified") {
  if (!(net_slot %in% names(colNets(ace)))) {
    warning(sprintf("net_slot does not exist in colNets(ace)."))
    return()
  } else {
    G <- colNets(ace)[[net_slot]]
  }

  features_use <- .get_feature_vec(ace, features_use)
  marker_mat_full <- .preprocess_annotation_markers(markers, features_use)
  mask <- Matrix::rowSums(abs(marker_mat_full)) != 0
  marker_mat <- marker_mat_full[mask, ]

  if (pre_imputation_algorithm == "none") {
    S <- assays(ace)[[assay_name]]
    sub_S <- S[mask, ]
  } else {
    sub_S <- Matrix::t(imputeGenes(ace, rownames(marker_mat), assay_name = assay_name, thread_no = thread_no, alpha = pre_alpha, diffusion_it = diffusion_it, net_slot = net_slot, algorithm = pre_imputation_algorithm))
  }
  sub_S <- as(sub_S, "sparseMatrix")

  network_normalization_code <- 0
  if (network_normalization_method == "pagerank_sym") {
    network_normalization_code <- 2
  }
  if (algorithm == "gmm2") {
    marker_stats <- aggregate_genesets_mahalanobis_2gmm(G, sub_S, marker_mat, network_normalization_method = network_normalization_code, expression_normalization_method = TFIDF_prenorm, gene_scaling_method = gene_scaling_method, pre_alpha = pre_alpha, post_alpha = post_alpha)
  } else if (algorithm == "arch2") {
    marker_stats <- aggregate_genesets_mahalanobis_2archs(G, sub_S, marker_mat, network_normalization_method = network_normalization_code, expression_normalization_method = TFIDF_prenorm, gene_scaling_method = gene_scaling_method, pre_alpha = pre_alpha, post_alpha = post_alpha)
  } else {
    warning(sprintf("Algorithm %s not found. Reverting back to gmm2", algorithm))
    marker_stats <- aggregate_genesets_mahalanobis_2gmm(G, sub_S, sub_marker_mat, network_normalization_method = network_normalization_method, expression_normalization_method = TFIDF_prenorm, gene_scaling_method = gene_scaling_method, pre_alpha = pre_alpha, post_alpha = post_alpha)
  }

  colnames(marker_stats) <- colnames(marker_mat)
  marker_stats[!is.finite(marker_stats)] <- 0
  annots <- colnames(marker_mat)[apply(marker_stats, 1, which.max)]
  conf <- apply(marker_stats, 1, max)

  out <- list(Label = annots, Confidence = conf, Enrichment = marker_stats)

  return(out)
}


annotateArchetypes <- function(ace, markers = NULL, labels = NULL, scores = NULL, archetype_slot = "H_unified", archetype_specificity_slot = "unified_feature_specificity") {
  annotations.count <- is.null(markers) + is.null(labels) + is.null(scores)
  if (annotations.count != 2) {
    stop("Exactly one of the `markers`, `labels`, or `scores` can be provided.")
  }

  if (!is.null(markers)) {
    features_use <- ACTIONet:::.get_feature_vec(ace, NULL)
    marker_mat <- as(ACTIONet:::.preprocess_annotation_markers(markers, features_use), "sparseMatrix")

    archetype_feature_specificity <- as.matrix(rowMaps(ace)[[archetype_specificity_slot]])
    colnames(archetype_feature_specificity) <- paste("A", 1:ncol(archetype_feature_specificity), sep = "")

    archetype_enrichment <- Matrix::t(assess_enrichment(archetype_feature_specificity, marker_mat)$logPvals)
    rownames(archetype_enrichment) <- colnames(archetype_feature_specificity)
    colnames(archetype_enrichment) <- colnames(marker_mat)
  } else if (!is.null(labels)) {
    X1 <- as.matrix(colMaps(ace)[[archetype_slot]])
    colnames(X1) <- paste("A", 1:ncol(X1), sep = "")

    if (length(labels) == 1) {
      l2 <- colData(ace)[[labels]]
    } else {
      l2 <- labels
    }
    f2 <- factor(l2)
    X2 <- model.matrix(~ .0 + f2)

    xi.out <- ACTIONet::XICOR(X1, X2)
    Z_pos <- xi.out$Z
    Z_pos[Z_pos < 0] <- 0
    dir <- sign(cor(X1, X2))
    archetype_enrichment <- dir * Z_pos

    rownames(archetype_enrichment) <- colnames(X1)
    colnames(archetype_enrichment) <- levels(f2)
  } else if (!is.null(scores)) {
    X1 <- as.matrix(colMaps(ace)[[archetype_slot]])
    colnames(X1) <- paste("A", 1:ncol(X1), sep = "")

    if (length(scores) == 1) {
      X2 <- as.matrix(colMaps(ace)[[scores]])
    } else {
      X2 <- as.matrix(scores)
    }

    xi.out <- ACTIONet::XICOR(X1, X2)
    Z_pos <- xi.out$Z
    Z_pos[Z_pos < 0] <- 0
    dir <- sign(cor(X1, X2))
    archetype_enrichment <- dir * Z_pos

    rownames(archetype_enrichment) <- colnames(X1)
    colnames(archetype_enrichment) <- colnames(X2)
  }
  archetype_enrichment[!is.finite(archetype_enrichment)] <- 0
  annots <- colnames(archetype_enrichment)[apply(archetype_enrichment, 1, which.max)]
  conf <- apply(archetype_enrichment, 1, max)

  out <- list(
    Label = annots,
    Confidence = conf,
    Enrichment = archetype_enrichment
  )

  return(out)
}


annotateClusters <- function(ace, markers = NULL, labels = NULL, scores = NULL, cluster_name = "leiden") {
  annotations.count <- is.null(markers) + is.null(labels) + is.null(scores)
  if (annotations.count != 2) {
    stop("Exactly one of the `markers`, `labels`, or `scores` can be provided.")
  }

  if (!is.null(markers)) {
    features_use <- .get_feature_vec(ace, NULL)
    marker_mat <- as(.preprocess_annotation_markers(markers, features_use), "sparseMatrix")

    scores <- as.matrix(rowMaps(ace)[[sprintf("%s_feature_specificity", cluster_name)]])
    cluster_enrichment <- Matrix::t(assess_enrichment(scores, marker_mat)$logPvals)
    rownames(cluster_enrichment) <- colnames(scores)
    colnames(cluster_enrichment) <- colnames(marker_mat)
  } else if (!is.null(labels)) {
    l1 <- colData(ace)[[cluster_name]]
    if (length(labels) == 1) {
      l2 <- colData(ace)[[labels]]
    } else {
      l2 <- labels
    }

    f1 <- factor(l1)
    f2 <- factor(l2)

    X1 <- model.matrix(~ .0 + f1)
    X2 <- model.matrix(~ .0 + f2)

    xi.out <- ACTIONet::XICOR(X1, X2)
    Z_pos <- xi.out$Z
    Z_pos[Z_pos < 0] <- 0
    dir <- sign(cor(X1, X2))
    cluster_enrichment <- dir * Z_pos

    rownames(cluster_enrichment) <- levels(f1)
    colnames(cluster_enrichment) <- levels(f2)
  } else if (!is.null(scores)) {
    if (length(scores) == 1) {
      X2 <- as.matrix(colMaps(ace)[[scores]])
    } else {
      X2 <- as.matrix(scores)
    }

    l1 <- colData(ace)[[cluster_name]]
    f1 <- factor(l1)
    X1 <- model.matrix(~ .0 + f1)

    xi.out <- ACTIONet::XICOR(X1, X2)
    Z_pos <- xi.out$Z
    Z_pos[Z_pos < 0] <- 0
    dir <- sign(cor(X1, X2))
    cluster_enrichment <- dir * Z_pos

    rownames(cluster_enrichment) <- levels(l1)
    colnames(cluster_enrichment) <- colnames(X2)
  }
  cluster_enrichment[!is.finite(cluster_enrichment)] <- 0
  annots <- colnames(cluster_enrichment)[apply(cluster_enrichment, 1, which.max)]
  conf <- apply(cluster_enrichment, 1, max)

  out <- list(
    Label = annots,
    Confidence = conf,
    Enrichment = cluster_enrichment
  )

  return(out)
}

projectArchs <- function(ace, archtype_scores, archetype_slot = "H_unified", normalize = TRUE) {
  cell.enrichment.mat <- map.cell.scores.from.archetype.enrichment(
    ace = ace,
    enrichment_mat = archtype_scores,
    normalize = normalize,
    H.slot = archetype_slot
  )
  cell.annotations <- colnames(cell.enrichment.mat)[apply(
    cell.enrichment.mat, 1,
    which.max
  )]

  Labels <- colnames(cell.enrichment.mat)[apply(cell.enrichment.mat, 1, which.max)]
  Labels.confidence <- apply(cell.enrichment.mat, 1, max)

  res <- list(
    Label = Labels,
    Confidence = Labels.confidence,
    Enrichment = cell.enrichment.mat
  )

  return(res)
}