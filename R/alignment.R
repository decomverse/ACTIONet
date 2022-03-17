compute.pairwise.alignment <- function(reference_profile,
                                       query_profile,
                                       reduced_dim = 50,
                                       log_transform = FALSE,
                                       deflate = TRUE,
                                       sparse = F,
                                       seed = 0) {
  set.seed(seed)

  if (log_transform == T) {
    if (max(query_profile) > 100) {
      query_profile <- log1p(query_profile)
    }

    if (max(reference_profile) > 100) {
      reference_profile <- log1p(reference_profile)
    }
  }

  query_profile <- scale(query_profile, center = F)
  reference_profile <- scale(reference_profile, center = F)

  # Deflate?
  if (deflate == TRUE) {
    query_profile.red <- reduce_kernel_full(
      S = query_profile,
      reduced_dim = reduced_dim,
      iter = 1000
    )

    reference_profile.red <- reduce_kernel_full(
      S = reference_profile,
      reduced_dim = reduced_dim,
      iter = 1000
    )

    query_profile.deflated <- query_profile + query_profile.red$A %*% Matrix::t(query_profile.red$B)
    reference_profile.deflated <- reference_profile + reference_profile.red$A %*%
      Matrix::t(reference_profile.red$B)

    query_profile_centered <- ACTIONetExperiment:::.tscalet(query_profile.deflated, scale = FALSE)
    reference_profile_centered <- ACTIONetExperiment:::.tscalet(reference_profile.deflated, scale = FALSE)
  } else {
    query_profile_centered <- ACTIONetExperiment:::.tscalet(query_profile, scale = FALSE)
    reference_profile_centered <- ACTIONetExperiment:::.tscalet(reference_profile, scale = FALSE)
  }

  reduced_dim <- min(min(min(dim(query_profile)), min(dim(reference_profile))) -
    1, reduced_dim)

  query_profile_centered[is.na(query_profile_centered)] <- 0
  reference_profile_centered[is.na(reference_profile_centered)] <- 0

  query_profile.red <- IRLB_SVD_full(
    A = Matrix::t(query_profile_centered),
    dim = reduced_dim,
    iters = 1000
  )

  reference_profile.red <- IRLB_SVD_full(
    A = Matrix::t(reference_profile_centered),
    dim = reduced_dim,
    iters = 1000
  )

  V.query <- query_profile.red$v[, 1:reduced_dim]
  V.reference <- reference_profile.red$v[, 1:reduced_dim]

  V.alignment <- V.query %*% (Matrix::t(V.query) %*% V.reference)
  S_r.query <- Matrix::t(V.alignment) %*% query_profile_centered
  S_r.reference <- Matrix::t(V.reference) %*% reference_profile_centered

  if(sparse == TRUE) { 
    X <- Matrix::t(run_simplex_regression(
      A = S_r.query,
      B = S_r.reference
    ))
  } else {
    X = cor (S_r.query, S_r.reference)
  }
  return(X)
}


plot_pairwise_alignment <- function(alignment,
                                    reference_labels = NULL,
                                    query_labels = NULL,
                                    normalize = FALSE) {
  W <- alignment - min(alignment)
  Wm <- as(MWM_hungarian(W), "dgTMatrix")

  ii <- Wm@i + 1
  jj <- Wm@j + 1
  perm <- order(Wm@x, decreasing = T)
  subW <- alignment[ii, jj]
  rownames(subW) <- reference_labels[ii]
  colnames(subW) <- query_labels[jj]

  library(seriation)
  CC <- cor(subW)
  CC[is.na(CC)] <- 0
  perm <- seriation::get_order(seriation::seriate(stats::as.dist(1 - CC), "OLO"))
  subW <- subW[perm, perm]

  gradPal <- (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu"))))(100)

  if (length(alignment) < 10000) {
    border.color <- "black"
  } else {
    border.color <- NA
  }

  X <- subW
  if (normalize == TRUE) {
    X <- doubleNorm(X)
  }

  ht <- ComplexHeatmap::Heatmap(
    matrix = X,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    rect_gp = grid::gpar(col = border.color),
    name = "Alignment score",
    col = gradPal,
    row_names_side = "left",
    row_title = "Reference",
    column_title = "Query",
    column_names_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    row_names_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    column_title_gp = grid::gpar(fontsize = 18, fontface = "bold"),
    row_title_gp = grid::gpar(fontsize = 18, fontface = "bold"),
    row_names_max_width = grid::unit(100, "cm"),
    column_names_max_height = grid::unit(100, "cm")
  )

  return(ht)
}


compute.cell.alignments.from.archetype.alignments <- function(reference_ace,
                                                              query_ace,
                                                              alignment,
                                                              alignment_threshold = 0.1,
                                                              footprint_threshold = 0.1,
                                                              reference_slot_name = "unified",
                                                              query_slot_name = "unified") {
  W <- alignment
  W[W < alignment_threshold] <- 0
  Wm <- as(MWM_hungarian(W), "dgTMatrix")

  ii <- Wm@i + 1
  jj <- Wm@j + 1
  w <- Wm@x

  reference_footprint <- colMaps(reference_ace)[[sprintf("H_%s", reference_slot_name)]]
  reference_footprint <- reference_footprint[, ii]

  query_footprint <- colMaps(query_ace)[[sprintf("H_%s", query_slot_name)]]
  query_footprint <- query_footprint[, jj]


  cell_to_cell_alignments <- lapply(1:length(ii), function(j) {
    u <- reference_footprint[, j]
    u <- u / max(u)
    v <- query_footprint[, j]
    v <- v / max(v)
    match.out <- MWM_rank1(u = u, v, footprint_threshold, footprint_threshold)

    match.ii <- match.out[1, ]
    match.jj <- match.out[2, ]
    match.w <- (w[j] * u[match.out[1, ]] * v[match.out[2, ]] * w[j])^(1 / 3)

    M <- Matrix::sparseMatrix(
      i = match.ii,
      j = match.jj,
      x = match.w,
      dims = c(length(u), length(v))
    )

    return(M)
  })

  cell_to_cell_alignment <- Reduce("+", cell_to_cell_alignments)

  cell_to_cell_alignment@x[cell_to_cell_alignment@x > 1] <- 1

  return(cell_to_cell_alignment)
}


annotate_cells_from_alignment <- function(ace,
                                          alignment,
                                          unification.slot = "H_unified") {
  cell.scores.mat <- colMaps(ace)[[unification.slot]]
  cell.enrichment.mat <- cell.scores.mat %*% Matrix::t(alignment)

  newLabels <- colnames(cell.enrichment.mat)[apply(cell.enrichment.mat, 1, which.max)]
  newLabels.confidence <- apply(cell.enrichment.mat, 1, max)

  newLabels.corrected <- correct.cell.labels(ace, newLabels)

  newLabels.annot <- annotate.archetypes.using.labels(ace, newLabels.corrected)

  M <- as(MWM_hungarian(newLabels.annot$Enrichment), "dgTMatrix")

  newLabels.CPal <- colorspace::lighten(grDevices::rgb(metadata(ace)$backbone$colors[M@i +
    1, ]), 0.25)
  names(newLabels.CPal) <- colnames(newLabels.annot$Enrichment)[M@j + 1]

  out <- list(
    Label = newLabels.corrected,
    Confidence = newLabels.confidence,
    Enrichment = cell.enrichment.mat,
    Pal = newLabels.CPal
  )

  return(out)
}