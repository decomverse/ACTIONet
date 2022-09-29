#' @export
normalize.scuttle <- function(ace,
                              batch_attr = NULL,
                              assay_name = "counts",
                              assay_out = "logcounts",
                              BPPARAM = SerialParam()) {
  ACTIONetExperiment:::.check_and_load_package(c("scran", "scater"))
  batch_attr <- ACTIONetExperiment::get.data.or.split(ace, attr = batch_attr, to_return = "data")

  sf <- scuttle::pooledSizeFactors(
    x = SummarizedExperiment::assays(ace)[[assay_name]],
    clusters = batch_attr,
    ref.clust = NULL,
    max.cluster.size = 3000,
    positive = TRUE,
    scaling = NULL,
    min.mean = NULL,
    subset.row = NULL,
    BPPARAM = BPPARAM
  )

  SummarizedExperiment::assays(ace)[[assay_out]] <- scuttle::normalizeCounts(
    x = SummarizedExperiment::assays(ace)[[assay_name]],
    size.factors = sf,
    log = NULL,
    transform = "log",
    pseudo.count = 1,
    center.size.factors = TRUE,
    subset.row = NULL,
    normalize.all = FALSE,
    downsample = FALSE,
    down.target = NULL,
    down.prop = 0.01,
    BPPARAM = BPPARAM,
    size_factors = NULL,
    pseudo_count = NULL,
    center_size_factors = NULL,
    subset_row = NULL,
    down_target = NULL,
    down_prop = NULL
  )

  return(ace)
}

#' @export
normalize.multiBatchNorm <- function(ace,
                                     batch_attr,
                                     assay_name = "counts",
                                     assay_out = "logcounts",
                                     BPPARAM = SerialParam()) {
  ACTIONetExperiment:::.check_and_load_package(c("scran", "batchelor"))
  batch_attr <- ACTIONetExperiment::get.data.or.split(ace, attr = batch_attr, to_return = "data")
  sce_temp <- as.SingleCellExperiment(ace)

  sce_temp <- batchelor::multiBatchNorm(
    sce_temp,
    batch = batch_attr,
    assay.type = assay_name,
    BPPARAM = BPPARAM
  )

  SummarizedExperiment::assays(ace)[[assay_out]] <- SummarizedExperiment::assays(sce_temp)[["logcounts"]]

  return(ace)
}

#' @export
normalize.Linnorm <- function(ace,
                              assay_name = "counts",
                              assay_out = "logcounts") {
  ACTIONetExperiment:::.check_and_load_package("Linnorm")
  SummarizedExperiment::assays(ace)[[assay_out]] <- Linnorm(SummarizedExperiment::assays(ace)[[assay_name]])
  return(ace)
}

#' @export
normalize.default <- function(ace,
                              assay_name = "counts",
                              assay_out = "logcounts",
                              top_features_frac = 1.0,
                              scale_param = median,
                              transformation = "log",
                              anchor_features = NULL) {
  S <- SummarizedExperiment::assays(ace)[[assay_name]]

  Snorm <- normalize.matrix(S, top_features_frac = top_features_frac, scale_param = scale_param, transformation = transformation, anchor_features = anchor_features)

  rownames(Snorm) <- rownames(ace)
  colnames(Snorm) <- colnames(ace)
  SummarizedExperiment::assays(ace)[[assay_out]] <- Snorm
  return(ace)
}

#' @export
normalize.ace <- function(ace,
                          norm_method = "default",
                          batch_attr = NULL,
                          assay_name = "counts",
                          assay_out = "logcounts",
                          top_features_frac = 1.0,
                          scale_param = median,
                          transformation = "log",
                          anchor_features = NULL,
                          BPPARAM = SerialParam(),
                          ...) {
  if (norm_method == "scran") {
    ace <- normalize.scuttle(
      ace = ace,
      batch_attr = batch_attr,
      assay_name = assay_name,
      assay_out = assay_out,
      BPPARAM = BPPARAM
    )
  } else if (norm_method == "multiBatchNorm") {
    ace <- normalize.multiBatchNorm(
      ace = ace,
      batch_attr = batch_attr,
      assay_name = assay_name,
      assay_out = assay_out,
      BPPARAM = BPPARAM
    )
  } else if (norm_method == "linnorm") {
    ace <- normalize.Linnorm(
      ace = ace,
      assay_name = assay_name,
      assay_out = assay_out
    )
  } else {
    ace <- normalize.default(
      ace = ace,
      assay_name = assay_name,
      assay_out = assay_out,
      top_features_frac = top_features_frac, scale_param = scale_param, transformation = transformation, anchor_features = anchor_features
    )
    norm_method <- sprintf("default_top%0.2f_%s", top_features_frac, transformation)
  }
  metadata(ace)$input_assay <- assay_name
  metadata(ace)$norm_method <- norm_method
  metadata(ace)$default_assay <- assay_out

  return(ace)
}


#' @export
post.normalize.ace <- function(ace, net_slot = "ACTIONet", counts_slot = "counts", normcounts_slot = "normcounts", log_transform = T, alpha_val = 0.99, lib.size = 10^6) {
  G <- colNets(ace)[[net_slot]]
  C <- assays(ace)[[counts_slot]]
  umis <- Matrix::colSums(C)
  umis.sp <- as(as.matrix(umis), "sparseMatrix")
  umis.norm <- compute_network_diffusion_approx(G, umis, alpha = alpha_val)

  x <- log(as.numeric(umis.norm[, 1]))
  x <- x - min(x, na.rm = T)
  x[is.na(x)] <- 0
  scale.factor <- lib.size * (x / max(x))

  denom <- umis
  denom[denom == 0] <- 1
  w <- scale.factor / denom

  if (is.matrix(C)) {
    B <- Matrix::t(Matrix::t(C) * w)
    if (log_transform == TRUE) {
      B <- log1p(B)
    }
  } else {
    A <- as(C, "dgTMatrix")

    x <- A@x * w[A@j + 1]
    if (log_transform == TRUE) {
      x <- log1p(x)
    }
    B <- Matrix::sparseMatrix(
      i = A@i + 1, j = A@j + 1, x = x,
      dims = dim(A)
    )
  }

  rownames(B) <- rownames(ace)
  colnames(B) <- colnames(ace)
  assays(ace)[[normcounts_slot]] <- B

  return(ace)
}