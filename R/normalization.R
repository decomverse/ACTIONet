#' @export
normalize.scran <- function(
  ace,
  batch_attr = NULL,
  assay_name = "counts",
  BPPARAM = SerialParam()
) {

    ACTIONetExperiment:::.check_and_load_package(c("scran", "scater"))
    batch_attr = ACTIONetExperiment:::.get_attr_or_split_idx(ace, batch_attr, return_vec = TRUE)

    ace = scran::computeSumFactors(
      ace,
      clusters = batch_attr,
      assay.type = assay_name,
      BPPARAM = BPPARAM
    )

    ace = scater::logNormCounts(ace, exprs_values = assay_name)

    return(ace)
}

#' @export
normalize.multiBatchNorm <- function(
  ace,
  batch_attr,
  assay_name = "counts",
  BPPARAM = SerialParam()
) {

    ACTIONetExperiment:::.check_and_load_package(c("scran", "batchelor"))
    batch_attr = ACTIONetExperiment:::.get_attr_or_split_idx(ace, batch_attr, return_vec = TRUE)
    sce_temp = as.SingleCellExperiment(ace)

    sce_temp = batchelor::multiBatchNorm(
      sce_temp,
      batch = batch_attr,
      assay.type = assay_name,
      BPPARAM = BPPARAM
    )

    SummarizedExperiment::assays(ace)[["logcounts"]] = SummarizedExperiment::assays(sce_temp)[["logcounts"]]

    return(ace)
}

#' @export
normalize.Linnorm <- function(
  ace,
  assay_name = "counts"
) {

    ACTIONetExperiment:::.check_and_load_package("Linnorm")
    SummarizedExperiment::assays(ace)[["logcounts"]] = Linnorm(SummarizedExperiment::assays(ace)[[assay_name]])
    return(ace)
}

#' @export
normalize.default <- function(
  ace,
  assay_name = "counts",
  log_scale = TRUE,
  median_scale = TRUE
) {

    S = SummarizedExperiment::assays(ace)[[assay_name]]
    B = rescale.matrix(S, log_scale, median_scale)
    rownames(B) = rownames(ace)
    colnames(B) = colnames(ace)
    SummarizedExperiment::assays(ace)[["logcounts"]] = B
    return(ace)
}

#' @export
normalize.ace <- function(
  ace,
  norm_method = "default",
  batch_attr = NULL,
  assay_name = "counts",
  BPPARAM = SerialParam()
) {

    if (norm_method == "scran") {

        ace = normalize.scran(
          ace = ace,
          batch_attr = batch_attr,
          assay_name = assay_name,
          BPPARAM = BPPARAM
        )

    } else if (norm_method == "multiBatchNorm") {

        ace = normalize.multiBatchNorm(
          ace = ace,
          batch_attr = batch_attr,
          assay_name = assay_name,
          BPPARAM = BPPARAM
        )

    } else if (norm_method == "linnorm") {

        ace = normalize.Linnorm(
          ace = ace,
          assay_name = assay_name
        )

    } else {

        ace = normalize.default(
          ace = ace,
          assay_name = assay_name,
          log_scale = TRUE,
          median_scale = TRUE
        )
        norm_method = "default"
    }

    metadata(ace)$normalization.method = norm_method

    return(ace)
}


#' @export
post.normalize.ace <- function(ace, net_slot = "ACTIONet", counts_slot = "counts", normcounts_slot = "normcounts", log_scale = T, alpha_val = 0.99, lib.size = 10^6) {
  G = colNets(ace)[[net_slot]]
  C = assays(ace)[[counts_slot]]
  umis = Matrix::colSums(C)
  umis.sp = as(as.matrix(umis), "sparseMatrix")
  umis.norm = compute_network_diffusion_fast(G, umis, alpha = alpha_val)

  x = log(as.numeric(umis.norm[, 1]))
  x = x - min(x, na.rm = T)
  x[is.na(x)] = 0
  scale.factor = lib.size * (x / max(x))

  denom = umis
  denom[denom == 0] = 1
  w = scale.factor / denom

  if(is.matrix(C)) {
    B =  Matrix::t(Matrix::t(C) * w)
    if (log_scale == TRUE) {
        B = log1p(B)
    }
  } else {
      A = as(C, "dgTMatrix")

      x = A@x*w[A@j + 1]
      if (log_scale == TRUE) {
          x = log1p(x)
      }
      B = Matrix::sparseMatrix(i = A@i + 1, j = A@j + 1, x = x,
          dims = dim(A))
  }

  rownames(B) = rownames(ace)
  colnames(B) = colnames(ace)
  assays(ace)[[normcounts_slot]] = B

  return(ace)
}
