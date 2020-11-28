#' @export
normalize.scran <- function(ace, batch_attr = NULL, assay = "counts", BPPARAM = SerialParam()) {
    .check_and_load_package(c("scran", "scater"))
    batch_attr = .get_ace_split_IDX(ace, batch_attr, return_split_vec = TRUE)
    ace = scran::computeSumFactors(ace, clusters = batch_attr, assay.type = assay, BPPARAM = BPPARAM)
    ace = scater::logNormCounts(ace, exprs_values = assay)
    return(ace)
}

#' @export
normalize.multiBatchNorm <- function(ace, batch_attr, assay = "counts", BPPARAM = SerialParam()) {
    .check_and_load_package(c("scran", "batchelor"))
    batch_attr = .get_ace_split_IDX(ace, batch_attr, return_split_vec = TRUE)
    ace = batchelor::multiBatchNorm(ace, batch = batch_attr, assay.type = assay, BPPARAM = BPPARAM)
    return(ace)
}

#' @export
normalize.Linnorm <- function(ace, assay = "counts") {
    .check_and_load_package("Linnorm")
    SummarizedExperiment::assays(ace)[["logcounts"]] = Linnorm(SummarizedExperiment::assays(ace)[[assay]])
    return(ace)
}

#' @export
normalize.default <- function(ace, assay = "counts", log_scale = TRUE){
  S = SummarizedExperiment::assays(ace)[[assay]]
  B = rescale.matrix(S, log_scale)
  rownames(B) = rownames(ace)
  colnames(B) = colnames(ace)
  SummarizedExperiment::assays(ace)[["logcounts"]] = B
  return(ace)
}

rescale.matrix <- function(S, log_scale = FALSE){
  if(is.matrix(S)) {
    cs = ACTIONet::fastColSums(S)
    cs[cs == 0] = 1
    B = median(cs)*scale(S, center = F, scale = cs)
    if(log_scale == TRUE){
      B = log1p(B)
    }
  } else {
    A = as(S, "dgTMatrix")
    cs = ACTIONet::fastColSums(S)
    cs[cs == 0] = 1
    x = median(cs) * (A@x/cs[A@j + 1])
    if(log_scale == TRUE){
      x = log1p(x)
    }
    B = Matrix::sparseMatrix(i = A@i + 1, j = A@j + 1, x = x, dims = dim(A))
  }
  return(B)
}

#' @export
normalize.ace <- function(ace, norm_method = "default", batch_attr = NULL, assay = "counts", BPPARAM = SerialParam()) {

    if (norm_method == "scran") {
        ace = normalize.scran(ace, batch_attr, assay = assay, BPPARAM = BPPARAM)
    } else if (norm_method == "multiBatchNorm"){
        ace = normalize.multiBatchNorm(ace, batch_attr, assay = assay, BPPARAM = BPPARAM)
    } else if (norm_method == "linnorm") {
        ace = normalize.Linnorm(ace, assay = assay)
    } else {
        ace = normalize.default(ace, assay = assay, log_scale = TRUE)
        norm_method = "default"
    }

    metadata(ace)$normalization.method = norm_method
    return(ace)
}
