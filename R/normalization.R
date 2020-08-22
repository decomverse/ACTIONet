normalize.scran <- function(ace, BPPARAM = SerialParam()) {
    .check_and_load_package(c("scran", "scater"))
    ace = scran::computeSumFactors(ace, BPPARAM = BPPARAM)
    ace = scater::logNormCounts(ace)
    return(ace)
}


normalize.Linnorm <- function(ace) {
    .check_and_load_package("Linnorm")
    SummarizedExperiment::assays(ace)[["logcounts"]] = Linnorm(counts(ace))
    return(ace)
}

normalize.default <- function(ace, log_scale = TRUE){
  S = SummarizedExperiment::assays(ace)[["counts"]]
  B = rescale.matrix(S, log_scale)
  rownames(B) = rownames(ace)
  colnames(B) = colnames(ace)
  SummarizedExperiment::assays(ace)[["logcounts"]] = B
  return(ace)
}

rescale.matrix <- function(S, log_scale = FALSE){
  if(is.matrix(S)) {
    cs = Matrix::colSums(S)
    cs[cs == 0] = 1
    B = median(cs)*scale(S, center = F, scale = cs)
    if(log_scale == TRUE){
      B = log1p(B)
    }
  } else {
    A = as(S, "dgTMatrix")
    cs = Matrix::colSums(A)
    cs[cs == 0] = 1
    x = median(cs) * (A@x/cs[A@j + 1])
    if(log_scale == TRUE){
      x = log1p(x)
    }
    B = Matrix::sparseMatrix(i = A@i + 1, j = A@j + 1, x = x, dims = dim(A))
  }
  return(B)
}

normalize.ace <- function(ace, norm.method = "default", BPPARAM = SerialParam()) {

    if (norm.method == "scran") {
        ace = normalize.scran(ace, BPPARAM = BPPARAM)
    } else if (norm.method == "linnorm") {
        ace = normalize.Linnorm(ace)
    } else {
        ace = normalize.default(ace, log_scale = TRUE)
    }
    metadata(ace)$normalization.method = norm.method
    return(ace)
}
