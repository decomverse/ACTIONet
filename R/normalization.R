normalize.scran <- function(ace, BPPARAM = SerialParam()) {
    .check_and_load_package("scran", "scater")
    ace = scran::computeSumFactors(ace, BPPARAM = BPPARAM)
    ace = scater::logNormCounts(ace)
    return(ace)
}


normalize.Linnorm <- function(ace) {
    .check_and_load_package("Linnorm")
    SummarizedExperiment::assays(ace)[["logcounts"]] = Linnorm(counts(ace))
    return(ace)
}


normalize.ace <- function(ace, norm.method = "default", BPPARAM = SerialParam()) {
    
    if (norm.method == "scran") {
        ace.norm = normalize.scran(ace, BPPARAM = BPPARAM)
    } else if (norm.method == "linnorm") {
        ace.norm = normalize.Linnorm(ace)
    } else {
        ace.norm = ace
        A = as(SummarizedExperiment::assays(ace.norm)[["counts"]], "dgTMatrix")
        cs = Matrix::colSums(A)
        cs[cs == 0] = 1
        B = Matrix::sparseMatrix(i = A@i + 1, j = A@j + 1, x = log1p(median(cs) * 
            (A@x/cs[A@j + 1])), dims = dim(A))
        rownames(B) = rownames(ace.norm)
        colnames(B) = colnames(ace.norm)
        SummarizedExperiment::assays(ace.norm)[["logcounts"]] = B
    }
    
    metadata(ace.norm)$normalization.method = norm.method
    # metadata(ace.norm)$normalization.time = Sys.time()
    
    # ace.norm = add.count.metadata(ace.norm)
    
    return(ace.norm)
}
