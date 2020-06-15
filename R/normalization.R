scran.normalize <- function(ace, BPPARAM = SerialParam()) {
    require(scran)
    require(scater)

    ace = scran::computeSumFactors(ace, BPPARAM = BPPARAM)
    ace = scater::logNormCounts(ace)
    return(ace)
}

# DESeq2.normalize <- function(ace, BPPARAM = SerialParam()) {
#     library(DESeq2)
#     require(scater)
#
#     sizeFactors(ace) <- estimateSizeFactorsForMatrix(counts(ace))
#     ace <- scater::logNormCounts(ace)
#
#     return(ace)
# }
#
# TMM.normalize <- function(ace, BPPARAM = SerialParam()) {
#     library(edgeR)
#     require(scater)
#
#     sizeFactors(ace) <- calcNormFactors(counts(ace), method = "TMM")
#     ace <- scater::logNormCounts(ace)
#
#     return(ace)
# }

# logCPM.normalize <- function(ace) {
#     library(edgeR)
#
#     logcounts(ace) = log2(edgeR::cpm(counts(ace)) + 1)
#
#     return(ace)
# }
#
linnorm.normalize <- function(ace) {
    require(Linnorm)

    SummarizedExperiment::assays(ace)[["logcounts"]] = Linnorm(counts(ace))
    return(ace)
}
#
# SCnorm.normalize <- function(ace) {
#     SCnorm_out = SCnorm(Data = counts(ace), Conditions = rep(1, ncol(ace)), FilterCellNum = 10, NCores = NUM_OF_THREAD)
#     logcounts(ace) = log2(normcounts(SCnorm_out) + 1)
#
#     return(ace)
# }
#
# scone.normalize <- function(ace) {
#     library(scone)
#
#     scaling = list(none = identity, sum = SUM_FN, tmm = TMM_FN, uq = UQ_FN, fq = FQT_FN, deseq = DESEQ_FN)
#     results = scone(SconeExperiment(counts(ace)), scaling = scaling, run = TRUE, k_qc = 0, k_ruv = 0, return.normalize = "in_memory",
#         zero = "postadjust", bpparam = BiocParallel::SerialParam())
#     out.normalize = get.normalizealized(results, method = rownames(get_params(results))[1])
#     logcounts(ace) = log2(out.normalize + 1)
#
#     return(ace)
# }

normalize.ace <- function(ace, norm.method = "default", BPPARAM = SerialParam()) {

  if (norm.method == "scran") {
      ace.norm = scran.normalize(ace,  BPPARAM = BPPARAM)
  # } else if (norm.method == "DESeq2") {
  #     ace.norm = DESeq2.normalize(ace, BPPARAM = BPPARAM)
  # } else if (norm.method == "TMM") {
  #     ace.norm = TMM.normalize(ace, BPPARAM = BPPARAM)
  #   } else if (norm.method == "logCPM") {
  #       ace.norm = logCPM.normalize(ace)
    } else if (norm.method == "linnorm") {
        ace.norm = linnorm.normalize(ace)
    # } else if (norm.method == "SCnorm") {
    #     ace.norm = SCnorm.normalize(ace)
    # } else if (norm.method == "scone") {
    #     ace.norm = scone.normalize(ace)
    } else {
        ace.norm = ace
        A = as(SummarizedExperiment::assays(ace.norm)[["counts"]], "dgTMatrix")
        cs = Matrix::colSums(A)
        cs[cs == 0] = 1
        B = Matrix::sparseMatrix(i = A@i + 1, j = A@j + 1, x = log1p(median(cs) * (A@x/cs[A@j + 1])), dims = dim(A))
        rownames(B) = rownames(ace.norm)
        colnames(B) = colnames(ace.norm)
        SummarizedExperiment::assays(ace.norm)[["logcounts"]] = B
    }

    metadata(ace.norm)$normalization.method = norm.method
    # metadata(ace.norm)$normalization.time = Sys.time()

	   ace.norm = add.count.metadata(ace.norm)

    return(ace.norm)
}

# renormalize.ace <- function(ace) {
# 	library(scater)
#
# 	system.time({ace <- computeSumFactors(ace, clusters=ace$unification.out$assignments.core)})
# 	summary(sizeFactors(ace))
#
# 	final.ace = normalize(ace)
#
# 	final.ace = add.count.metadata(final.ace)
#
#     metadata(final.ace)$normalization.method = 'renormalized'
#     metadata(final.ace)$sizeFactors = sizeFactors(final.ace)
#     metadata(final.ace)$normalization.time = Sys.time()
#
# 	return(final.ace)
# }
