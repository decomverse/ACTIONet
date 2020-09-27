#' Converts a SummarizedExperiment (SE) object to an ACTIONetExperiment (ACE) object
#'
#' @param from SummarizedExperiment object
#'
#' @exportMethod coerce
setAs("SummarizedExperiment", "ACTIONetExperiment", function(from) {
    # ace = ACTIONetExperiment(from, rowNets = S4Vectors::SimpleList(), colNets = S4Vectors::SimpleList(),
    #     rowMaps = S4Vectors::SimpleList(), colMaps = S4Vectors::SimpleList())

    ace = ACTIONet::ACTIONetExperiment(
      assays = SummarizedExperiment::assays(from),
      rowData = from@elementMetadata,
      colData = from@colData,
      metadata = from@metadata)

    rownames(ace) = rownames(from)

    if(class(from) == "RangedSummarizedExperiment") {
		rowRanges(ace) = rowRanges(from)
	}

    return(ace)
})

#' Converts an ACTIONetExperiment (ACE) method to a SummarizedExperiment (SE) object
#'
#' @param from ACTIONetExperiment object
#'
#' @exportMethod coerce
setAs("ACTIONetExperiment", "SummarizedExperiment", function(from) {

    SE = SummarizedExperiment::SummarizedExperiment(
      assays = SummarizedExperiment::assays(from),
      rowData = from@elementMetadata,
      colData = from@colData,
      metadata = from@metadata)

    rownames(SE) = rownames(from)

    return(SE)
})

#' Converts an ACTIONetExperiment (ACE) object to a RangedSummarizedExperiment (RSE) object
#'
#' @param from ACTIONetExperiment object
#'
#' @exportMethod coerce
setAs( "ACTIONetExperiment", "RangedSummarizedExperiment", function(from) {

    SE = SummarizedExperiment::SummarizedExperiment(
      assays = SummarizedExperiment::assays(from),
      rowRanges = from@rowRanges,
      colData = from@colData,
      metadata = from@metadata)

    rownames(SE) = rownames(from)

    return(SE)
})

#' Converts an ACTIONetExperiment (ACE) object to a SingleCellExperiment (SCE) object
#'
#' @param from ACTIONetExperiment object
#'
#' @exportMethod coerce
setAs("ACTIONetExperiment", "SingleCellExperiment", function(from) {
    SE = as(from, "RangedSummarizedExperiment")
    sce = as(SE, "SingleCellExperiment")

    Xs = colMaps(from)
    Xs = Xs[colMapTypes(from) != "internal"]

    transposed_factors = as(lapply(Xs, function(X) X), "SimpleList")
    SingleCellExperiment::reducedDims(sce) = transposed_factors

    # rowData(sce) = DataFrame(as.data.frame(rowData(from))) colData(sce) =
    # DataFrame(as.data.frame(colData(from)))

    return(sce)
})

#' Converts a SingleCellExperiment (SCE) object to an ACTIONetExperiment (ACE) object
#'
#' @param from SingleCellExperiment object
#'
#' @exportMethod coerce
setAs("SingleCellExperiment", "ACTIONetExperiment", function(from) {
    SE = as(from, "RangedSummarizedExperiment")
    rownames(SE) = rownames(from)
    rowData(SE) = rowData(from)


    ace = as(SE, "ACTIONetExperiment")

    transposed_factors = as(lapply(SingleCellExperiment::reducedDims(from), function(x) SummarizedExperiment(assays = list(X = x))), "SimpleList")
    colMaps(ace) = transposed_factors

    rowData(ace) = DataFrame(as.data.frame(rowData(ace)))
    colData(ace) = DataFrame(as.data.frame(colData(ace)))

    return(ace)
})
