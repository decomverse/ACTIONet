#' Converts a SummarizedExperiment object to an ACTIONetExperiment object.
#'
#' @param from SummarizedExperiment object
#'
#' @exportMethod coerce
setAs("SummarizedExperiment", "ACTIONetExperiment", function(from) {

    ace = ACTIONetExperiment(
      assays = SummarizedExperiment::assays(from),
      rowData = from@elementMetadata,
      colData = from@colData,
      metadata = from@metadata
    )

    rownames(ace) = rownames(from)

    if (class(from) == "RangedSummarizedExperiment") {
        rowRanges(ace) = rowRanges(from)
    }

    return(ace)
})


#' Converts an ACTIONetExperiment object to a SingleCellExperiment object.
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

    return(sce)
})


#' Converts a SingleCellExperiment object to an ACTIONetExperiment object.
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

#' Converts a SummarizedExperiment, RangedSummarizedExperiment, or SingleCellExperiment object to an ACTIONetExperiment object.
#'
#' @param object SummarizedExperiment, RangedSummarizedExperiment, or SingleCellExperiment object
#'
#' @export
as.ACTIONetExperiment <- function(object) {
    if(is(object, "ACTIONetExperiment")) {
        return(object)
    }
    
    if (class(object) %in% c("SummarizedExperiment", "RangedSummarizedExperiment",
        "SingleCellExperiment")) {
        ace = as(object, "ACTIONetExperiment")
    } else {
        err = sprintf("'object' must be class 'SummarizedExperiment', 'RangedSummarizedExperiment', or 'SingleCellExperiment'.\n")
        stop(err)
    }
    return(ace)
}


#' Converts an ACTIONetExperiment method to a SummarizedExperiment object.
#'
#' @param ace ACTIONetExperiment object
#'
#' @export
as.SummarizedExperiment <- function(ace) {

    SE = SummarizedExperiment::SummarizedExperiment(
      assays = SummarizedExperiment::assays(ace),
      rowRanges = ace@rowRanges,
      colData = ace@colData,
      metadata = ace@metadata
    )

    rownames(SE) = rownames(ace)

    return(SE)
}

#' Converts an ACTIONetExperiment object to a SingleCellExperiment object.
#'
#' @param ace ACTIONetExperiment object
#'
#' @export
as.SingleCellExperiment <- function(ace) {
    ace = as(ace, "SingleCellExperiment")
    return(ace)
}
