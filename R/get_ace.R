#' Get row-associated networks

#' @return List of adjacency matrices
#' @rdname rowNets
#' @export
setMethod("rowNets", "ACTIONetExperiment", function(object) {
    out <- object@rowNets

    out
})

#' Get column-associated networks

#' @return List of adjacency matrices
#' @rdname colNets
#' @export
setMethod("colNets", "ACTIONetExperiment", function(object) {
    out <- object@colNets

    out
})

#' Get row-associated factors

#' @return List of matrices
#' @rdname rowMaps
#' @export
setMethod("rowMaps", "ACTIONetExperiment", function(object, all = T) {
    out = as(lapply(object@rowMaps, function(M) assays(M)$X), "SimpleList")

    if (all == F & length(out) > 0) {
        mask = sapply(object@rowMaps, function(M) metadata(M)$type != "internal")
        out = out[mask]
    }

    out
})

#' Get column-associated factors

#' @return List of matrices
#' @rdname colMaps
#' @export
setMethod("colMaps", "ACTIONetExperiment", function(object, all = T) {
    out = as(lapply(object@colMaps, function(M) assays(M)$X), "SimpleList")
    if (all == F & length(out) > 0) {
        mask = sapply(object@colMaps, function(M) metadata(M)$type != "internal")
        out = out[mask]
    }
    out
})

#' Get row-associated factor types

#' @return List of types
#' @rdname rowMapTypes
#' @export
setMethod("rowMapTypes", "ACTIONetExperiment", function(object, all = T) {
    out = lapply(object@rowMaps, function(M) metadata(M)$type)
    if (all == F & length(out) > 0) {
        mask = sapply(object@rowMaps, function(M) metadata(M)$type != "internal")
        out = out[mask]
    }
    out
})

#' Get column-associated factor types

#' @return List of types
#' @rdname colMapTypes
#' @export
setMethod("colMapTypes", "ACTIONetExperiment", function(object, all = T) {
    out = lapply(object@colMaps, function(M) metadata(M)$type)
    if (all == F & length(out) > 0) {
        mask = sapply(object@colMaps, function(M) metadata(M)$type != "internal")
        out = out[mask]
    }
    out
})

#' Get row-associated factor metadata

#' @return List of types
#' @rdname rowMapTypes
#' @export
setMethod("rowMapMeta", "ACTIONetExperiment", function(object, all = T) {
    out = lapply(object@rowMaps, function(M) colData(M))
    if (all == F & length(out) > 0) {
        mask = sapply(object@rowMaps, function(M) metadata(M)$type != "internal")
        out = out[mask]
    }
    out
})

#' Get column-associated factor metadata

#' @return List of types
#' @rdname colMapTypes
#' @export
setMethod("colMapMeta", "ACTIONetExperiment", function(object, all = T) {
    out = lapply(object@colMaps, function(M) colData(M))
    if (all == F & length(out) > 0) {
        mask = sapply(object@colMaps, function(M) metadata(M)$type != "internal")
        out = out[mask]
    }
    out
})

setMethod("reducedDims", "ACTIONetExperiment", function(object) {
    Xs = colMaps(object)
    Xs = Xs[colMapTypes(object) %in% c("embedding", "reduction")]

    transposed_factors = as(lapply(Xs, function(object) Matrix::t(object)), "SimpleList")
    return(transposed_factors)
})

setMethod("reducedDimNames", "ACTIONetExperiment", function(object) {
    Xs = colMaps(object)
    Xs = Xs[colMapTypes(object) %in% c("embedding", "reduction")]

    return(names(Xs))
})


#' @export
setMethod("rowEmbeddings", "ACTIONetExperiment", function(object) {
    Xs = rowMaps(object)
    Xs = Xs[rowMapTypes(object) %in% c("embedding")]

    return(Xs)
})

#' @export
setMethod("colEmbeddings", "ACTIONetExperiment", function(object) {
    Xs = colMaps(object)
    Xs = Xs[colMapTypes(object) %in% c("embedding")]

    return(Xs)
})


#' @export
setMethod("rowReductions", "ACTIONetExperiment", function(object) {
    Maps = rowMaps(object)
    Maps = Maps[rowMapTypes(object) %in% c("reduction")]

    return(Maps)
})


#' @export
setMethod("colReductions", "ACTIONetExperiment", function(object) {
    Maps = colMaps(object)
    Maps = Maps[colMapTypes(object) %in% c("reduction")]

    return(Maps)
})

#' @rdname sizeFactors
#' @importFrom SummarizedExperiment colData
#' @importFrom BiocGenerics sizeFactors
#' @export
setMethod("sizeFactors", "ACTIONetExperiment", function(object) {
    output <- colData(object)[["sizeFactors"]]
    output
})

#' @importFrom BiocGenerics counts
#' @export
setMethod("counts", "ACTIONetExperiment", function(object) {
    (object)
    assays(object)$counts
})

#' @export
setMethod("logcounts", "ACTIONetExperiment", function(object) {
    (object)
    assays(object)$logcounts
})

#' @export
setMethod("normcounts", "ACTIONetExperiment", function(object) {
    (object)
    assays(object)$normcounts
})
