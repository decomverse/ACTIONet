#' Get row-associated networks
#'
#' @return List of adjacency matrices
#'
#' @rdname rowNets
setMethod("rowNets", "ACTIONetExperiment", function(x) {
    out <- x@rowNets
    
    out
})


#' Get column-associated networks
#'
#' @return List of adjacency matrices
#'
#' @rdname colNets
setMethod("colNets", "ACTIONetExperiment", function(x) {
    out <- x@colNets
    
    out
})


#' Get row-associated factors
#'
#' @return List of matrices
#'
#' @rdname rowMaps
setMethod("rowMaps", "ACTIONetExperiment", function(x, all = T) {
    out = as(lapply(x@rowMaps, function(M) assays(M)$X), "SimpleList")
    
    if (all == F & length(out) > 0) {
        mask = sapply(x@rowMaps, function(M) metadata(M)$type != "internal")
        out = out[mask]
    }
    
    out
})

#' Get column-associated factors
#'
#' @return List of matrices
#'
#' @rdname colMaps
setMethod("colMaps", "ACTIONetExperiment", function(x, all = T) {
    out = as(lapply(x@colMaps, function(M) assays(M)$X), "SimpleList")
    if (all == F & length(out) > 0) {
        mask = sapply(x@colMaps, function(M) metadata(M)$type != "internal")
        out = out[mask]
    }
    
    out
})

#' Get row-associated factor types
#'
#' @return List of types
#'
#' @rdname rowMapTypes
setMethod("rowMapTypes", "ACTIONetExperiment", function(x, all = T) {
    out = lapply(x@rowMaps, function(M) metadata(M)$type)
    if (all == F & length(out) > 0) {
        mask = sapply(x@rowMaps, function(M) metadata(M)$type != "internal")
        out = out[mask]
    }
    
    out
})

#' Get column-associated factor types
#'
#' @return List of types
#'
#' @rdname colMapTypes
setMethod("colMapTypes", "ACTIONetExperiment", function(x, all = T) {
    out = lapply(x@colMaps, function(M) metadata(M)$type)
    if (all == F & length(out) > 0) {
        mask = sapply(x@colMaps, function(M) metadata(M)$type != "internal")
        out = out[mask]
    }
    
    out
})



#' Get row-associated factor metadata
#'
#' @return List of types
#'
#' @rdname rowMapTypes
setMethod("rowMapMeta", "ACTIONetExperiment", function(x, all = T) {
    out = lapply(x@rowMaps, function(M) colData(M))
    if (all == F & length(out) > 0) {
        mask = sapply(x@rowMaps, function(M) metadata(M)$type != "internal")
        out = out[mask]
    }
    
    out
})

#' Get column-associated factor metadata
#'
#' @return List of types
#'
#' @rdname colMapTypes
setMethod("colMapMeta", "ACTIONetExperiment", function(x, all = T) {
    out = lapply(x@colMaps, function(M) rowData(M))
    if (all == F & length(out) > 0) {
        mask = sapply(x@colMaps, function(M) metadata(M)$type != "internal")
        out = out[mask]
    }
    
    out
})





#' Get row-associated factor metadata
#'
#' @return List of types
#'
#' @rdname rowMapTypes
setMethod("rowMapMeta", "ACTIONetExperiment", function(x, all = T) {
    out = lapply(x@rowMaps, function(M) colData(M))
    if (all == F) {
        mask = sapply(x@rowMaps, function(M) metadata(M)$type != "internal")
        out = out[mask]
    }
    
    out
})

#' Get column-associated factor metadata
#'
#' @return List of types
#'
#' @rdname colMapTypes
setMethod("colMapMeta", "ACTIONetExperiment", function(x, all = T) {
    out = lapply(x@colMaps, function(M) rowData(M))
    if (all == F) {
        mask = sapply(x@colMaps, function(M) metadata(M)$type != "internal")
        out = out[mask]
    }
    
    out
})



#' @export
#' @importFrom BiocGenerics counts
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


#' @export
setMethod("reducedDims", "ACTIONetExperiment", function(x) {
    Xs = colMaps(x)
    Xs = Xs[colMapTypes(x) %in% c("embedding", "reduction")]
    
    transposed_factors = as(lapply(Xs, function(X) Matrix::t(X)), "SimpleList")
    
    return(transposed_factors)
})


#' @export
setMethod("reducedDimNames", "ACTIONetExperiment", function(x) {
    Xs = colMaps(x)
    Xs = Xs[colMapTypes(x) %in% c("embedding", "reduction")]
    
    return(names(Xs))
})





#' @export
setMethod("rowEmbeddings", "ACTIONetExperiment", function(x) {
    Xs = rowMaps(x)
    Xs = Xs[rowMapTypes(x) %in% c("embedding")]
    
    return(Xs)
})

#' @export
setMethod("colEmbeddings", "ACTIONetExperiment", function(x) {
    Xs = colMaps(x)
    Xs = Xs[colMapTypes(x) %in% c("embedding")]
    
    return(Xs)
})


#' @export
setMethod("rowReductions", "ACTIONetExperiment", function(x) {
    Xs = rowMaps(x)
    Xs = Xs[rowMapTypes(x) %in% c("reduction")]
    
    return(Xs)
})


#' @export
setMethod("colReductions", "ACTIONetExperiment", function(x) {
    Xs = colMaps(x)
    Xs = Xs[colMapTypes(x) %in% c("reduction")]
    
    return(Xs)
})

#' @export
#' @rdname sizeFactors
#' @importFrom SummarizedExperiment colData
#' @importFrom BiocGenerics sizeFactors
setMethod("sizeFactors", "ACTIONetExperiment", function(object) {
    output <- colData(object)[["sizeFactors"]]
    output
})
