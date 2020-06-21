#' Set row-associated networks
#'
#' @return List of adjacency matrices
#'
#' @rdname rowNets
setReplaceMethod("rowNets", "ACTIONetExperiment", function(x, value) {
    value <- as(value, "SimpleList")
    x@rowNets <- value
    validObject(x)
    x
})

#' Set column-associated networks
#'
#' @return List of adjacency matrices
#'
#' @rdname colNets
setReplaceMethod("colNets", "ACTIONetExperiment", function(x, value) {
    value <- as(value, "SimpleList")
    x@colNets <- value
    validObject(x)
    x
})


#' Set row-associated factors
#'
#' @return List of matrices
#'
#' @rdname rowMaps
setReplaceMethod("rowMaps", "ACTIONetExperiment", function(x, value) {

    if(length(value) == 0){
      x@rowMaps = SimpleList()
      validObject(x)
      return(x)
    }

    value = as(value, "SimpleList")
    value = value[names(value) != "", drop = F]

    if (length(value) > 0) {
        value = lapply(value, function(M) {
            if (length(which(is(M) == "SummarizedExperiment")) != 0) {
                if ("X" %in% names(assays(M))) {
                  assays(M)$X
                } else {
                  return(NULL)
                }
            } else if (is.matrix(M) | is.sparseMatrix(M)) {
                return(M)
            } else {
                return(NULL)
            }
        })
        value = value[sapply(value, function(SE) !is.null(SE)), drop = F]

        if (length(value) > 0) {
            nn = intersect(names(value), names(x@rowMaps))
			mismatched_dims = names(which(!sapply(nn, function(n) all(dim(rowMaps(x)[[n]]) == dim(value[[n]])))))			
            nn = setdiff(nn, mismatched_dims)
            
            
            for (n in nn) {
                assays(x@rowMaps[[n]])$X = value[[n]]
            }

            nn = union(setdiff(names(value), names(x@rowMaps)), mismatched_dims)
            for (n in nn) {
                X = value[[n]]
                if (is.null(colnames(X)))
                  colnames(X) = 1:ncol(X)
                rownames(X) = rownames(x)

                SE = SummarizedExperiment(assays = list(X = X))
                if (nrow(SE) <= 3)
                  metadata(SE)$type = "embedding" else metadata(SE)$type = "generic"

                x@rowMaps[[n]] = SE
            }

            nn = setdiff(names(x@rowMaps), names(value))
            x@rowMaps = x@rowMaps[! (names(x@rowMaps) %in% nn) ]
        }
    }

    validObject(x)
    return(x)
})



#' Set column-associated factors
#'
#' @return List of matrices
#'
#' @rdname colMaps
setReplaceMethod("colMaps", "ACTIONetExperiment", function(x, value) {

    if(length(value) == 0){
      x@colMaps = SimpleList()
      validObject(x)
      return(x)
    }

    value = as(value, "SimpleList")
    value = value[names(value) != "", drop = F]
    if (length(value) > 0) {
        value = lapply(value, function(M) {
            if (length(which(is(M) == "SummarizedExperiment")) != 0) {
                if ("X" %in% names(assays(M))) {
                  assays(M)$X
                } else {
                  return(NULL)
                }
            } else if (is.matrix(M) | is.sparseMatrix(M)) {
                return(M)
            } else {
                return(NULL)
            }
        })
        value = value[sapply(value, function(SE) !is.null(SE)), drop = F]

        if (length(value) > 0) {
            nn = intersect(names(value), names(x@colMaps))
			mismatched_dims = names(which(!sapply(nn, function(n) all(dim(colMaps(x)[[n]]) == dim(value[[n]])))))			
            nn = setdiff(nn, mismatched_dims)
            
            for (n in nn) {
                assays(x@colMaps[[n]])$X = value[[n]]
            }

            nn = union(setdiff(names(value), names(x@colMaps)), mismatched_dims)
            for (n in nn) {
                X = value[[n]]
                if (is.null(rownames(X)))
                  rownames(X) = 1:nrow(X)
                colnames(X) = colnames(x)

                SE = SummarizedExperiment(assays = list(X = X))
                if (nrow(SE) <= 3)
                  metadata(SE)$type = "embedding" else metadata(SE)$type = "generic"

                x@colMaps[[n]] = SE
            }


            nn = setdiff(names(x@colMaps), names(value))
            x@colMaps = x@colMaps[! (names(x@colMaps) %in% nn) ]

        }
    }

    validObject(x)
    x
})




#' Set row-associated factor types
#'
#' @return List of matrices
#'
#' @rdname rowMapTypes
setReplaceMethod("rowMapTypes", "ACTIONetExperiment", function(x, value) {
    common_names = intersect(names(value)[sapply(value, function(x) is.character(x) & length(x) == 1)], names(x@rowMaps))

    for (n in common_names) {
        metadata(x@rowMaps[[n]])$type = value[[n]]
    }
    validObject(x)
    x
})


#' Set column-associated factor annotations
#'
#' @return List of matrices
#'
#' @rdname colMapTypes
setReplaceMethod("colMapTypes", "ACTIONetExperiment", function(x, value) {
    common_names = intersect(names(value)[sapply(value, function(x) is.character(x) & length(x) == 1)], names(x@colMaps))

    for (n in common_names) {
        metadata(x@colMaps[[n]])$type = value[[n]]
    }
    validObject(x)
    x
})


#' Set column-associated factor annotations
#'
#' @return List of matrices
#'
#' @rdname colMapMeta
setReplaceMethod("colMapMeta", "ACTIONetExperiment", function(x, value) {
    value = value[names(value) != ""]

    for (n in names(value)) {
        DF = value[[n]]
        if (is.data.frame(DF))
            DF = DataFrame(DF)

        if ((length(which(is(DF) == "DataFrame")) != 0)) {
            if (nrow(DF) == nrow(x@colMaps[[n]])) {
                mask = (n == names(x@colMaps))
                if (sum(mask) == 1) {
                  rowData(x@colMaps[[which(mask)]]) = DF
                }
            }
        }
    }

    validObject(x)
    x
})



#' Set row-associated factor annotations
#'
#' @return List of matrices
#'
#' @rdname rowMapMeta
setReplaceMethod("rowMapMeta", "ACTIONetExperiment", function(x, value) {
    value = value[names(value) != ""]


    for (n in names(value)) {
        DF = value[[n]]
        if (is.data.frame(DF))
            DF = DataFrame(DF)

        if ((length(which(is(DF) == "DataFrame")) != 0)) {
            if (nrow(DF) == ncol(x@rowMaps[[n]])) {
                mask = (n == names(x@rowMaps))
                if (sum(mask) == 1) {
                  colData(x@rowMaps[[which(mask)]]) = DF
                }
            }
        }
    }

    validObject(x)
    x
})

#' Set column-associated reductions
#'
#' @return List of matrices
#'
#' @rdname colReductions
setReplaceMethod("colReductions", "ACTIONetExperiment", function(x, value) {
    (x)
    for(i in seq_along(value)){
      colMaps(x)[[names(value)[i]]] = value[[i]]
      colMapTypes(x)[[names(value)[i]]] = "reduction"
    }
    x
})

#' Set row-associated reductions
#'
#' @return List of matrices
#'
#' @rdname rowReductions
setReplaceMethod("rowReductions", "ACTIONetExperiment", function(x, value) {
    (x)
    for(i in seq_along(value)){
      rowMaps(x)[[names(value)[i]]] = value[[i]]
      rowMapTypes(x)[[names(value)[i]]] = "reduction"
    }
    x
})

#' Set column-associated embeddings
#'
#' @return List of matrices
#'
#' @rdname colEmbeddings
setReplaceMethod("colEmbeddings", "ACTIONetExperiment", function(x, value) {
    (x)
    for(i in seq_along(value)){
      colMaps(x)[[names(value)[i]]] = value[[i]]
      colMapTypes(x)[[names(value)[i]]] = "embedding"
    }
    x
})

#' Set row-associated embeddings
#'
#' @return List of matrices
#'
#' @rdname colEmbeddings
setReplaceMethod("rowEmbeddings", "ACTIONetExperiment", function(x, value) {
    (x)
    for(i in seq_along(value)){
      rowMaps(x)[[names(value)[i]]] = value[[i]]
      rowMapTypes(x)[[names(value)[i]]] = "embedding"
    }
    x
})

#' @export
setReplaceMethod("counts", "ACTIONetExperiment", function(object, value) {
    (object)
    assays(object)$counts = value
    object
})

#' @export
setReplaceMethod("logcounts", "ACTIONetExperiment", function(object, value) {
    (object)
    assays(object)$logcounts = value
    object
})


#' @export
setReplaceMethod("normcounts", "ACTIONetExperiment", function(object, value) {
    (object)
    assays(object)$normcounts = value
    object
})
