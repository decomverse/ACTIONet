#' An extension of the SummarizedExperiment class to store
#' the results of ACTIONet method
#'
#' @slot rowNets,colNets gene-gene and cell-cell networks, respectively
#' @slot rowMaps,colMaps Factorization results (W and H matrices)
#'
#'
#' @return an ACTIONetExperiment (ACE) object
#' slot.

#' @rdname ACTIONetExperiment
#' @export
#' @import methods
#' @importFrom stats setNames
#' @importClassesFrom S4Vectors SimpleList
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.ACTIONetExperiment <- setClass("ACTIONetExperiment", slots = representation(rowNets = "SimpleList", colNets = "SimpleList", rowMaps = "SimpleList", 
    colMaps = "SimpleList"), contains = "SummarizedExperiment")


#' Creates an ACTIONetExperiment (ACE) object
#'
#' @param ... SummarizedExperiment and SummarizedExperiment components
#' @param rowNets,colNets gene-gene and cell-cell networks, respectively
#' @param rowMaps,colMaps Factorization results (W and H matrices)
#'
#' @return An ACTIONetExperiment (ACE) object, derived from SummarizedExperiment, with additional slots to store ACTIONet results
#'
#' @export
#' @import SummarizedExperiment SummarizedExperiment
ACTIONetExperiment <- function(rowNets = S4Vectors::SimpleList(), colNets = S4Vectors::SimpleList(), rowMaps = S4Vectors::SimpleList(), colMaps = S4Vectors::SimpleList(), 
    ...) {
    SE <- SummarizedExperiment::SummarizedExperiment(...)
    out <- .ACTIONetExperiment(SE, rowNets = rowNets, colNets = colNets, rowMaps = rowMaps, colMaps = colMaps)
    return(out)
}




#' @export
#' @S3method .DollarNames ACTIONetExperiment
.DollarNames.ACTIONetExperiment <- function(x, pattern = "") {
    ll = c(names(colData(x)), names(rowMaps(x, all = F)), names(colMaps(x, all = F)), names(colNets(x)), names(rowNets(x)))
    grep(pattern, ll, value = TRUE)
}

#' @export
setMethod(".DollarNames", "ACTIONetExperiment", .DollarNames.ACTIONetExperiment)



#' @export
setMethod("$", "ACTIONetExperiment", function(x, name) {
    if (name %in% names(colData(x))) {
        colData(x)[[name]]
    } else if (name %in% names(rowMaps(x, all = F))) {
        rowMaps(x)[[name]]
    } else if (name %in% names(colMaps(x, all = F))) {
        colMaps(x)[[name]]
    } else if (name %in% names(colNets(x))) {
        colNets(x)[[name]]
    } else if (name %in% names(rowNets(x))) {
        rowNets(x)[[name]]
    } else {
        message(sprintf("Attribute %s not found", name))
    }
    
})

#' @export
setReplaceMethod("$", "ACTIONetExperiment", function(x, name, value) {
    if (name %in% names(colData(x))) {
        colData(x)[[name]] <- value
    } else if (name %in% names(rowMaps(x))) {
        rowMaps(x)[[name]] <- value
    } else if (name %in% names(colMaps(x))) {
        colMaps(x)[[name]] <- value
    } else if (name %in% names(colNets(x))) {
        colNets(x)[[name]] <- value
    } else if (name %in% names(rowNets(x))) {
        rowNets(x)[[name]] <- value
    } else {
        colData(x)[[name]] <- value
    }
    
    x
})
