#' An extension of the SummarizedExperiment class to store
#' the results of ACTIONet method
#'
#' @slot rowNets,colNets gene-gene and cell-cell networks, respectively
#' @slot rowMaps,colMaps Factorization results (W and H matrices)
#'
#'
#' @return an ACTIONetExperiment (ACE) object
#' @rdname ACTIONetExperiment
#' @import methods
#' @importFrom stats setNames
#' @importClassesFrom S4Vectors SimpleList
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
setClass("ACTIONetExperiment",
  slots = c(rowNets = "SimpleList",
  colNets = "SimpleList",
  rowMaps = "SimpleList",
  colMaps = "SimpleList"),
  contains = "RangedSummarizedExperiment")


#' Creates an ACTIONetExperiment (ACE) object
#'
#' @param ... SummarizedExperiment and SummarizedExperiment components
#' @param rowNets,colNets gene-gene and cell-cell networks, respectively
#' @param rowMaps,colMaps Factorization results (W and H matrices)
#'
#' @return An ACTIONetExperiment (ACE) object, derived from SummarizedExperiment, with additional slots to store ACTIONet results
#' @importFrom methods new
#' @importFrom S4Vectors SimpleList
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @export
ACTIONetExperiment <- function(...,
  rowNets = S4Vectors::SimpleList(),
  colNets = S4Vectors::SimpleList(),
  rowMaps = S4Vectors::SimpleList(),
  colMaps = S4Vectors::SimpleList()) {

  SE <- SummarizedExperiment::SummarizedExperiment(...)
  if(!is(SE, "RangedSummarizedExperiment")) {
    SE <- as(SE, "RangedSummarizedExperiment")
  }
  out = new("ACTIONetExperiment", SE, rowNets = rowNets, colNets = colNets)
  out = .insert_mapping(out, rowMaps, 1)
  out = .insert_mapping(out, colMaps, 2)


  if(NROW(out) > 0){
    # if(NCOL(SummarizedExperiment::rowData(out)) == 0){
    #   SummarizedExperiment::rowData(out) = .default_rowData(NROW(out))
    #   rownames(out) = SummarizedExperiment::rowData(out)[,1] %>% .make_chars_unique
    # }
    # if(is.null(rownames(SummarizedExperiment::rowData(out)))){
    #   rownames(out) = SummarizedExperiment::rowData(out)[,1] %>% .make_chars_unique
    # } else {
    #   rownames(out) = rownames(SummarizedExperiment::rowData(out)) %>% .make_chars_unique
    # }
    if(is.null(rownames(out))){
      rownames(out) = .default_rownames(NROW(out))
    }
  }

  if(NCOL(out) > 0){
    # if(NCOL(SummarizedExperiment::colData(out)) == 0){
    #   SummarizedExperiment::colData(out) = .default_colData(NCOL(out))
    #   colnames(out) = SummarizedExperiment::colData(out)[,1] %>% .make_chars_unique
    # }
    # if(is.null(rownames(SummarizedExperiment::colData(out)))){
    #   colnames(out) = SummarizedExperiment::colData(out)[,1] %>% .make_chars_unique
    # } else {
    #   colnames(out) = rownames(SummarizedExperiment::colData(out)) %>% .make_chars_unique
    # }
    if(is.null(colnames(out))){
      colnames(out) = .default_colnames(NCOL(out))
    }
  }

  validObject(out)
  return(out)
}


# @S3method .DollarNames ACTIONetExperiment
#' @method .DollarNames ACTIONetExperiment
#' @export
.DollarNames.ACTIONetExperiment <- function(x, pattern = "") {
    ll = c(names(colData(x)), names(rowMaps(x, all = F)), names(colMaps(x, all = F)),
        names(colNets(x)), names(rowNets(x)))
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
