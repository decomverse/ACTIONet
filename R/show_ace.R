#' Shows slots of ACTIONetExperiment (ACE) object
#'
#' @param object ACTIONetExperiment object
#'
#' @export
#' @importMethodsFrom SummarizedExperiment show
setMethod("show", "ACTIONetExperiment", function(object) {
    callNextMethod()
    cat(		
        "rowNets(", (length(rowNets(object, withDimnames=F))), "): ", paste(names(rowNets(object, withDimnames=F)), collapse = ' '), "\n",
        "colNets(", (length(colNets(object, withDimnames=F))), "): ", paste(names(colNets(object, withDimnames=F)), collapse = ' '), "\n",
        "rowMaps(", (length(rowMaps(object, withDimnames=F))), "): ", paste(names(rowMaps(object, withDimnames=F)), collapse = ' '), "\n",
        "colMaps(", (length(colMaps(object, withDimnames=F))), "): ", paste(names(colMaps(object, withDimnames=F)), collapse = ' '), "\n",
        sep=""
    )
})
