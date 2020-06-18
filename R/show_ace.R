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
        "rowFactors(", (length(rowFactors(object, withDimnames=F))), "): ", paste(names(rowFactors(object, withDimnames=F)), collapse = ' '), "\n",
        "colFactors(", (length(colFactors(object, withDimnames=F))), "): ", paste(names(colFactors(object, withDimnames=F)), collapse = ' '), "\n",
        sep=""
    )
})
