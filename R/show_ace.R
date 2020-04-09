#' Shows slots of ACTIONetExperiment (ACE) object
#'
#' @param object ACTIONetExperiment object
#'
#' @examples
#' @export
#' @importMethodsFrom SingleCellExperiment show
setMethod("show", "ACTIONetExperiment", function(object) {
    callNextMethod()
    cat(		
        "rowNets(", (length(rowNets(object, withDimnames=F))), "): ", names(rowNets(object, withDimnames=F)), "\n",
        "colNets(", (length(colNets(object, withDimnames=F))), "): ", names(colNets(object, withDimnames=F)), "\n",
        "rowFactors(", (length(rowFactors(object, withDimnames=F))), "): ", names(rowFactors(object, withDimnames=F)), "\n",
        "colFactors(", (length(colFactors(object, withDimnames=F))), "): ", names(colFactors(object, withDimnames=F)), "\n",
        sep=""
    )
})
