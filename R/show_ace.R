#' Shows slots of ACTIONetExperiment (ACE) object
#'
#' @param object ACTIONetExperiment object
#'
#' @export
#' @importMethodsFrom SummarizedExperiment show
setMethod("show", "ACTIONetExperiment", function(object) {
    
    validObject(object)
    
    callNextMethod()
    
    cat("rowMaps(", (length(rowMaps(object, all = F))), "): ", paste(names(rowMaps(object, 
        all = F)), collapse = " "), "\n", "colMaps(", (length(colMaps(object, all = F))), 
        "): ", paste(names(colMaps(object, all = F)), collapse = " "), "\n", "rowNets(", 
        (length(rowNets(object))), "): ", paste(names(rowNets(object)), collapse = " "), 
        "\n", "colNets(", (length(colNets(object))), "): ", paste(names(colNets(object)), 
            collapse = " "), "\n", sep = "")
    
})
