#' An extension of the RangedSummarizedExperiment class to store
#' the results of ACTIONet method
#'
#' @slot rowNets,colNets gene-gene and cell-cell networks, respectively
#' @slot rowFactors,colFactors Factorization results (W and H matrices)
#'
#'
#' @return an ACTIONetExperiment (ACE) object
#' slot.

#' @rdname ACTIONetExperiment
#' @export
#' @import methods
#' @importFrom stats setNames
#' @importClassesFrom S4Vectors SimpleList
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
.ACTIONetExperiment <- setClass("ACTIONetExperiment",
		slots= representation(rowNets = "SimpleList", colNets = "SimpleList", rowFactors = "SimpleList", colFactors = "SimpleList"),
		contains = "RangedSummarizedExperiment"        
)
         
         
#' Creates an ACTIONetExperiment (ACE) object
#'
#' @param ... RangedSummarizedExperiment and SummarizedExperiment components
#' @param rowNets,colNets gene-gene and cell-cell networks, respectively
#' @param rowFactors,colFactors Factorization results (W and H matrices)
#'
#' @return An ACTIONetExperiment (ACE) object, derived from RangedSummarizedExperiment, with additional slots to store ACTIONet results
#'
#' @export
#' @import RangedSummarizedExperiment SummarizedExperiment
ACTIONetExperiment <- function(rowNets=S4Vectors::SimpleList(), 
    colNets=S4Vectors::SimpleList(), 
    rowFactors=S4Vectors::SimpleList(), 
    colFactors=S4Vectors::SimpleList(),
    ...)
{
	SE <- SummarizedExperiment::RangedSummarizedExperiment(...)
	out <- .ACTIONetExperiment(SE, rowNets=rowNets, colNets=colNets, rowFactors=rowFactors, colFactors=colFactors)
	return(out)
}





#' @S3method .DollarNames ACTIONetExperiment
.DollarNames.ACTIONetExperiment <- function(x, pattern = "") {
	ll = c( names(colData(x)), names(rowFactors(x)), names(colFactors(x)), names(colNets(x)), names(rowNets(x)))
    grep(pattern, ll, value=TRUE)
}

#' @export
setMethod('.DollarNames', 'ACTIONetExperiment', .DollarNames.ACTIONetExperiment)



#' @export
setMethod("$", "ACTIONetExperiment",
    function(x, name)
{
	if(name %in% names(colData(x))) {
		colData(x)[[name]]
	} else if (name %in% names(rowFactors(x))) {
		rowFactors(x)[[name]]
	} else if (name %in% names(colFactors(x))) {
		colFactors(x)[[name]]
	} else if (name %in% names(colNets(x))) {
		colNets(x)[[name]]
	} else if(name %in% names(rowNets(x))) {
		rowNets(x)[[name]]
	} else {
		message(sprintf("Attribute %s not found", name))
	}
	
})

#' @export
setReplaceMethod("$", "ACTIONetExperiment",
    function(x, name, value)
{
    colData(x)[[name]] <- value
    x
})
