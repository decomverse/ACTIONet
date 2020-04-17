#' An extension of the SingleCellExperiment class to store
#' the results of ACTIONet method
#'
#' @slot rowNets,colNets gene-gene and cell-cell networks, respectively
#' @slot rowFactors,colFactors Factorization results (W and H matrices)
#'
#'
#' @return an ACTIONetExperiment (ACE) object
#' slot.

#' @rdname ACTIONetExperiment
#' @import methods
#' @importFrom stats setNames
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom S4Vectors SimpleList
#' @export
.ACTIONetExperiment <- setClass("ACTIONetExperiment",
		slots= representation(rowNets = "SimpleList", colNets = "SimpleList", rowFactors = "SimpleList", colFactors = "SimpleList"),
		contains = "SingleCellExperiment"        
)
         
         
#' Creates an ACTIONetExperiment (ACE) object
#'
#' @param ... SingleCellExperiment and SummarizedExperiment components
#' @param rowNets,colNets gene-gene and cell-cell networks, respectively
#' @param rowFactors,colFactors Factorization results (W and H matrices)
#'
#' @return An ACTIONetExperiment (ACE) object, derived from SingleCellExperiment, with additional slots to store ACTIONet results
#'
#' @import SingleCellExperiment SummarizedExperiment
#'
#' @export
ACTIONetExperiment <- function(...,
    rowNets=S4Vectors::SimpleList(), 
    colNets=S4Vectors::SimpleList(), 
    rowFactors=S4Vectors::SimpleList(), 
    colFactors=S4Vectors::SimpleList())
{
	sce <- SingleCellExperiment::SingleCellExperiment(...)
	out <- .ACTIONetExperiment(sce, rowNets=rowNets, colNets=colNets, rowFactors=rowFactors, colFactors=colFactors)
	return(out)
}


