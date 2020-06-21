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
	value=as(value, "SimpleList")
	value = value[names(value) != "", drop = F]
	if(length(value) > 0) {
		value = lapply(value, function(M) {
			if(length(which(is(M) == "SummarizedExperiment")) != 0) {
				if("X" %in% names(assays(M))) {
					assays(M)$X	
				}
				else {
					return(NULL)
				}
			} else if(is.matrix(M) | is.sparseMatrix(M)) {
				return(M)
			} else {
				return(NULL)
			}
		})
		value = value[sapply(value, function(SE) !is.null(SE)), drop = F]
		
		if(length(value) > 0) {
			nn = intersect(names(value), names(x@rowMaps))
			for(n in nn) {
				assays(x@rowMaps[[n]])$X = value[[n]]
			}
		
			nn = setdiff(names(value), names(x@rowMaps))
			for(n in nn) {
				SE = SummarizedExperiment(assays = list(X = value[[n]]))
				if(nrow(SE) <= 3)
					metadata(SE)$type = "embedding"
				else
					metadata(SE)$type = "generic"
				
				x@rowMaps[[n]] = SE
			}
		}
	}

  validObject(x)
	x
})



#' Set column-associated factors
#'
#' @return List of matrices
#'
#' @rdname colMaps
setReplaceMethod("colMaps", "ACTIONetExperiment", function(x, value) {
	value=as(value, "SimpleList")
	value = value[names(value) != "", drop = F]
	if(length(value) > 0) {
		value = lapply(value, function(M) {
			if(length(which(is(M) == "SummarizedExperiment")) != 0) {
				if("X" %in% names(assays(M))) {
					assays(M)$X	
				}
				else {
					return(NULL)
				}
			} else if(is.matrix(M) | is.sparseMatrix(M)) {
				return(M)
			} else {
				return(NULL)
			}
		})
		value = value[sapply(value, function(SE) !is.null(SE)), drop = F]
		
		if(length(value) > 0) {
			nn = intersect(names(value), names(x@colMaps))
			for(n in nn) {
				assays(x@colMaps[[n]])$X = value[[n]]
			}
		
			nn = setdiff(names(value), names(x@colMaps))
			for(n in nn) {
				SE = SummarizedExperiment(assays = list(X = value[[n]]))
				if(nrow(SE) <= 3)
					metadata(SE)$type = "embedding"
				else
					metadata(SE)$type = "generic"
				
				x@colMaps[[n]] = SE
			}
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
	common_names = intersect(names(value)[sapply(value, function(x) is.character(x)& length(x) == 1 )], names(x@rowMaps))

	for(n in common_names) {
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
	common_names = intersect(names(value)[sapply(value, function(x) is.character(x)& length(x) == 1 )], names(x@colMaps))
	
	for(n in common_names) {
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

	for(n in names(value)) {
		DF = value[[n]]
		if(is.data.frame(DF))
			DF = DataFrame(DF)

		if((length(which(is(DF) == "DataFrame")) != 0)) {
			if(nrow(DF) == nrow(x@colMaps[[n]])) {
				mask = (n == names(x@colMaps))
				if(sum(mask) == 1) {
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
#' @rdname rowMapTypes
setReplaceMethod("rowMapMeta", "ACTIONetExperiment", function(x, value) {
	value = value[names(value) != ""]


	for(n in names(value)) {
		DF = value[[n]]
		if(is.data.frame(DF))
			DF = DataFrame(DF)

		if((length(which(is(DF) == "DataFrame")) != 0)) {
			if(nrow(DF) == ncol(x@rowMaps[[n]])) {
				mask = (n == names(x@rowMaps))
				if(sum(mask) == 1) {
					colData(x@rowMaps[[which(mask)]]) = DF
				}
			}
		}
	}

	validObject(x)
	x
})





#' @export
setReplaceMethod("counts", "ACTIONetExperiment", function(object, value) {
	(x)	
	assays(x)$counts = value	
	x
})

#' @export
setReplaceMethod("logcounts", "ACTIONetExperiment", function(object, value) {
	(x)	
	assays(x)$logcounts = value	
	x
})


#' @export
setReplaceMethod("normcounts", "ACTIONetExperiment", function(object, value) {
	(x)	
	assays(x)$normcounts = value	
	x
})
