#' Set row-associated networks
#'
#' @return List of adjacency matrices
#'
#' @rdname rowNets
setReplaceMethod("rowNets", "ACTIONetExperiment", function(x, value) {
  value <- .check_if_mapping_list(value)
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
  value <- .check_if_mapping_list(value)
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
  value = .check_if_mapping_list(value)
	SEs = lapply(value, function(X) {
		SE = SummarizedExperiment(assays=list(X=X))
		metadata(SE)$type = "internal"
		return(SE)
	}
	x@rowMaps <- SEs

  validObject(x)
  x
})


#' Set column-associated factors
#'
#' @return List of matrices
#'
#' @rdname colMaps
setReplaceMethod("colMaps", "ACTIONetExperiment", function(x, value) {
  value <- .check_if_mapping_list(value)

  .insert_and_validate_mapping(x, value, 2)
  input_names = names(value)


  value = lapply(value, .coerce_mapping_to_SE)
  dropped_vals = sapply(value, function(v) is.null(v) | ncol(x) != ncol(x)) %>%
    names(.) == ""
  dropped_names = input_names[dropped_vals]
  value = value[!dropped]

  x = .insert_SE_to_mapping(x, value, 2)

  validObject(x)
	x
})

.check_if_mapping_list <- function(value){
  err = sprintf("New mappings must be a named list.\n")
  if(!(class(value) %in% c("list", "SimpleList")))
    stop(err)
  if(is.null(names(value)))
    stop(value)

  value = as(value, "SimpleList")
  return(value)
}

.coerce_mapping_to_SE <- function(value){
  if(class(value) == "SummarizedExperiment"){
    SE = value
  } else if(is.matrix(value) | is.sparseMatrix(value)){
    SE = SummarizedExperiment(assays=list(X=X))
  } else{
      X = as.matrix(value)
      if(is.numeric(X))
        SE = SummarizedExperiment(assays=list(X=X))
      else
        return(SummarizedExperiment())
  }

  mdata = S4Vectors::metadata(SE)
  if(!("type" %in% names(mdata))){
    mdata$type = "generic"
    S4Vectors::metadata(SE) <- mdata
  }

  return(SE)
}

.insert_SE_to_mapping <- function(x, value, dim){
  for(i in seq_along(value)){
    if(dim == 1)
      x@rowMaps[[names(value)[i]]] <- value[[i]]
    if(dim == 2)
      x@colMaps[[names(value)[i]]] <- value[[i]]
  }
  return(x)
}

.insert_and_validate_mapping <- function(x, value, 2){\
  input_names = names(value)
  value = lapply(value, .coerce_mapping_to_SE)
  dropped_vals = sapply(value, function(v) is.null(v) | ncol(x) != ncol(x)) %>%
    names(.) == ""
  dropped_names = input_names[dropped_vals]
  .dropped_vals_warning(x)
  value = value[!dropped]
  x = .insert_SE_to_mapping(x, value, 2)

}







.dropped_vals_warning <- function(x){


}



# par_func = as.character(sys.call(-1)[1])
# w = paste(sprintf("In %s: ", par_func), "Non-concatable slot <(",
#           used_slots, sprintf(")> will not be preserved.\n"), sep ="")
# warning(w, call. = FALSE)
