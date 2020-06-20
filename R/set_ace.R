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
  value <- .check_if_mapping_list(value)
  x <- .insert_and_validate_mapping(x, value, 1)
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
  x <- .insert_and_validate_mapping(x, value, 2)
  validObject(x)
	x
})

.check_if_mapping_list <- function(value){
  err = sprintf("New mappings must be a named list.\n")
  if( !(class(value) %in% c("list", "SimpleList")) )
    stop(err)
  if(is.null(names(value)))
    stop(value)

  value = as(value, "SimpleList")
  return(value)
}

.insert_and_validate_mapping <- function(x, value, d){
  input_names = names(value)
  value = lapply(value, .coerce_mapping_to_SE)
  dropped_vals = sapply(value, function(v){
      is.null(v) | (dim(v)[d] != dim(x)[d])
  }) | names(value) == ""
  value = value[!dropped_vals]
  dropped_names = setdiff(input_names, names(values))
  if(length(dropped_names) > 0)
    .dropped_vals_warning(dropped_names)
  x = .insert_SE_to_mapping(x, value, d)
  return(x)
}

.coerce_mapping_to_SE <- function(value){
  if(class(value) == "SummarizedExperiment"){
    SE = value
  } else if(is.matrix(value) | is.sparseMatrix(value)){
    SE = SummarizedExperiment(assays=list(X=value))
  } else{
      X = as.matrix(value)
      if(is.numeric(X))
        SE = SummarizedExperiment(assays=list(X=X))
      else
        return(SummarizedExperiment())
  }

  mdata = S4Vectors::metadata(SE)
  if( !("type" %in% names(mdata)) ){
    mdata$type = "generic"
    S4Vectors::metadata(SE) <- mdata
  }
  return(SE)
}

.insert_SE_to_mapping <- function(x, value, insert_dim){
  for(i in seq_along(value)){
    if(insert_dim == 1)
      x@rowMaps[[names(value)[i]]] <- value[[i]]
    else if(insert_dim == 2)
      x@colMaps[[names(value)[i]]] <- value[[i]]
  }
  return(x)
}

.dropped_vals_warning <- function(value){
  if(length(value) == 0)
    return
  else{
    sapply(value, function(v){
      par_func = as.character(sys.call(-2)[1])
      w = sprintf("In %s: Object '%s' has incompatible format and will be dropped.\n", par_func, v)
      warning(w, call. = FALSE)
    })
  }
}
