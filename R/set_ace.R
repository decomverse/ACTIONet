#' Set row-associated networks

#' @return List of adjacency matrices
#' @rdname rowNets
#' @export
setReplaceMethod("rowNets", "ACTIONetExperiment", function(object, value) {
    value <- as(value, "SimpleList")
    object@rowNets <- value
    validObject(object)
    object
})

#' Set column-associated networks

#' @return List of adjacency matrices
#' @rdname colNets
#' @export
setReplaceMethod("colNets", "ACTIONetExperiment", function(object, value) {
    value <- as(value, "SimpleList")
    object@colNets <- value
    validObject(object)
    object
})

#' Set row-associated factors

#' @return List of matrices
#' @rdname rowMaps
#' @export
setReplaceMethod("rowMaps", "ACTIONetExperiment", function(object, value) {

     if (length(value) == 0) {
         object@rowMaps = SimpleList()
         validObject(object)
         return(object)
     }

    object <- .insert_mapping(object, value, 1)
    validObject(object)
    object
})

#' Set column-associated factors

#' @return List of matrices
#' @rdname colMaps
#' @export
setReplaceMethod("colMaps", "ACTIONetExperiment", function(object, value) {

     if (length(value) == 0) {
         object@colMaps = SimpleList()
         validObject(object)
         return(object)
     }

    object <- .insert_mapping(object, value, 2)
    validObject(object)
    object
})

#' Set row-associated factor types

#' @return List of matrices
#' @rdname rowMapTypes
#' @export
setReplaceMethod("rowMapTypes", "ACTIONetExperiment", function(object, value) {
    common_names = intersect(names(value)[sapply(value, function(object) is.character(object) &
        length(object) == 1)], names(object@rowMaps))

    for (n in common_names) {
        metadata(object@rowMaps[[n]])$type = value[[n]]
    }
    validObject(object)
    object
})

#' Set column-associated factor annotations

#' @return List of matrices
#' @rdname colMapTypes
#' @export
setReplaceMethod("colMapTypes", "ACTIONetExperiment", function(object, value) {
    common_names = intersect(names(value)[sapply(value, function(object) is.character(object) &
        length(object) == 1)], names(object@colMaps))

    for (n in common_names) {
        metadata(object@colMaps[[n]])$type = value[[n]]
    }
    validObject(object)
    object
})

#' Set column-associated factor annotations

#' @return List of matrices
#' @rdname colMapMeta
#' @export
setReplaceMethod("colMapMeta", "ACTIONetExperiment", function(object, value) {
    object <- .insert_MapMeta(object, value, 2)
    validObject(object)
    object
})

#' Set row-associated factor annotations

#' @return List of matrices
#' @rdname rowMapMeta
#' @export
setReplaceMethod("rowMapMeta", "ACTIONetExperiment", function(object, value) {
    object <- .insert_MapMeta(object, value, 1)
    validObject(object)
    object
})

#' Set column-associated reductions

#' @return List of matrices
#' @rdname colReductions
#' @export
setReplaceMethod("colReductions", "ACTIONetExperiment", function(object, value) {
    (object)
    for (i in seq_along(value)) {
        colMaps(object)[[names(value)[i]]] = value[[i]]
        colMapTypes(object)[[names(value)[i]]] = "reduction"
    }
    validObject(object)
    object
})

#' Set row-associated reductions

#' @return List of matrices
#' @rdname rowReductions
#' @export
setReplaceMethod("rowReductions", "ACTIONetExperiment", function(object, value) {
    (object)
    for (i in seq_along(value)) {
        rowMaps(object)[[names(value)[i]]] = value[[i]]
        rowMapTypes(object)[[names(value)[i]]] = "reduction"
    }
    validObject(object)
    object
})

#' Set column-associated embeddings

#' @return List of matrices
#' @rdname colEmbeddings
#' @export
setReplaceMethod("colEmbeddings", "ACTIONetExperiment", function(object, value) {
    (object)
    for (i in seq_along(value)) {
        colMaps(object)[[names(value)[i]]] = value[[i]]
        colMapTypes(object)[[names(value)[i]]] = "embedding"
    }
    validObject(object)
    object
})

#' Set row-associated embeddings

#' @return List of matrices
#' @rdname colEmbeddings
#' @export
setReplaceMethod("rowEmbeddings", "ACTIONetExperiment", function(object, value) {
    (object)
    for (i in seq_along(value)) {
        rowMaps(object)[[names(value)[i]]] = value[[i]]
        rowMapTypes(object)[[names(value)[i]]] = "embedding"
    }
    validObject(object)
    object
})

#' @export
setReplaceMethod("reducedDims", "ACTIONetExperiment", function(object, value) {
  if (length(value) == 0) {
    err = sprintf("value passed to 'reducedDims' cannot be empty. To clear column-associated reductions use 'colReductions'.\n")
    stop(err)
  }

  value = as(lapply(value, function(object) Matrix::t(object)), "SimpleList")
  value = .coerce_input_to_SE(value)
  for(i in seq_along(value)){
    value[[i]] = .set_map_type(value[[i]], "reduction", force_embed = TRUE)
  }

  object <- .insert_mapping(object, value, 2)

  validObject(object)
  object
})

#' @export
setReplaceMethod("reducedDimNames", "ACTIONetExperiment", function(object, value) {
  .validate_names(value)

  mask = colMapTypes(object) %in% c("embedding", "reduction")
  names(object@colMaps)[mask] <- value

  validObject(object)
  object
})

#' @importFrom BiocGenerics counts
#' @export
setReplaceMethod("counts", "ACTIONetExperiment", function(object, value) {
    (object)
    SummarizedExperiment::assays(object)$counts = value
    object
})

#' @export
setReplaceMethod("logcounts", "ACTIONetExperiment", function(object, value) {
    (object)
    SummarizedExperiment::assays(object)$logcounts = value
    object
})


#' @export
setReplaceMethod("normcounts", "ACTIONetExperiment", function(object, value) {
    (object)
    SummarizedExperiment::assays(object)$normcounts = value
    object
})

#' Set column-associated size factors
#' @rdname sizeFactors
#' @importFrom SummarizedExperiment colData<-
#' @importFrom BiocGenerics sizeFactors<-
#' @export
setReplaceMethod("sizeFactors", "ACTIONetExperiment", function(object, ..., value) {
    colData(object)[["sizeFactors"]] <- value
    object
})

.insert_MapMeta <- function(object, value, d){
  # value = value[names(value) != ""]
  # .check_for_duplicates(value)
  valid_names = switch(d, names(object@rowMaps), names(object@colMaps))
  .validate_names(value, valid_names)
  value = lapply(value, function(v){
    if(is.data.frame(v))
      v = DataFrame(v)
    return(v)
  })

  # maps <- switch(d, rowMapTypes(object), colMapTypes(object))

  for (n in names(value)) {
      DF = value[[n]]
      if(d == 1){
        if (nrow(DF) == nrow(object)) {
            colData(object@rowMaps[[n]]) = DF
        }
      } else if(d == 2){
        if (nrow(DF) == ncol(object)) {
            colData(object@colMaps[[n]]) = DF
        }
      }
        # if(d == 1){
        #   if (nrow(DF) == ncol(object@colMaps[[n]])) {
        #       mask = (n == names(object@colMaps))
        #       if (sum(mask) == 1) {
        #         colData(object@colMaps[[which(mask)]]) = DF
        #       }
        #   }
        # } else if(d == 2){
        #   if (nrow(DF) == ncol(object@rowMaps[[n]])) {
        #       mask = (n == names(object@rowMaps))
        #       if (sum(mask) == 1) {
        #         colData(object@rowMaps[[which(mask)]]) = DF
        #       }
        #   }
        # }
      }
  return(object)
}

.insert_mapping <- function(object, value, d){
  value = .check_if_mapping_list(value)

  map_types <- switch(d, rowMapTypes(object), colMapTypes(object))
  .validate_names(value)

  value = .coerce_input_to_SE(value)

  if (length(value) == 0) {
      value = SimpleList()
  } else{

    value <- sapply(names(value), function(n){
      v = value[[n]]
      if(dim(v)[2] != dim(object)[d]){
        err = sprintf("ncol(value) must equal %s.\n", dim(object)[d])
        stop(err)
      }
      colnames(v) <- dimnames(object)[[d]]
      if(is.null(rownames(v)))
        rownames(v) <- 1:NROW(v)

      if(is.null(S4Vectors::metadata(v)$type))
          v <- .set_map_type(v, map_types[[n]])

      return(v)
    }, simplify = FALSE)
  }

  if(d == 1){
    object@rowMaps <- as(value, "SimpleList")
  } else if(d == 2){
    object@colMaps <- as(value, "SimpleList")
  }

  return(object)
}

.validate_names <- function(value, valid_names = NULL){
  par_func = as.character(sys.call(-1)[1])
  if(any(names(value) == "")){
    err = sprintf("Values passed to '%s' cannot be unnamed.\n", par_func)
    stop(err)
  }

  if(any(duplicated(names(value)))){
    err = sprintf("Values passed to '%s' have duplicate names.\n", par_func)
    stop(err)
  }

  if(!is.null(valid_names)){
    not_in_object = setdiff(names(value), valid_names)
    if(length(not_in_object) > 0){
      err = sprintf("No element named '%s'.\n", not_in_object)
      stop(err)
    }
  }

  return
}

.coerce_input_to_SE <- function(value){
  value = value[sapply(value, function(v){!is.null(v)})]
  value = lapply(value, function(M) {
      if (any(is(M) == "SummarizedExperiment")) {
          return(M)
      } else if (is.matrix(M) | is.sparseMatrix(M)) {
          M = SummarizedExperiment(assays=list(X=M))
          return(M)
      } else {
        M = as.matrix(value)
        if(is.numeric(M)){
          M = SummarizedExperiment(assays=list(X=M))
          return(M)
        } else{
          par_func = as.character(sys.call(-1)[1])
          err = sprintf("Values passed to '%s' must be coercible to matrix, of class 'SummarizedExperiment', or NULL.\n", par_func)
          stop(err)
        }
      }
    })

  return(value)
}

.check_if_mapping_list <- function(value){
  err = sprintf("New mappings must be a named list.\n")
  if( !(class(value) %in% c("list", "SimpleList")) )
    stop(err)
  if(is.null(names(value)))
    stop(err)

  value = as(value, "SimpleList")
  return(value)
}

.set_map_type <- function(value, map_type = NULL, force_embed = FALSE){

  if(is.null(map_type))
    S4Vectors::metadata(value)$type <- ifelse(NROW(value) <= 3, "embedding", "generic")
  else{
    S4Vectors::metadata(value)$type <- map_type
    if(force_embed)
      S4Vectors::metadata(value)$type <- ifelse(NROW(value) <= 3, "embedding", map_type)
  }

  return(value)
}
