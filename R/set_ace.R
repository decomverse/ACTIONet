#' Set row-associated networks

#' @return List of adjacency matrices
#' @rdname rowNets
#' @export
setReplaceMethod("rowNets", "ACTIONetExperiment", function(x, value) {
    value <- as(value, "SimpleList")
    x@rowNets <- value
    validObject(x)
    x
})

#' Set column-associated networks

#' @return List of adjacency matrices
#' @rdname colNets
#' @export
setReplaceMethod("colNets", "ACTIONetExperiment", function(x, value) {
    value <- as(value, "SimpleList")
    x@colNets <- value
    validObject(x)
    x
})

#' Set row-associated factors

#' @return List of matrices
#' @rdname rowMaps
#' @export
setReplaceMethod("rowMaps", "ACTIONetExperiment", function(x, value) {

     if (length(value) == 0) {
         x@rowMaps = SimpleList()
         validObject(x)
         return(x)
     }

    x <- .insert_mapping(x, value, 1)
    validObject(x)
    x
})

#' Set column-associated factors

#' @return List of matrices
#' @rdname colMaps
#' @export
setReplaceMethod("colMaps", "ACTIONetExperiment", function(x, value) {

     if (length(value) == 0) {
         x@colMaps = SimpleList()
         validObject(x)
         return(x)
     }

    x <- .insert_mapping(x, value, 2)
    validObject(x)
    x
})

#' Set row-associated factor types

#' @return List of matrices
#' @rdname rowMapTypes
#' @export
setReplaceMethod("rowMapTypes", "ACTIONetExperiment", function(x, value) {
    common_names = intersect(names(value)[sapply(value, function(x) is.character(x) &
        length(x) == 1)], names(x@rowMaps))

    for (n in common_names) {
        metadata(x@rowMaps[[n]])$type = value[[n]]
    }
    validObject(x)
    x
})

#' Set column-associated factor annotations

#' @return List of matrices
#' @rdname colMapTypes
#' @export
setReplaceMethod("colMapTypes", "ACTIONetExperiment", function(x, value) {
    common_names = intersect(names(value)[sapply(value, function(x) is.character(x) &
        length(x) == 1)], names(x@colMaps))

    for (n in common_names) {
        metadata(x@colMaps[[n]])$type = value[[n]]
    }
    validObject(x)
    x
})

#' Set column-associated factor annotations

#' @return List of matrices
#' @rdname colMapMeta
#' @export
setReplaceMethod("colMapMeta", "ACTIONetExperiment", function(obj, val) {
    obj <- .insert_MapMeta(obj, val, 2)
    validObject(obj)
    obj
})

#' Set row-associated factor annotations

#' @return List of matrices
#' @rdname rowMapMeta
#' @export
setReplaceMethod("rowMapMeta", "ACTIONetExperiment", function(obj, val) {
    obj <- .insert_MapMeta(obj, val, 1)
    validObject(obj)
    obj
})

#' Set column-associated reductions

#' @return List of matrices
#' @rdname colReductions
#' @export
setReplaceMethod("colReductions", "ACTIONetExperiment", function(obj, val) {
    (obj)
    for (i in seq_along(val)) {
        colMaps(obj)[[names(val)[i]]] = val[[i]]
        colMapTypes(obj)[[names(val)[i]]] = "reduction"
    }
    obj
})

#' Set row-associated reductions

#' @return List of matrices
#' @rdname rowReductions
#' @export
setReplaceMethod("rowReductions", "ACTIONetExperiment", function(obj, val) {
    (obj)
    for (i in seq_along(val)) {
        rowMaps(obj)[[names(val)[i]]] = val[[i]]
        rowMapTypes(obj)[[names(val)[i]]] = "reduction"
    }
    obj
})

#' Set column-associated embeddings

#' @return List of matrices
#' @rdname colEmbeddings
#' @export
setReplaceMethod("colEmbeddings", "ACTIONetExperiment", function(x, value) {
    (x)
    for (i in seq_along(value)) {
        colMaps(x)[[names(value)[i]]] = value[[i]]
        colMapTypes(x)[[names(value)[i]]] = "embedding"
    }
    x
})

#' Set row-associated embeddings

#' @return List of matrices
#' @rdname colEmbeddings
#' @export
setReplaceMethod("rowEmbeddings", "ACTIONetExperiment", function(x, value) {
    (x)
    for (i in seq_along(value)) {
        rowMaps(x)[[names(value)[i]]] = value[[i]]
        rowMapTypes(x)[[names(value)[i]]] = "embedding"
    }
    x
})

#' @export
setReplaceMethod("reducedDims", "ACTIONetExperiment", function(x, value) {
  if (length(value) == 0) {
    err = sprintf("value passed to 'reducedDims' cannot be empty. To clear column-associated reductions use 'colReductions'.\n")
    stop(err)
  }

  value = as(lapply(value, function(X) Matrix::t(X)), "SimpleList")
  value = .coerce_input_to_SE(value)
  value = .set_map_type(value, "reduction", force_embed = TRUE)

  x <- .insert_mapping(x, value, 2)

  validObject(x)
  x
})

#' @export
setReplaceMethod("reducedDimNames", "ACTIONetExperiment", function(x, value) {
  .validate_names(value)

  mask = colMapTypes(x) %in% c("embedding", "reduction")
  names(colMaps(x)[mask]) <- value

  validObject(x)
  x
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

.insert_MapMeta <- function(obj, val, d){
  # val = val[names(val) != ""]
  # .check_for_duplicates(val)
  .validate_names(val)

  for (n in names(val)) {
      DF = val[[n]]
      if (is.data.frame(DF))
          DF = DataFrame(DF)

      if (any(is(DF) == "DataFrame")) {

        if(d == 1){
          if (nrow(DF) == ncol(obj@colMaps[[n]])) {
              mask = (n == names(obj@colMaps))
              if (sum(mask) == 1) {
                colData(obj@colMaps[[which(mask)]]) = DF
              }
          }
        } else if(d == 2){
          if (nrow(DF) == ncol(obj@rowMaps[[n]])) {
              mask = (n == names(obj@rowMaps))
              if (sum(mask) == 1) {
                colData(obj@rowMaps[[which(mask)]]) = DF
              }
          }
        }
      }
  }
  return(obj)
}

.insert_mapping <- function(obj, val, d){
  val = .check_if_mapping_list(val)
  .validate_names(val)

  map_types <- switch(d, rowMapTypes(obj), colMapTypes(obj))

  val = .coerce_input_to_SE(val)

  if (length(val) == 0) {
      val = SimpleList()
  } else{

    val <- sapply(names(val), function(n){
      v = val[[n]]
      if(dim(v)[2] != dim(obj)[d]){
        err = sprintf("ncol(val) must equal %s.\n", dim(obj)[d])
        stop(err)
      }
      colnames(v) <- dimnames(obj)[[d]]
      if(is.null(rownames(v)))
        rownames(v) <- 1:NROW(v)

      if(is.null(S4Vectors::metadata(v)$type))
          v <- .set_map_type(v, map_types[[n]])

      return(v)
    }, simplify = FALSE)
  }

  if(d == 1){
    obj@rowMaps <- as(val, "SimpleList")
  } else if(d == 2){
    obj@colMaps <- as(val, "SimpleList")
  }

  return(obj)
}

.validate_names <- function(val){
  if(any(names(val) == "")){
    par_func = as.character(sys.call(-1)[1])
    err = sprintf("Values passed to '%s' cannot be unnamed.\n", par_func)
    stop(err)
  }

  if(any(duplicated(names(val)))){
    par_func = as.character(sys.call(-1)[1])
    err = sprintf("Values passed to '%s' have duplicate names.\n", par_func)
    stop(err)
  }
  return
}

.coerce_input_to_SE <- function(val){
  val = val[sapply(val, function(v){!is.null(v)})]
  val = lapply(val, function(M) {
      if (any(is(M) == "SummarizedExperiment")) {
          return(M)
      } else if (is.matrix(M) | is.sparseMatrix(M)) {
          M = SummarizedExperiment(assays=list(X=M))
          return(M)
      } else {
        M = as.matrix(val)
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

  return(val)
}

.check_if_mapping_list <- function(val){
  err = sprintf("New mappings must be a named list.\n")
  if( !(class(val) %in% c("list", "SimpleList")) )
    stop(err)
  if(is.null(names(val)))
    stop(err)

  val = as(val, "SimpleList")
  return(val)
}

.set_map_type <- function(val, map_type = NULL, force_embed = FALSE){

  if(is.null(map_type))
    S4Vectors::metadata(val)$type <- ifelse(NROW(val) <= 3, "embedding", "generic")
  else{
    S4Vectors::metadata(val)$type <- map_type
    if(force_embed)
      S4Vectors::metadata(val)$type <- ifelse(NROW(val) <= 3, "embedding", map_type)
  }

  return(val)
}
