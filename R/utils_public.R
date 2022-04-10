#' @export
orthoProject <- function(A, S, prenorm = F, postnorm = F) {
    if (prenorm) {
        A <- scale(A)
        S <- scale(S)
    }
    A_r <- A - S %*% MASS::ginv(t(S) %*% S) %*% (t(S) %*% A)
    if (postnorm) {
        A_r <- scale(A_r)
    }
    return(A_r)
}

#' @export
normalize.matrix <- function(
  S,
  log_transform = FALSE,
  scale_param = NULL
) {

  if (!is.matrix(S) && !ACTIONetExperiment:::is.sparseMatrix(S)) {
    err = sprintf("`S` must be `matrix` or `sparseMatrix`.\n")
    stop(err)
  }

  if (is.null(scale_param)) {
    scale_param = 1
  } else if ( !is.function(scale_param) && !is.numeric(scale_param) ) {
    err = sprintf("`scale_param` must be `function` or `numeric`.\n")
    stop(err)
  } else if ( !(length(scale_param) == NCOL(S)) &&  !(length(scale_param) == 1) ) {
    err = sprintf("`scale_param` must be of length 1 or `NCOL(S)`.\n")
  }

  if(ACTIONetExperiment:::is.sparseMatrix(S) && !is(S, "dgCMatrix")) {
    S = as(S, "dgCMatrix")
  }

  cs = Matrix::colSums(S)
  cs[cs == 0] = 1
  B = Matrix::t(Matrix::t(S) / cs)

  if (is.function(scale_param)){
    B = B * scale_param(cs)
  } else {
    if(length(scale_param) > 1){
      B = Matrix::t(Matrix::t(B) * scale_param)
    } else {
      B = B * scale_param
    }
  }

  if (log_transform == TRUE) {
    B = log1p(B)
  }

  if(ACTIONetExperiment:::is.sparseMatrix(S) && !is(S, "dgCMatrix")) {
    S = as(S, "dgCMatrix")
  }

  return(B)
}


.groupedRowSums <- function(S, group_vec) {

  if(ACTIONetExperiment:::is.sparseMatrix(S)) {
    mat = compute_grouped_rowsums(S, sample_assignments = group_vec)
  } else {
    mat = compute_grouped_rowsums_full(S, sample_assignments = group_vec)
  }
  return(mat)
}


.groupedRowMeans <- function(S, group_vec) {

  if(ACTIONetExperiment:::is.sparseMatrix(S)) {
    mat = compute_grouped_rowmeans(S, sample_assignments = group_vec)
  } else {
    mat = compute_grouped_rowmeans_full(S, sample_assignments = group_vec)
  }
  return(mat)
}


#' @export
aggregateMatrix <- function(S, group_vec, method = c("sum", "mean")) {

  method = match.arg(method, several.ok = FALSE)

  if(any(is.na(group_vec))) {
    err = sprintf("'NA' values in 'group_vec'.\n")
    stop(err)
  }

  lf <- factor(.validate_attr(S, attr = group_vec, obj_name = "S", attr_name = "group_vec"))
  labels = as.numeric(lf)
  keys = levels(lf)

  if(ACTIONetExperiment:::is.sparseMatrix(S) && !is(S, "dgCMatrix")) {
    S = as(S, "dgCMatrix")
  } 

  if (method == "sum") {
    mat = .groupedRowSums(S, labels)
  } else if (method == "mean") {
    mat = .groupedRowMeans(S, labels)
  }

  colnames(mat) <- keys
  rownames(mat) = rownames(S)
  return(mat)
}