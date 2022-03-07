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

  if (is.null(scale_param)) {
    scale_param = 1
  } else if ( !is.function(scale_param) && !is.numeric(scale_param) ) {
    err = sprintf("`scale_param` must be `function` or `numeric`.\n")
    stop(err)
  } else if ( !(length(scale_param) == NCOL(S)) &&  !(length(scale_param) == 1) ) {
    err = sprintf("`scale_param` must be of length 1 or `NCOL(S)`.\n")
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

  return(B)
}
