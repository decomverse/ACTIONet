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

rescale.matrix <- function(
  S,
  log_scale = FALSE,
  median_scale = FALSE
) {

    cs = Matrix::colSums(S)
    cs[cs == 0] = 1
    B = Matrix::t(Matrix::t(S) / cs)

    if (median_scale == TRUE){
      B = B * median(cs)
    }

    if (log_scale == TRUE) {
        B = log1p(B)
    }

    return(B)
}
