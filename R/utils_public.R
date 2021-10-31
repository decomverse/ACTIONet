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

    if (is.matrix(S)) {
        cs = Matrix::colSums(S)
        cs[cs == 0] = 1
        B = Matrix::t(Matrix::t(S) / cs)

        if (median_scale == TRUE){
          B = B * median(cs)
        }

        if (log_scale == TRUE) {
            B = log1p(B)
        }

    } else {
        A = as(S, "dgTMatrix")
        cs = Matrix::colSums(S)
        cs[cs == 0] = 1
        x = A@x/cs[A@j + 1]

        if (median_scale == TRUE){
          x = x * median(cs)
        }

        if (log_scale == TRUE) {
            x = log1p(x)
        }
        B = Matrix::sparseMatrix(
          i = A@i + 1,
          j = A@j + 1,
          x = x,
          dims = dim(A)
        )
    }

    return(B)
}
