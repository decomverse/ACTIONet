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

normalize.matrix <- function(
  S,
  log_transform = FALSE,
  scale_factor = NULL
) {

    if (is.null(scale_factor)) {
      scale_factor = 1
    } else if ( (scale_factor != "median") && !is.numeric(scale_factor) ) {
      err = sprintf("'scale_factor' must be 'median' or numeric.\n")
      stop(err)
    }

    cs = Matrix::colSums(S)
    cs[cs == 0] = 1
    B = Matrix::t(Matrix::t(S) / cs)

    if (scale_factor == "median"){
      B = B * median(cs)
    } else {
      B = B * scale_factor
    }
    
    if (log_transform == TRUE) {
        B = log1p(B)
    }

    return(B)
}

# #' @export
# fastRowVars <- function(mat) {
#     mat <- as(mat, "dgTMatrix")
#     E <- fastRowMeans(mat)
#     V <- computeSparseRowVariances(mat@i + 1, mat@x, E, ncol(mat))
#     return(V)
# }
