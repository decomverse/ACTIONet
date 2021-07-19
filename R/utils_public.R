#' @export
orthoProject <- function(A, S) {
    A = scale(A)
    S = scale(S)
    A_r = A - S %*% MASS::ginv(t(S) %*% S) %*% (t(S) %*% A)
    A_r = scale(A_r)
    return(A_r)
}
