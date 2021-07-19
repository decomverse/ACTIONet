.check_and_load_package <- function(pkg_names) {
    for (pk in pkg_names) {
        if (!require(pk, character.only = T)) {
            err = sprintf("Package '%s' is not installed.\n", pk)
            stop(err)
        }
    }
}


    return(labels)
}

.tscalet <- function(A, center = TRUE, scale = TRUE){
  A = Matrix::t(scale(Matrix::t(A), center = center, scale = scale))
  return(A)
}

.preprocess_design_matrix_and_var_names <- function(design_mat, variable_name = NULL){
  if (is.null(variable_name)) {
    variable_name = colnames(design_mat)[ncol(design_mat)]
  }
  vn_idx = which(variable_name == colnames(design_mat))[1]
  colnames(design_mat) = make.names(colnames(design_mat), unique = TRUE, allow_ = FALSE)
  variable_name = colnames(design_mat)[vn_idx]

  out = list(design_mat = design_mat, variable_name = variable_name)
  return(out)
}
