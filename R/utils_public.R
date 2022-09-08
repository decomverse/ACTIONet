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
normalize.matrix <- function(S,
                             log_transform = FALSE,
                             scale_param = NULL) {
  if (!is.matrix(S) && !ACTIONetExperiment:::is.sparseMatrix(S)) {
    err <- sprintf("`S` must be `matrix` or `sparseMatrix`.\n")
    stop(err)
  }

  if (is.null(scale_param)) {
    scale_param <- 1
  } else if (!is.function(scale_param) && !is.numeric(scale_param)) {
    err <- sprintf("`scale_param` must be `function` or `numeric`.\n")
    stop(err)
  } else if (!(length(scale_param) == NCOL(S)) && !(length(scale_param) == 1)) {
    err <- sprintf("`scale_param` must be of length 1 or `NCOL(S)`.\n")
  }

  if (ACTIONetExperiment:::is.sparseMatrix(S) && !is(S, "dgCMatrix")) {
    S <- as(S, "dgCMatrix")
  }

  # cs = Matrix::colSums(S)
  # cs[cs == 0] = 1
  # B = Matrix::t(Matrix::t(S) / cs)
  if (is.matrix(S)) {
    B <- normalize_mat(S, 1)
  } else {
    B <- normalize_spmat(S, 1)
  }

  if (is.function(scale_param)) {
    lib_sizes <- Matrix::colSums(S)
    B <- B * scale_param(lib_sizes)
  } else {
    if (length(scale_param) > 1) {
      B <- Matrix::t(Matrix::t(B) * scale_param)
    } else {
      B <- B * scale_param
    }
  }

  if (log_transform == TRUE) {
    B <- log1p(B)
  }

  if (ACTIONetExperiment:::is.sparseMatrix(S) && !is(S, "dgCMatrix")) {
    S <- as(S, "dgCMatrix")
  }

  return(B)
}


.groupedRowSums <- function(S, group_vec) {
  if (ACTIONetExperiment:::is.sparseMatrix(S)) {
    mat <- compute_grouped_rowsums(S, sample_assignments = group_vec)
  } else {
    mat <- compute_grouped_rowsums_full(S, sample_assignments = group_vec)
  }
  return(mat)
}


.groupedRowMeans <- function(S, group_vec) {
  if (ACTIONetExperiment:::is.sparseMatrix(S)) {
    mat <- compute_grouped_rowmeans(S, sample_assignments = group_vec)
  } else {
    mat <- compute_grouped_rowmeans_full(S, sample_assignments = group_vec)
  }
  return(mat)
}


.groupedRowVars <- function(S, group_vec) {
  if (ACTIONetExperiment:::is.sparseMatrix(S)) {
    mat <- compute_grouped_rowvars(S, sample_assignments = group_vec)
  } else {
    mat <- compute_grouped_rowvars_full(S, sample_assignments = group_vec)
  }
  return(mat)
}


#' @export
aggregateMatrix <- function(S, group_vec, method = c("sum", "mean", "var")) {
  method <- match.arg(method, several.ok = FALSE)

  if (any(is.na(group_vec))) {
    err <- sprintf("'NA' values in 'group_vec'.\n")
    stop(err)
  }

  lf <- factor(.validate_attr(S, attr = group_vec, obj_name = "S", attr_name = "group_vec"))
  labels <- as.numeric(lf)
  keys <- levels(lf)

  if (ACTIONetExperiment:::is.sparseMatrix(S) && !is(S, "dgCMatrix")) {
    S <- as(S, "dgCMatrix")
  }

  if (method == "sum") {
    mat <- .groupedRowSums(S, labels)
  } else if (method == "mean") {
    mat <- .groupedRowMeans(S, labels)
  } else if (method == "var") {
    mat <- .groupedRowVars(S, labels)
  }

  colnames(mat) <- keys
  rownames(mat) <- rownames(S)
  return(mat)
}

verify_aces <- function(ace1, ace2) {
  ###############################################################
  ###############################################################
  ############## Check ACTION decomposition #####################
  ###############################################################
  ###############################################################

  # check multi-level decomposition
  C1 <- colMaps(ace1)$C_stacked
  H1 <- colMaps(ace1)$H_stacked
  C2 <- colMaps(ace2)$C_stacked
  H2 <- colMaps(ace2)$H_stacked

  n1 <- ncol(H1)
  n2 <- ncol(H2)
  stopifnot(n1 == n2)
  print(sprintf("Same number of multi-level archetypes [%d] (Passed)", n1))

  deltaH <- sum(abs(H1 - H2)) / length(H1)
  stopifnot(deltaH < 1e-5)
  print(sprintf("Delta H (multi-level) = %.2e (Passed)", deltaH))

  deltaC <- sum(abs(C1 - C2)) / length(H1)
  stopifnot(deltaC < 1e-5)
  print(sprintf("Delta C (multi-level) = %.2e (Passed)", deltaC))

  # check multi-level decomposition
  C1 <- colMaps(ace1)$C_unified
  H1 <- colMaps(ace1)$H_unified
  C2 <- colMaps(ace2)$C_unified
  H2 <- colMaps(ace2)$H_unified

  n1 <- ncol(H1)
  n2 <- ncol(H2)
  stopifnot(n1 == n2)
  print(sprintf("Same number of multi-resolution archetypes [%d] (Passed)", n1))

  deltaH <- sum(abs(H1 - H2)) / length(H1)
  stopifnot(deltaH < 1e-5)
  print(sprintf("Delta H (multi-resolution) = %.2e (Passed)", deltaH))

  deltaC <- sum(abs(C1 - C2)) / length(H2)
  stopifnot(deltaC < 1e-5)
  print(sprintf("Delta C (multi-resolution) = %.2e (Passed)", deltaC))

  archs1 <- ace1$assigned_archetype
  archs2 <- ace2$assigned_archetype
  mismatch_perc <- 100 * sum(archs1 != archs2) / ncol(ace)
  stopifnot(mismatch_perc < 0.5)
  sprintf("%.02d %% archetype assignment mismatch (Passed)", mismatch_perc)


  ###############################################################
  ###############################################################
  ########## Check archetype feature specificity  ###############
  ###############################################################
  ###############################################################
  spec1 <- round(ace1$unified_feature_specificity, 3)
  spec2 <- round(ace2$unified_feature_specificity, 3)

  deltaSpec <- sum(abs(spec1 - spec2)) / length(spec1)
  stopifnot(deltaSpec < 1e-5)
  print(sprintf("Delta archetype feature specificity = %.2e (Passed)", deltaSpec))


  ###############################################################
  ###############################################################
  ############## Check network construction #####################
  ###############################################################
  ###############################################################
  net1 <- round(ace1$ACTIONet, 3)
  net2 <- round(ace2$ACTIONet, 3)

  mismatch.edges <- 100 * sum(net1 != net2) / length((net1@i))
  stopifnot(mismatch.edges < 0.5)
  print(sprintf("%.02f %% ACTIONet edges mismatch (Passed)", mismatch.edges))


  ###############################################################
  ###############################################################
  ############# Check network visualization  ####################
  ###############################################################
  ###############################################################
  ## 2D
  coor2D1 <- round(ace1$ACTIONet2D, 3)
  coor2D2 <- round(ace2$ACTIONet2D, 3)

  mismatch.2D <- 100 * sum(coor2D1 != coor2D2) / length(coor2D1)
  stopifnot(mismatch.2D < 0.5)
  print(sprintf("%.02f %% 2D mismatch (Passed)", mismatch.2D))


  ## 3D
  coor3D1 <- round(ace1$ACTIONet3D, 3)
  coor3D2 <- round(ace2$ACTIONet3D, 3)

  mismatch.3D <- 100 * sum(coor3D1 != coor3D2) / length(coor3D1)
  stopifnot(mismatch.3D < 0.5)
  print(sprintf("%.02d %% 3D mismatch (Passed)", mismatch.3D))

  ## Colors
  colors1 <- round(ace1$denovo_color, 3)
  colors2 <- round(ace2$denovo_color, 3)

  mismatch.colors <- 100 * sum(colors1 != colors2) / length(colors1)
  stopifnot(mismatch.colors < 0.5)
  print(sprintf("%.02f %% colors mismatch (Passed)", mismatch.colors))

  print("Congratulations! Your ace objects match perfectly!!")
}

export_minimal_sce <- function(ace, export_logcounts = FALSE) {
  if (export_logcounts) {
    sce <- SingleCellExperiment(assays = list(counts = counts(ace), logcounts = logcounts(ace)))
  } else {
    sce <- SingleCellExperiment(assays = list(counts = counts(ace)))
  }
  colData(sce) <- colData(ace)
  rowData(sce) <- rowData(ace)

  return(sce)
}