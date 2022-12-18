#' @export
imputeGenes <- function(
  ace,
  genes,
  algorithm = c("actionet", "action", "pca"),
  features_use = NULL,
  alpha = 0.9,
  diffusion_algorithm = "pagerank_sym",
  thread_no = 0,
  diffusion_it = 5,
  force_reimpute = FALSE,
  assay_name = "logcounts",
  reduction_slot = "ACTION",
  net_slot = "ACTIONet"
) {

    features_use <- .get_feature_vec(ace, features_use = features_use)
    matched_feat <- intersect(unique(genes), features_use)
    idx_feat <- match(matched_feat, features_use)

    algorithm <- tolower(algorithm)
    algorithm <- match.arg(algorithm)

    expr_raw = SummarizedExperiment::assays(ace)[[assay_name]][idx_feat, , drop = FALSE]

    if (algorithm == "pca") {

        V_slot <- sprintf("%s_V", reduction_slot)
        if (!(V_slot %in% names(rowMaps(ace)))) {
            err <- sprintf("`%s` does not exist in rowMaps(ace). Run `reduce.ace()`.", V_slot)
            stop(err)
        }

        smooth_red_name = sprintf("%s_smooth", reduction_slot)
        smooth_U_name = sprintf("%s_U", reduction_slot)
        if ( !(smooth_red_name %in% names(colMaps(ace)) || smooth_U_name %in% names(rowMaps(ace)) || force_reimpute == TRUE) ){
            out <- .smoothPCs(
              ace = ace,
              diffusion_algorithm = diffusion_algorithm,
              alpha = alpha,
              diffusion_it = diffusion_it,
              reduction_slot = reduction_slot,
              net_slot = net_slot,
              thread_no = thread_no,
              return_raw = TRUE
            )

            H <- out$H
            W <- out$SVD.out$u
            W <- W[idx_feat, , drop = FALSE]
        } else {
            W <- rowMaps(ace)[[smooth_U_name]][idx_feat, , drop = FALSE]
            H <- colMaps(ace)[[smooth_red_name]]
        }

        expr_imp <- W %*% Matrix::t(H)
        expr_imp[expr_imp < 0] <- 0

    } else if (algorithm == "action") { # TODO: Fix this!! We need to also impute C. What alpha values?
        if (!("archetype_footprint" %in% names(colMaps(ace))) | (force_reimpute == TRUE)) {

            H <- networkDiffusion(
              obj = ace,
              scores = colMaps(ace)[["H_unified"]],
              algorithm = diffusion_algorithm,
              alpha = alpha,
              thread_no = thread_no,
              max_it = diffusion_it,
              net_slot = net_slot
            )

        } else {
            H <- ace$archetype_footprint
        }
        C <- colMaps(ace)$C_unified
        W <- as.matrix(expr_raw %*% C)
        expr_imp <- W %*% Matrix::t(H)

    } else {
        expr_imp <- networkDiffusion(
          obj = ace,
          scores = Matrix::t(expr_raw),
          algorithm = diffusion_algorithm,
          alpha = alpha,
          thread_no = thread_no,
          max_it = diffusion_it,
          net_slot = net_slot
        )
        expr_imp = Matrix::t(expr_imp)
    }

    # Re-scale expression of genes
    m1 <- apply(expr_raw, 1, max)
    m2 <- apply(expr_imp, 1, max)
    ratio <- m1 / m2
    ratio[m2 == 0] <- 1
    D <- Matrix::Diagonal(nrow(expr_imp), ratio)
    expr_imp <- Matrix::t(as.matrix(D %*% expr_imp))

    colnames(expr_imp) <- matched_feat
    return(expr_imp)
}


#' Imputing expression of genes by interpolating over archetype profile
#'
#' @param ace ACTIONet output
#' @param genes List of genes to impute
#' @param features_use A vector of features of length NROW(ace) or the name of a column of rowData(ace) containing the genes given in 'genes'.
#'
#' @return A matrix of imputed expression values
#'
#' @examples
#' expression_imputed <- impute.genes.using.archetype(ace, genes)
#' @export
impute.genes.using.archetypes <- function(ace, genes, features_use = NULL) {

    features_use <- .get_feature_vec(ace, features_use = features_use)
    matched_feat <- intersect(unique(genes), features_use)
    idx_feat <- match(matched_feat, features_use)

    Z <- rowMaps(ace)[["archetype_gene_profile"]][idx_feat, , drop = FALSE]
    H <- Matrix::t(colMaps(ace)[["H_unified"]])

    expression_imputed <- Matrix::t(Z %*% H)
    colnames(expression_imputed) <- matched_feat

    return(expression_imputed)
}


#' Imputing expression specificity of genes by interpolating over archetype profile
#'
#' @param ace ACTIONet output
#' @param genes List of genes to impute
#' @param features_use A vector of features of length NROW(ace) or the name of a column of rowData(ace) containing the genes given in 'genes'.
#'
#' @return A matrix of imputed expression values
#'
#' @examples
#' expression_imputed <- impute.genes.using.archetype(ace, genes)
#' @export
impute.specific.genes.using.archetypes <- function(ace, genes, features_use = NULL) {

    features_use <- .get_feature_vec(ace, features_use = features_use)
    matched_feat <- intersect(unique(genes), features_use)
    idx_feat <- match(matched_feat, features_use)

    Z <- log1p(rowMaps(ace)[["unified_feature_specificity"]][idx_feat, , drop = FALSE])
    H <- Matrix::t(colMaps(ace)[["H_unified"]])

    expression_imputed <- Matrix::t(Z %*% H)
    colnames(expression_imputed) <- matched_feat

    return(expression_imputed)
}

#' Gene expression imputation using network diffusion.
#'
#' @param ace ACTIONetExperiment object containing output of 'run.ACTIONet()'.
#' @param genes The list of genes to perform imputation on.
#' @param features_use A vector of features of length NROW(ace) or the name of a column of rowData(ace) containing the genes given in 'genes'.
#' @param alpha Depth of diffusion between (0, 1).
#' The larger it is, the deeper the diffusion, which results in less nonzeros (default = 0.85).
#' @param thread_no Number of parallel threads
#' @param diffusion_it Number of diffusion iterations (default = 5)
#' @param assay_name Slot in the ace object with normalized counts.
#'
#' @return Imputed gene expression matrix. Column names are set with imputed genes names and rows are cells.
#'
#' @examples
#' imputed.genes <- impute.genes.using.ACTIONet(ace, c("CD14", "CD19", "CD3G"))
#' plot.ACTIONet.gradient(ace, imputed.genes[, 1])
#' @export
impute.genes.using.ACTIONet <- function(
  ace,
  genes,
  features_use = NULL,
  alpha = 0.85,
  thread_no = 0,
  diffusion_it = 5,
  assay_name = "logcounts"
) {

    features_use <- .get_feature_vec(ace, features_use = features_use)

    matched_feat <- intersect(unique(genes), features_use)
    idx_feat <- match(matched_feat, features_use)

    # Smooth/impute gene expressions
    if (!(assay_name %in% names(SummarizedExperiment::assays(ace)))) {
        err <- sprintf("%s is not in assays of ace\n", assay_name)
        stop(err)
    }

    expression_raw <- SummarizedExperiment::assays(ace)[[assay_name]][idx_feat, ,
        drop = FALSE
    ]
    if (ACTIONetExperiment:::is.sparseMatrix(expression_raw)) {
        expression_raw <- U <- Matrix::t(as(expression_raw, "dgTMatrix"))
        U[U < 0] <- 0
        cs <- Matrix::colSums(U)
        U <- Matrix::sparseMatrix(
            i = U@i + 1,
            j = U@j + 1,
            x = U@x / cs[U@j + 1],
            dims = dim(U)
        )
    } else {
        expression_raw <- U <- Matrix::t(expression_raw)
        U[U < 0] <- 0
        cs <- Matrix::colSums(U)
        U <- expression_raw / Matrix::colSums(expression_raw)
    }
    U <- U[, cs > 0]
    gg <- matched_feat[cs > 0]

    # Perform network-diffusion
    G <- colNets(ace)$ACTIONet

    expression_imputed <- compute_network_diffusion_approx(
        G = G,
        X0 = Matrix::as.matrix(U),
        thread_no = thread_no,
        alpha = alpha,
        max_it = diffusion_it
    )

    expression_imputed[is.na(expression_imputed)] <- 0

    # Rescale the baseline expression of each gene
    expression_imputed <- sapply(1:dim(expression_imputed)[2], function(col) {
        x <- expression_raw[, col]
        y <- expression_imputed[, col]

        x.Q <- quantile(x, 1)
        y.Q <- quantile(y, 1)

        if (y.Q == 0) {
            return(array(0, length(x)))
        }

        y <- y * x.Q / y.Q
        y[y > max(x)] <- max(x)

        return(y)
    })

    colnames(expression_imputed) <- gg

    return(expression_imputed)
}
