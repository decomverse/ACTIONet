#' @export
imputeGenes <- function(
  ace,
  genes,
  algorithm = c("actionet", "action", "pca"),
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
    matched.genes <- intersect(unique(genes), features_use)
    matched.idx <- match(matched.genes, features_use)
    # genes <- intersect(genes, rownames(ace))

    algorithm <- tolower(algorithm)
    algorithm <- match.arg(algorithm)

    # S <- assays(ace)[[assay_name]]
    # exp_raw <- S[genes, ]

    exp_raw = SummarizedExperiment::assays(ace)[[assay_name]][matched.idx, , drop = FALSE]

    G <- .validate_net(ace, net_slot)

    if (algorithm == "pca") {
        # V_slot <- sprintf("%s_V", reduction_slot)
        # if (!(V_slot %in% names(rowMaps(ace)))) {
        #     err <- sprintf("%s does not exist in rowMaps(ace). Run `reduce.ace()`.", V_slot)
        #     stop(err)
        # } else {
        #     if (!("SVD_V_smooth" %in% names(colMaps(ace))) | (force_reimpute == TRUE)) {
        #         S_r <- colMaps(ace)[[reduction_slot]]
        #         A <- rowMaps(ace)[[sprintf("%s_A", reduction_slot)]]
        #         B <- colMaps(ace)[[sprintf("%s_B", reduction_slot)]]
        #         V <- rowMaps(ace)[[V_slot]]
        #         sigma <- S4Vectors::metadata(ace)[[sprintf("%s_sigma", reduction_slot)]]
        #         U <- as.matrix(S_r %*% Diagonal(length(sigma), 1 / sigma))
        #         SVD.out <- ACTIONet::perturbedSVD(V, sigma, U, -A, B)
        #         V_svd <- SVD.out$v
        #
        #         # This can also be done with network-regularized SVD directly
        #         V.smooth <- networkDiffusion(
        #           G = G,
        #           scores = V_svd,
        #           algorithm = diffusion_algorithm,
        #           alpha = alpha,
        #           thread_no = thread_no,
        #           max_it = diffusion_it,
        #           res_threshold = 1e-8,
        #           net_slot = NULL
        #         )
        #
        #         H <- V.smooth %*% diag(SVD.out$d)
        #         U <- SVD.out$u
        #         rownames(U) <- rownames(ace)
        #         W <- U[matched.idx, , drop = FALSE]
        #     } else {
        #         W <- rowMaps(ace)[["SVD_U"]][matched.idx, , drop = FALSE]
        #         H <- colMaps(ace)[["SVD_V_smooth"]]
        #     }
        # }
        # exp_imp <- W %*% Matrix::t(H)
        # exp_imp[exp_imp < 0] <- 0

        V_slot <- sprintf("%s_V", reduction_slot)
        if (!(V_slot %in% names(rowMaps(ace)))) {
            err <- sprintf("`%s` does not exist in rowMaps(ace). Run `reduce.ace()`.", V_slot)
            stop(err)
        }

        if (!("SVD_V_smooth" %in% names(colMaps(ace))) | (force_reimpute == TRUE)) {
            # S_r <- colMaps(ace)[[reduction_slot]]
            # A <- rowMaps(ace)[[sprintf("%s_A", reduction_slot)]]
            # B <- colMaps(ace)[[sprintf("%s_B", reduction_slot)]]
            # V <- rowMaps(ace)[[V_slot]]
            # sigma <- S4Vectors::metadata(ace)[[sprintf("%s_sigma", reduction_slot)]]
            # U <- as.matrix(S_r %*% Diagonal(length(sigma), 1 / sigma))
            # SVD.out <- ACTIONet::perturbedSVD(V, sigma, U, -A, B)
            # V_svd <- SVD.out$v
            #
            # # This can also be done with network-regularized SVD directly
            # V.smooth <- networkDiffusion(
            #   G = G,
            #   scores = V_svd,
            #   algorithm = diffusion_algorithm,
            #   alpha = alpha,
            #   thread_no = thread_no,
            #   max_it = diffusion_it,
            #   res_threshold = 1e-8,
            #   net_slot = NULL
            # )

            out <- .smoothPCs(
              ace = ace,
              G = G,
              diffusion_algorithm = diffusion_algorithm,
              alpha = alpha,
              diffusion_it = diffusion_it,
              reduction_slot = reduction_slot,
              thread_no = thread_no,
              return_raw = TRUE
            )

            H <- out$H
            W <- out$SVD.out$u
            W <- W[idx_feat, , drop = FALSE]
        } else {
            W <- rowMaps(ace)[["SVD_U"]][idx_feat, , drop = FALSE]
            H <- colMaps(ace)[["SVD_V_smooth"]]
        }

        exp_imp <- W %*% Matrix::t(H)
        exp_imp[exp_imp < 0] <- 0


    } else if (algorithm == "action") { # TODO: Fix this!! We need to also impute C. What alpha values?
        if (!("archetype_footprint" %in% names(colMaps(ace))) | (force_reimpute == TRUE)) {

            H <- networkDiffusion(
              G = G,
              scores = colMaps(ace)[["H_unified"]],
              algorithm = diffusion_algorithm,
              alpha = alpha,
              thread_no = thread_no,
              max_it = diffusion_it
            )

        } else {
            H <- ace$archetype_footprint
        }
        C <- colMaps(ace)$C_unified
        W <- as.matrix(exp_raw %*% C)
        exp_imp <- W %*% Matrix::t(H)
    } else if (algorithm == "actionet") {
        exp_imp <- networkDiffusion(
          G = G,
          scores = Matrix::t(exp_raw),
          algorithm = diffusion_algorithm,
          alpha = alpha,
          thread_no = thread_no,
          max_it = diffusion_it,
          res_threshold = 1e-8,
          net_slot = NULL
        )
        exp_imp = Matrix::t(exp_imp)
    }

    rownames(exp_imp) <- matched.genes

    # Re-scaling expresion of genes
    m1 <- apply(exp_raw, 1, max)
    m2 <- apply(exp_imp, 1, max)
    ratio <- m1 / m2
    ratio[m2 == 0] <- 1
    D <- Diagonal(nrow(exp_imp), ratio)
    exp_imp <- as.matrix(D %*% exp_imp)

    return(exp_imp)
}


# .impute.features.ACTIONet <- function(){
#
# }
#
# .impute.features.ACTION <- function(){
#
# }
#
# .impute.features.PCA <- function(
#   ace,
#   G,
#   idx_feat,
#   reduction_slot = "ACTION",
#   force_reimpute = FALSE,
#   diffusion_algorithm = "pagerank_sym",
#   alpha = 0.9,
#   diffusion_it = 5,
#   thread_no = 0
# ){
#
#   V_slot <- sprintf("%s_V", reduction_slot)
#   if (!(V_slot %in% names(rowMaps(ace)))) {
#       err <- sprintf("`%s` does not exist in rowMaps(ace). Run `reduce.ace()`.", V_slot)
#       stop(err)
#   }
#
#   if (!("SVD_V_smooth" %in% names(colMaps(ace))) | (force_reimpute == TRUE)) {
#       # S_r <- colMaps(ace)[[reduction_slot]]
#       # A <- rowMaps(ace)[[sprintf("%s_A", reduction_slot)]]
#       # B <- colMaps(ace)[[sprintf("%s_B", reduction_slot)]]
#       # V <- rowMaps(ace)[[V_slot]]
#       # sigma <- S4Vectors::metadata(ace)[[sprintf("%s_sigma", reduction_slot)]]
#       # U <- as.matrix(S_r %*% Diagonal(length(sigma), 1 / sigma))
#       # SVD.out <- ACTIONet::perturbedSVD(V, sigma, U, -A, B)
#       # V_svd <- SVD.out$v
#       #
#       # # This can also be done with network-regularized SVD directly
#       # V.smooth <- networkDiffusion(
#       #   G = G,
#       #   scores = V_svd,
#       #   algorithm = diffusion_algorithm,
#       #   alpha = alpha,
#       #   thread_no = thread_no,
#       #   max_it = diffusion_it,
#       #   res_threshold = 1e-8,
#       #   net_slot = NULL
#       # )
#
#       out <- .smoothPCs(
#         ace = ace,
#         G = G,
#         diffusion_algorithm = diffusion_algorithm,
#         alpha = alpha,
#         diffusion_it = diffusion_it,
#         reduction_slot = reduction_slot,
#         thread_no = thread_no,
#         return_raw = TRUE
#       )
#
#       # out <- list(U = U, SVD.out = SVD.out, V.smooth = V.smooth, H = H)
#
#       H <- out$H
#       W <- out$SVD.out$u
#       rownames(W) <- rownames(ace)
#       W <- W[idx_feat, , drop = FALSE]
#   } else {
#       W <- rowMaps(ace)[["SVD_U"]][idx_feat, , drop = FALSE]
#       H <- colMaps(ace)[["SVD_V_smooth"]]
#   }
#
#   exp_imp <- W %*% Matrix::t(H)
#   exp_imp[exp_imp < 0] <- 0
#
#   return(exp_imp)
# }

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
    # features_use <- .get_feature_vec(ace, features_use = features_use)
    # genes <- intersect(unique(genes), features_use)

    features_use <- .get_feature_vec(ace, features_use = features_use)
    matched.genes <- intersect(unique(genes), features_use)
    matched.idx <- match(matched.genes, features_use)

    Z <- rowMaps(ace)[["archetype_gene_profile"]][matched.idx, , drop = FALSE]
    H <- Matrix::t(colMaps(ace)[["H_unified"]])

    expression_imputed <- Matrix::t(Z %*% H)
    colnames(expression_imputed) <- matched.genes

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
    # features_use <- .get_feature_vec(ace, features_use = features_use)
    # genes <- intersect(unique(genes), features_use)

    features_use <- .get_feature_vec(ace, features_use = features_use)
    matched.genes <- intersect(unique(genes), features_use)
    matched.idx <- match(matched.genes, features_use)

    Z <- log1p(rowMaps(ace)[["unified_feature_specificity"]][matched.idx, , drop = FALSE])
    H <- Matrix::t(colMaps(ace)[["H_unified"]])

    expression_imputed <- Matrix::t(Z %*% H)
    colnames(expression_imputed) <- matched.genes

    return(expression_imputed)
}

#' Gene expression imputation using network diffusion.
#'
#' @param ace Input results to be clustered
#' (alternatively it can be the ACTIONet igraph object)
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
    # genes <- unique(genes)

    matched.genes <- intersect(unique(genes), features_use)
    matched.idx <- match(matched.genes, features_use)

    # Smooth/impute gene expressions
    if (!(assay_name %in% names(SummarizedExperiment::assays(ace)))) {
        err <- sprintf("%s is not in assays of ace\n", assay_name)
        stop(err)
    }

    expression_raw <- SummarizedExperiment::assays(ace)[[assay_name]][matched.idx, ,
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
    gg <- matched.genes[cs > 0]

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
