#' Imputing expression of genes by interpolating over archetype profile
#'
#' @param ace ACTIONet output
#' @param genes List of genes to impute
#'
#' @return A matrix of imputed expression values
#'
#' @examples
#' expression_imputed = impute.genes.using.archetype(ace, genes)
#' @export
impute.genes.using.archetypes <- function(ace, genes) {

    features_use = .preprocess_annotation_features(ace, features_use = features_use)
    genes = intersect(unique(genes), rownames(ace))


    Z = rowMaps(ace)[["archetype_gene_profile"]][genes, ]
    H = Matrix::t(colMaps(ace)[["H_unified"]])

    expression_imputed = t(Z %*% H)
    colnames(expression_imputed) = genes

    return(expression_imputed)
}


#' Imputing expression specificity of genes by interpolating over archetype profile
#'
#' @param ace ACTIONet output
#' @param genes List of genes to impute
#'
#' @return A matrix of imputed expression values
#'
#' @examples
#' expression_imputed = impute.genes.using.archetype(ace, genes)
#' @export
impute.specific.genes.using.archetypes <- function(ace, genes) {

    features_use = .preprocess_annotation_features(ace, features_use = features_use)
    genes = intersect(unique(genes), rownames(ace))


    Z = log1p(rowMaps(ace)[["unified_feature_specificity"]][genes, ])
    H = Matrix::t(colMaps(ace)[["H_unified"]])

    expression_imputed = t(Z %*% H)
    colnames(expression_imputed) = genes

    return(expression_imputed)
}

#' Gene expression imputation using network diffusion.
#'
#' @param ace Input results to be clustered
#' (alternatively it can be the ACTIONet igraph object)
#' @param genes The list of genes to perform imputation for.
#' @param alpha_val Depth of diffusion between (0, 1).
#' The larger it is, the deeper the diffusion, which results in less nonzeros (default = 0.85).
#' @param thread_no Number of parallel threads
#' @param diffusion_iters Number of diffusion iterations (default = 5)
#' @param assay_name Slot in the ace object with normalized counts.
#'
#' @return Imputed gene expression matrix. Column names are set with imputed genes names and rows are cells.
#'
#' @examples
#' imputed.genes = impute.genes.using.ACTIONet(ace, c('CD14', 'CD19', 'CD3G'))
#' plot.ACTIONet.gradient(ace, imputed.genes[, 1])
#' @export
impute.genes.using.ACTIONet <- function(ace, genes, features_use = NULL, alpha_val = 0.85, thread_no = 0,
    diffusion_iters = 5, assay_name = "logcounts") {

    features_use = .preprocess_annotation_features(ace, features_use = features_use)

    genes = unique(genes)

    matched.genes = intersect(genes, features_use)
    matched.idx = match(matched.genes, features_use)

    # Smooth/impute gene expressions
    if (!(assay_name %in% names(SummarizedExperiment::assays(ace)))) {
      err = sprintf("%s is not in assays of ace\n", assay_name)
      stop(err)
    }

    expression_raw = SummarizedExperiment::assays(ace)[[assay_name]][matched.idx, , drop = FALSE]
    if(is.sparseMatrix(expression_raw)){
      expression_raw = U = Matrix::t(as(expression_raw, "dgTMatrix"))
      U[U < 0] = 0
      cs = fastColSums(U)
      U = Matrix::sparseMatrix(i = U@i + 1, j = U@j + 1, x = U@x/cs[U@j + 1], dims = dim(U))
    } else{
      expression_raw = U = Matrix::t(expression_raw)
      U[U < 0] = 0
      cs = fastColSums(U)
      U = expression_raw / fastColSums(expression_raw)
    }
    U = U[, cs > 0]
    gg = matched.genes[cs > 0]

    # if (length(matched.idx) > 1) {
    #     expression_raw = Matrix::t(as(SummarizedExperiment::assays(ace)[[assay_name]][matched.idx,
    #         ], "dgTMatrix"))
    #     U = expression_raw
    #     U[U < 0] = 0
    #     cs = ACTIONet::fastColSums(U)
    #     U = Matrix::sparseMatrix(i = U@i + 1, j = U@j + 1, x = U@x/cs[U@j + 1], dims = dim(U))
    #     U = U[, cs > 0]
    #     gg = matched.genes[cs > 0]
    # } else {
    #     expression_raw = SummarizedExperiment::assays(ace)[[assay_name]][matched.idx,]
    #     U = expression_raw/sum(expression_raw)
    #     gg = matched.genes
    # }

    # Perform network-diffusion
    G = colNets(ace)$ACTIONet
    expression_imputed = compute_network_diffusion(G, as(U, "sparseMatrix"), thread_no = thread_no, alpha = alpha_val, max_it = diffusion_iters)

    expression_imputed[is.na(expression_imputed)] = 0


    # Rescale the baseline expression of each gene
    expression_imputed = sapply(1:dim(expression_imputed)[2], function(col) {
        x = expression_raw[, col]
        y = expression_imputed[, col]

        x.Q = quantile(x, 1)
        y.Q = quantile(y, 1)

        if (y.Q == 0) {
            return(array(0, length(x)))
        }

        y = y * x.Q/y.Q
        y[y > max(x)] = max(x)

        return(y)
    })

    colnames(expression_imputed) = gg

    return(expression_imputed)
}
