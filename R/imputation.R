#' Imputing expression of genes by interpolating over archetype profile
#'
#' @param ace ACTIONet output
#' @param genes List of genes to impute
#'
#' @return A matrix of imputed expression values
#'
#' @examples
#' imputed.gene.expression = impute.genes.using.archetype(ace, genes)
#' @export
impute.genes.using.archetypes <- function(ace, genes) {
    require(igraph)

    genes = intersect(unique(genes), rownames(ace))


    Z = rowMaps(ace)[["archetype_gene_profile"]][genes, ]
    H = Matrix::t(colMaps(ace)[["H_unified"]])

    imputed.gene.expression = t(Z %*% H)
    colnames(imputed.gene.expression) = genes

    return(imputed.gene.expression)
}


#' Imputing expression specificity of genes by interpolating over archetype profile
#'
#' @param ace ACTIONet output
#' @param genes List of genes to impute
#'
#' @return A matrix of imputed expression values
#'
#' @examples
#' imputed.gene.expression = impute.genes.using.archetype(ace, genes)
#' @export
impute.specific.genes.using.archetypes <- function(ace, genes) {
    require(igraph)

    genes = intersect(unique(genes), rownames(ace))


    Z = log1p(rowMaps(ace)[["unified_feature_specificity"]][genes, ])
    H = Matrix::t(colMaps(ace)[["H_unified"]])

    imputed.gene.expression = t(Z %*% H)
    colnames(imputed.gene.expression) = genes

    return(imputed.gene.expression)
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
#' @param data_slot Slot in the ace object with normalized counts.
#'
#' @return Imputed gene expression matrix. Column names are set with imputed genes names and rows are cells.
#'
#' @examples
#' imputed.genes = impute.genes.using.ACTIONet(ace, c('CD14', 'CD19', 'CD3G'))
#' plot.ACTIONet.gradient(ace, imputed.genes[, 1])
#' @export
impute.genes.using.ACTIONet <- function(ace, genes, alpha_val = 0.85, thread_no = 8,
    diffusion_iters = 5, data_slot = "logcounts") {
    genes = unique(genes)


    matched.genes = intersect(genes, rownames(ace))
    matched.idx = match(matched.genes, rownames(ace))

    # Smooth/impute gene expressions
    if (!(data_slot %in% names(SummarizedExperiment::assays(ace)))) {
        R.utils::printf("%s is not in assays of ace\n", data_slot)
    }

    if (length(matched.idx) > 1) {
        raw.gene.expression = Matrix::t(as(SummarizedExperiment::assays(ace)[[data_slot]][matched.idx,
            ], "dgTMatrix"))
        U = raw.gene.expression
        U[U < 0] = 0
        cs = fast_column_sums(U)
        U = Matrix::sparseMatrix(i = U@i + 1, j = U@j + 1, x = U@x/cs[U@j + 1], dims = dim(U))
        U = U[, cs > 0]
        gg = matched.genes[cs > 0]
    } else {
        raw.gene.expression = SummarizedExperiment::assays(ace)[[data_slot]][matched.idx,
            ]
        U = raw.gene.expression/sum(raw.gene.expression)
        gg = matched.genes
    }

    # Perform network-diffusion
    G = colNets(ace)$ACTIONet
    imputed.gene.expression = compute_network_diffusion(G, as(U, "sparseMatrix"),
        alpha = alpha_val, max_it = diffusion_iters)

    imputed.gene.expression[is.na(imputed.gene.expression)] = 0


    # Rescale the baseline expression of each gene
    imputed.gene.expression = sapply(1:dim(imputed.gene.expression)[2], function(col) {
        x = raw.gene.expression[, col]
        y = imputed.gene.expression[, col]

        x.Q = quantile(x, 1)
        y.Q = quantile(y, 1)

        if (y.Q == 0) {
            return(array(0, length(x)))
        }

        y = y * x.Q/y.Q

        y[y > max(x)] = max(x)

        return(y)
    })

    colnames(imputed.gene.expression) = gg

    return(imputed.gene.expression)
}
