variance.adjusted.DE.Limma.test <- function(
  se,
  design,
  slot_E = "counts",
  slot_V = NULL,
  W_mat = NULL,
  variable_name = NULL,
  min_covered_samples = 3
) {

  .check_and_load_package(c("SummarizedExperiment", "limma"))

  if (class(se) != "SummarizedExperiment")
    stop("se must be an object of type 'SummarizedExperiment'.")

  design_mat = make_design_mat(design, data = SummarizedExperiment::colData(se))

  design_list = .preprocess_design_matrix_and_var_names(design_mat, variable_name)
  design_mat = design_list$design_mat
  variable_name = design_list$variable_name
  selected_vars = fastColSums(design_mat > 0) >= min_covered_samples

  E = SummarizedExperiment::assays(se)[[slot_E]]

  if (!is.null(W_mat)) {
    W = W_mat
  } else if (!is.null(slot_V)) {
    V = SummarizedExperiment::assays(se)[[slot_V]]
    W = 1/V
  } else {
    W = NULL
  }

  if(!is.null(W)){
    W_masked = Matrix::t(apply(W, 1, function(w) {
      mask = (!is.na(w) & is.finite(w) & (w != 0))

      upper = stats::median(w[mask] + 3 * stats::mad(w[mask]))
      w[w > upper] = 0

      lower = stats::median(w[mask] - 3 * stats::mad(w[mask]))
      w[w < lower] = 0

      w[!mask] = 0
      return(w)
    }))
    W_masked[W_masked == 0] = 1e-16

    fil_dm = apply(design_mat[, selected_vars, drop = FALSE], 2, function(x) {
      mm = (x > 0)
      v = as.numeric((fastRowSums(W_masked[, mm, drop = FALSE]) == 0) > 0)
      return(v)
    })

    fil_idx = which(fastRowSums(fil_dm) > 0)

    selected_feats = setdiff(1:NROW(W_masked), fil_idx)
    lm_weights = W_masked[selected_feats, ]
  } else {
    selected_feats = 1:NROW(se)
    lm_weights = NULL
  }

  suppressWarnings({
    fit <- limma::lmFit(
      object = E[selected_feats, ],
      design = design_mat[, selected_vars, drop = FALSE],
      weights = lm_weights
    )
  })

  suppressWarnings({
    contrast_mat <- limma::makeContrasts(contrasts = variable_name, levels = design_mat)
  })

  cfit <- limma::contrasts.fit(fit, contrast_mat[selected_vars])
  suppressWarnings(efit <- limma::eBayes(cfit, trend = FALSE, robust = TRUE))
  tbl <- limma::topTable(efit, number = Inf, sort.by = "none")
  rownames(tbl) = rownames(se)[selected_feats]

  return(tbl)
}
