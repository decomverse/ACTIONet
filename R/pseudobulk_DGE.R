create_formula <- function(vars){
  fml = Reduce(function(x,y) paste(x,y,sep = " + "), vars)
  fml = paste("~", fml)
  fml = as.formula(fml)
  return(fml)
}

get.pseudobulk.SE <- function(ace, batch_attr, ensemble = FALSE, bins = 20, assay = "counts", colData = NULL, pseudocount = 0, with_S = FALSE, with_E = FALSE, with_V = FALSE, BPPARAM = SerialParam()){

  counts.mat = SummarizedExperiment::assays(ace)[[assay]]
  IDX = .get_attr_or_split_idx(ace, batch_attr)
  sample_names = names(IDX)
  counts.list = lapply(IDX, function(idx) counts.mat[, idx, drop = FALSE])

  # counts.list = lapply(counts.list, as, "dgCMatrix")
  se.assays = list()

  S0 = sapply(counts.list, fastRowSums) + pseudocount
  se.assays$counts = S0

  E0 = sapply(counts.list, fastRowMeans)

  se.assays$mean = E0

  V0 = sapply(counts.list, fastRowVars)
  se.assays$var = V0

  if(ensemble == TRUE){
    if(!any(with_S, with_E, with_V)){
      err = sprintf("No ensemble assays to make.\n")
      stop(err)
    }
    mr.assays = make.ensemble.assays(counts.list = counts.list, bins = bins, pseudocount = pseudocount, with_S = with_S, with_E = with_E, with_V = with_V, BPPARAM = BPPARAM)
    se.assays = c(se.assays, mr.assays$assays)
  }

  n_cells = sapply(counts.list, NCOL)
  nnz_feat_mean = sapply(counts.list, function(X) mean(fastColSums(X > 0)) )
  cd = data.frame(n_cells = n_cells, nnz_feat_mean = nnz_feat_mean, sample = factor(sample_names))

  if(!is.null(colData)){
    md = colData[colData[[batch_attr]] %in% sample_names, , drop = FALSE]
    md = md[match(sample_names, md[[batch_attr]]), ,  drop = FALSE]
    cd = data.frame(md, cd)
  }
  rownames(cd) = sample_names
  cd = droplevels(cd)
  se = SummarizedExperiment(assays = se.assays, colData = cd, rowData = rowData(ace))

  if(ensemble == TRUE){
    metadata(se)$bins = bins
  }
  invisible(gc())
  return(se)
}

make.ensemble.assays <- function(counts.list, bins, pseudocount, with_S = FALSE, with_E = FALSE, with_V = FALSE, BPPARAM = SerialParam()){

  mr.lists = bplapply(names(counts.list), function(n){
    S = Matrix::t(counts.list[[n]])
    bin_IDX = round(seq(1, (1 - bins^-1)*nrow(S), length.out = bins))

      out = list()

    if(with_S == TRUE){
      S.sorted = apply(S, 2, function(s) cumsum(sort(s)))
      cs = S.sorted[nrow(S.sorted), ]
      S.prefix_sum = sapply(bin_IDX, function(idx) cs - S.sorted[idx, ])
      out = list(Sp = S.prefix_sum)
    }

    if(with_E == TRUE){
      E.sorted = apply(S, 2, function(s) rev(cumsum(sort(s, decreasing = T)) / seq.int(length(s)) ) )
      E.prefix_sum = sapply(bin_IDX, function(idx) E.sorted[idx, ])
      out$Ep = E.prefix_sum
    }
    if(with_V == TRUE){
      V.sorted = apply(S, 2, function(s) rev(ACTIONet::roll_var(sort(s, decreasing = T))))
      V.sorted[is.na(V.sorted)] = 0
      V.prefix_sum = sapply(bin_IDX, function(idx) V.sorted[idx, ])
      out$Vp = V.prefix_sum
    }

    return(out)
  }, BPPARAM = BPPARAM)

  mr.assays = list()

  if(with_S == TRUE){
    S.list = lapply(1:bins, function(i){
      sapply(mr.lists, function(L) L$Sp[, i]) + pseudocount
    })
    names(S.list) = paste0("S",1:bins)
    mr.assays = c(mr.assays, S.list)
  }

  if(with_E == TRUE){
    E.list = lapply(1:bins, function(i){
      sapply(mr.lists, function(L) L$Ep[, i])
    })
    names(E.list) = paste0("E",1:bins)
    mr.assays = c(mr.assays, E.list)
  }

  if(with_V == TRUE){
    V.list = lapply(1:bins, function(i){
      sapply(mr.lists, function(L) L$Vp[, i])
    })
    names(V.list) = paste0("V",1:bins)
    mr.assays = c(mr.assays, V.list)
  }

  mr.out = list(assays = mr.assays)

  invisible(gc())
  return(mr.out)
}

run.ensemble.pseudobulk.DESeq <- function(se, design, bins = NULL, bins_use = NULL, slot_prefix = "S", p_adj_method = "fdr", BPPARAM = SerialParam()){
  ACTIONet:::.check_and_load_package("DESeq2")

  if(is.null(bins)){
    bins = metadata(se)$bins
  }

  if(is.null(bins_use))
    bins_use = 1:bins

  dds_out = bplapply(paste0(slot_prefix, 1:bins), function(S){
    cts = assays(se)[[S]]
    invisible({
    dds = DESeqDataSetFromMatrix(cts, design = design, colData = colData(se), rowData = rowData(se))
    dds = DESeq(dds, betaPrior = F)
    })
    return(dds)
  }, BPPARAM = BPPARAM)

  dds_res = lapply(dds_out, function(dds){
    dds_res = results(dds)
    out = data.frame(rowData(dds), dds_res)
    out = out[,c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "dispOutlier")]
    return(out)
  })

  outlier_count = ACTIONet::fast_row_sums(sapply(dds_res, function(dds) dds$dispOutlier))

  mat.bm = sapply(dds_res, function(dds) dds$baseMean)
  bm_mean = rowMeans(mat.bm, na.rm = T)

  mat.lfc = sapply(dds_res, function(dds) dds$log2FoldChange)
  lfc_mean = apply(mat.lfc, 1, function(r) mean(r, na.rm = T, trim = 0.25))

  mat.se = sapply(dds_res, function(dds) dds$lfcSE)
  se_mean = apply(mat.se, 1, function(r) mean(r, na.rm = T, trim = 0.25))

  mat.stat = sapply(dds_res, function(dds) dds$stat)
  stat_mean = rowMeans(mat.stat, na.rm = T)

  mat.pval = sapply(dds_res, function(dds) dds$pvalue)
  mat.pval[is.na(mat.pval)] = 1
  mat.pval = Matrix::t(-log10(mat.pval))
  corr.pvals = ACTIONet::combine.logPvals(mat.pval, base = 10)
  corr.pvals = 10^-corr.pvals

  res = data.frame(rowData(se), baseMean = bm_mean, log2FCMean = lfc_mean, lfcSEMean = se_mean, statMean = stat_mean, pvalue = corr.pvals, padj = p.adjust(corr.pvals, method = p_adj_method), dispOutlierFrac = outlier_count/bins)
  invisible(gc())
  return(res)
}

run.ensemble.pseudobulk.Limma <- function(se, design, bins = NULL, bins_use = NULL, phenotype.name = NULL, min.covered.samples = 2, p_adj_method = "fdr",BPPARAM = SerialParam()){
  ACTIONet:::.check_and_load_package(c("SummarizedExperiment","limma", "edgeR"))

  if(class(se) != "SummarizedExperiment")
    stop("se must be an object of type 'SummarizedExperiment'.")

  if(is.null(bins)){
    bins = metadata(se)$bins
  }

  if(is.null(bins_use))
    bins_use = 1:bins

  if(is(design, "formula")){
    design.mat = model.matrix(design, data = colData(se))
  } else if(is.matrix(design)){
    design.mat = design
  }

  if(is.null(phenotype.name)){
    phenotype.name = colnames(design.mat)[ncol(design.mat)]
  }

  out = bplapply(bins_use, function(i){
    print(i)
    E = assays(se)[[paste0("E",i)]]
    V = assays(se)[[paste0("V",i)]]
    W = 1/V

    slot_E = paste0("E",i)
    slot_V = paste0("V",i)

    tbl = variance.adjusted.DE.Limma(se, slot_E = slot_E, slot_V = slot_V, design = design.mat, phenotype.name = phenotype.name, min.covered.samples = min.covered.samples)

  }, BPPARAM = BPPARAM)

  common = lapply(out, function(tbl) rownames(tbl))
  common = Reduce(intersect, common)
  out = lapply(out, function(tbl) tbl[rownames(tbl) %in% common,])

  mat.bm = sapply(out, function(tbl) tbl$AveExpr)
  bm_mean = rowMeans(mat.bm, na.rm = T)

  mat.lfc = sapply(out, function(tbl) tbl$logFC)
  lfc_mean = apply(mat.lfc, 1, function(r) mean(r, na.rm = T, trim = 0.25))
  # lfc_max = apply(mat.lfc, 1, function(r) r[which.max(abs(r))])

  mat.t = sapply(out, function(tbl) tbl$t)
  t_mean = rowMeans(mat.t, na.rm = T)

  mat.B = sapply(out, function(tbl) tbl$B)
  B_mean = rowMeans(mat.B, na.rm = T)

  mat.pval = sapply(out, function(tbl) tbl$P.Value)
  mat.pval[is.na(mat.pval)] = 1
  mat.pval = Matrix::t(-log10(mat.pval))
  corr.pvals = ACTIONet::combine.logPvals(mat.pval, base = 10)
  corr.pvals = 10^-corr.pvals

  rdat = rowData(se)[rownames(se) %in% common, ]
  res = data.frame(rdat, aveExpr = bm_mean, log2FCMean = lfc_mean, tStatMean = t_mean, BStatMean = B_mean, pvalue = corr.pvals, padj = p.adjust(corr.pvals, method = p_adj_method))
  invisible(gc())
  return(res)
}

variance.adjusted.DE.Limma <- function(se, design, slot_E = "counts", slot_V = NULL, W.mat = NULL, phenotype.name = NULL, min.covered.samples = 3){

  ACTIONet:::.check_and_load_package(c("SummarizedExperiment","limma", "edgeR"))

  if(class(se) != "SummarizedExperiment")
    stop("se must be an object of type 'SummarizedExperiment'.")

  E = SummarizedExperiment::assays(se)[[slot_E]]

  if(!is.null(W.mat)){
    W = W.mat
  } else if(!is.null(slot_V)){
    V = SummarizedExperiment::assays(se)[[slot_V]]
    W = 1/V
  }

  if(is(design, "formula")){
    design.mat = model.matrix(design, data = colData(se))
  } else if(is.matrix(design)){
    design.mat = design
  }

  if(is.null(phenotype.name)){
    phenotype.name = colnames(design.mat)[ncol(design.mat)]
  }

  W.masked = t(apply(W, 1, function(w) {
    mask = (!is.na(w) & is.finite(w) & (w != 0))

    upper = median(w[mask] + 3*mad(w[mask]))
    w[w > upper] = 0

    lower = median(w[mask] - 3*mad(w[mask]))
    w[w < lower] = 0

    w[!mask] = 0
    return(w)
  }))
  W.masked[W.masked == 0] = 1e-16

  selected.vars = ACTIONet::fast_column_sums(design.mat > 0) >= min.covered.samples

  selected.genes = setdiff(1:nrow(W.masked), which(ACTIONet::fast_row_sums(apply(design.mat[, selected.vars], 2, function(x) {
    mm =  (x > 0)
    v = as.numeric( (ACTIONet::fast_row_sums(W.masked[, mm]) == 0) > 0 )
    return(v)
  })) > 0 ))

  suppressWarnings( {fit <- lmFit(E[selected.genes, ], design = design.mat[, selected.vars], weights = W.masked[selected.genes, ])} )

  suppressWarnings({contrast.mat <- makeContrasts(contrasts = phenotype.name, levels = design.mat)})

  cfit <- contrasts.fit(fit, contrast.mat[selected.vars])
  suppressWarnings(efit <- eBayes(cfit, trend = F, robust = T))
  tbl <- topTable(efit, number = Inf, sort.by = "none")
  rownames(tbl) = rownames(se)[selected.genes]
  return(tbl)
}
