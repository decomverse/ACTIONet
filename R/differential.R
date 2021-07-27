process.var.of.interest <- function(ace, var_of_interest, max.class = 100) {
  if(length(var_of_interest) == 1) {
    if(is.character(var_of_interest)) {
      if(var_of_interest %in% colnames(colData(ace))) {
        v = colData(ace)[[var_of_interest]]
      } else {
        warning(sprintf("Variable %s not found", var_of_interest))
        return(NULL)
      }
    } 
  } else if(length(var_of_interest) == ncol(ace)) {
    v = var_of_interest
  } else {
    warning(sprintf("Unknown var_of_interest"))
    return(NULL)
  }
  
  if( (class(v) == "numeric") | (class(v) == "character"))  {
    uv = sort(unique(v))
    if(max.class < length(uv)) {
      warning(sprintf("Variable %s has %d unique values (max.class = %d)", var_of_interest, length(uv), max.class))
      return(NULL)
    }
    f = factor(v, uv)
  } else if(class(v) == "factor") {
    f = v
  }
  
  f = droplevels(f)
  return(f)
}

findMarkers.ACTIONet <- function(ace, f, out.name = "cond", pos.only = T, blacklist.pattern = "^MT-|^MT[:.:]|^RPS|^RPL|^MALAT") {
  require(ACTIONet)
  
  if(class(f) != "factor") {
    warning("f must be a factor")
    return(ace)
  }
  f = droplevels(f)

  out = compute_cluster_feature_specificity(logcounts(ace), as.numeric(f))
  metadata(ace)[[sprintf("%s_markers_ACTIONet", out.name)]] = out

  scores = as.matrix(out$upper_significance - out$lower_significance)
  colnames(scores) = levels(f)
  rownames(scores) = rownames(ace)
  
  if(pos.only == T)
    scores[scores < 0] = 0 # Only "positive markers" [negative would be markers in other levels]
  
  blacklisted.rows = grep(blacklist.pattern, rownames(ace), ignore.case = T)
  scores[blacklisted.rows, ] = 0
  
  rowMaps(ace)[[sprintf("%s_markers_ACTIONet", out.name)]] = scores
  rowMapTypes(ace)[[sprintf("%s_markers_ACTIONet", out.name)]] = "reduction"
  
  return(ace)
}

findMarkers.wilcox <- function(ace, f, out.name = "cond", pos.only = T, blacklist.pattern = "^MT-|^MT[:.:]|^RPS|^RPL|^MALAT") {
  require(presto)
  
  if(class(f) != "factor") {
    warning("f must be a factor")
    return(ace)
  }
  f = droplevels(f)
  DF = presto::wilcoxauc(logcounts(ace), f)
  metadata(ace)[[sprintf("%s_markers_wilcox", out.name)]] = DF
  DFs = split(DF, DF$group)
  DFs = DFs[levels(f)]  
  scores = do.call(cbind, lapply(DFs, function(subDF) {
    x = subDF$auc[match(rownames(ace), subDF$feature)] - 0.5
    return(x)
  }))
  colnames(scores) = levels(f)
  rownames(scores) = rownames(ace)
  
  if(pos.only == T)
    scores[scores < 0] = 0 # Only "positive markers" [negative would be markers in other levels]
  
  blacklisted.rows = grep(blacklist.pattern, rownames(ace), ignore.case = T)
  scores[blacklisted.rows, ] = 0
  
  rowMaps(ace)[[sprintf("%s_markers_wilcox", out.name)]] = scores
  rowMapTypes(ace)[[sprintf("%s_markers_wilcox", out.name)]] = "reduction"
  
  return(ace)
}

findMarkers.scran <- function(ace, f, out.name = "cond", pos.only = T, blacklist.pattern = "^MT-|^MT[:.:]|^RPS|^RPL|^MALAT") {
  require(scran)
  
  if(class(f) != "factor") {
    warning("f must be a factor")
    return(ace)
  }
  f = droplevels(f)
  
  out = scran::findMarkers(logcounts(ace), f)
  metadata(ace)[[sprintf("%s_markers_scran", out.name)]] = out

  scores = do.call(cbind, lapply(out, function(DF) {
    idx = match(rownames(ace), rownames(DF))
    x = -log10(DF$p.value[idx]) * sign(DF$summary.logFC[idx])
  }))
  colnames(scores) = levels(f)
  rownames(scores) = rownames(ace)
  scores[scores == Inf] = max(scores[!is.infinite(scores)])
  scores[scores == -Inf] = min(scores[!is.infinite(scores)])
  
  if(pos.only == T)
    scores[scores < 0] = 0 # Only "positive markers" [negative would be markers in other levels]

  blacklisted.rows = grep(blacklist.pattern, rownames(ace), ignore.case = T)
  scores[blacklisted.rows, ] = 0
  
  rowMaps(ace)[[sprintf("%s_markers_scran", out.name)]] = scores
  rowMapTypes(ace)[[sprintf("%s_markers_scran", out.name)]] = "reduction"
  
  return(ace)
}

findMarkers.limma <- function(ace, f, out.name = "cond", pos.only = T, blacklist.pattern = "^MT-|^MT[:.:]|^RPS|^RPL|^MALAT") {
  require(limma)
  
  if(class(f) != "factor") {
    warning("f must be a factor")
    return(ace)
  }
  f = droplevels(f)

  design = model.matrix(~0. + f)
  fit0 = limma::lmFit(logcounts(ace), design)
  fit = eBayes(fit0, robust=TRUE, trend=TRUE)
  
  scores = matrix(0, nrow = nrow(ace), ncol = length(levels(f)))
  colnames(scores) = levels(f)
  rownames(scores) = rownames(ace)
  DFs = vector("list", length(levels(f)))
  names(DFs) = levels(f)
  for(i in 1:length(levels(f))) {
    contrast.mat = matrix(-1/(ncol(design)-1), nrow = ncol(design))
    contrast.mat[[i]] = 1
    rownames(contrast.mat) = colnames(design)
    
    fit <- limma::contrasts.fit(fit0, contrasts = contrast.mat)
    fit <- limma::eBayes(fit, robust=TRUE, trend=TRUE)#, trend = TRUE, proportion = 0.05)
    DF = limma::topTable(
      fit = fit,
      number = Inf,
      adjust.method = "BH",
      sort.by = "none"
    )
  
    DFs[[i]] = DF
    scores[rownames(DF), i] = DF$t
  }

  if(pos.only == T)
    scores[scores < 0] = 0 # Only "positive markers" [negative would be markers in other levels]
  
  blacklisted.rows = grep(blacklist.pattern, rownames(ace), ignore.case = T)
  scores[blacklisted.rows, ] = 0
  
  metadata(ace)[[sprintf("%s_markers_limma", out.name)]] = DFs
  rowMaps(ace)[[sprintf("%s_markers_limma", out.name)]] = scores
  rowMapTypes(ace)[[sprintf("%s_markers_limma", out.name)]] = "reduction"
  
  return(ace)
}



