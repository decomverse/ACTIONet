#' @export
findMarkers.ACTIONet <- function(ace,
                                 cluster_attr,
                                 top_genes = 10,
                                 most_specific = FALSE,
                                 features_use = NULL,
                                 feat_subset = NULL,
                                 assay_name = "logcounts",
                                 thread_no = 0,
                                 to_return = c("data.frame", "df", "list")) {
  to_return <- match.arg(to_return)

  sa <- ACTIONetExperiment::get.data.or.split(ace, attr = cluster_attr, to_return = "levels")
  features_use <- .get_feature_vec(ace, features_use = features_use)

  specificity <- clusterFeatureSpecificity(
    obj = ace,
    cluster_attr = sa$index,
    assay_name = assay_name,
    thread_no = thread_no,
    return_raw = TRUE
  )

  feat_spec <- specificity[["upper_significance"]] - specificity[["lower_significance"]]
  feat_spec[feat_spec < 0] <- 0
  rownames(feat_spec) <- features_use
  colnames(feat_spec) <- sa$keys

  if (!is.null(feat_subset)) {
    feat_spec <- feat_spec[rownames(feat_spec) %in% feat_subset, , drop = FALSE]
  }

  if (most_specific == TRUE) {
    W <- select.top.k.features(
      feat_spec,
      top_features = top_genes,
      normalize = FALSE,
      reorder_columns = FALSE
    )
    feat_spec_top <- apply(W, 2, function(v) rownames(W)[order(v, decreasing = TRUE)][1:top_genes])
  } else {
    feat_spec_top <- sapply(colnames(feat_spec), function(type) {
      c <- feat_spec[, type]
      names(head(sort(c, decreasing = TRUE), top_genes))
    })
  }

  df <- data.frame(feat_spec_top)

  if (to_return == "list") {
    return(as.list(df))
  } else {
    return(df)
  }
}

process.var.of.interest <- function(ace, var_of_interest, max_class = 100) {
  if (length(var_of_interest) == 1) {
    if (is.character(var_of_interest)) {
      if (var_of_interest %in% colnames(colData(ace))) {
        v <- colData(ace)[[var_of_interest]]
      } else {
        warning(sprintf("Variable %s not found", var_of_interest))
        return(NULL)
      }
    }
  } else if (length(var_of_interest) == ncol(ace)) {
    v <- var_of_interest
  } else {
    warning(sprintf("Unknown var_of_interest"))
    return(NULL)
  }

  if ((class(v) == "numeric") | (class(v) == "character")) {
    uv <- sort(unique(v))
    if (max_class < length(uv)) {
      warning(sprintf("Variable %s has %d unique values (max_class = %d)", var_of_interest, length(uv), max_class))
      return(NULL)
    }
    f <- factor(v, uv)
  } else if (class(v) == "factor") {
    f <- v
  }

  f <- droplevels(f)
  return(f)
}

computeGeneSpecifity.ACTIONet <- function(ace, f, out_name = "cond", pos_only = T, blacklist_pattern = "^MT-|^MT[:.:]|^RPS|^RPL|^MALAT") {
  require(ACTIONet)
  print("Running computeGeneSpecifity.ACTIONet()")

  if (class(f) != "factor") {
    warning("f must be a factor")
    return(ace)
  }
  f <- droplevels(f)

  S <- logcounts(ace)

  if (is.matrix(S)) {
    out <- compute_cluster_feature_specificity_full(S, as.numeric(f))
  } else {
    out <- compute_cluster_feature_specificity(S, as.numeric(f))
  }
  metadata(ace)[[sprintf("%s_feature_specificity_ACTIONet", out_name)]] <- out

  scores <- as.matrix(out$upper_significance - out$lower_significance)
  colnames(scores) <- levels(f)
  rownames(scores) <- rownames(ace)

  if (pos_only == T) {
    scores[scores < 0] <- 0
  } # Only "positive markers" [negative would be markers in other levels]

  blacklisted.rows <- grep(blacklist_pattern, rownames(ace), ignore.case = T)
  scores[blacklisted.rows, ] <- 0

  rowMaps(ace)[[sprintf("%s_feature_specificity", out_name)]] <- scores
  rowMapTypes(ace)[[sprintf("%s_feature_specificity", out_name)]] <- "reduction"

  return(ace)
}

computeGeneSpecifity.wilcoxon <- function(ace, f, out_name = "cond", pos_only = T, blacklist_pattern = "^MT-|^MT[:.:]|^RPS|^RPL|^MALAT") {
  require(presto)
  print("Running computeGeneSpecifity.wilcox()")

  if (class(f) != "factor") {
    warning("f must be a factor")
    return(ace)
  }
  f <- droplevels(f)
  DF <- presto::wilcoxauc(logcounts(ace), f)
  metadata(ace)[[sprintf("%s_feature_specificity_wilcox", out_name)]] <- DF
  DFs <- split(DF, DF$group)
  DFs <- DFs[levels(f)]
  scores <- do.call(cbind, lapply(DFs, function(subDF) {
    x <- subDF$auc[match(rownames(ace), subDF$feature)] - 0.5
    return(x)
  }))
  colnames(scores) <- levels(f)
  rownames(scores) <- rownames(ace)

  if (pos_only == T) {
    scores[scores < 0] <- 0
  } # Only "positive markers" [negative would be markers in other levels]

  blacklisted.rows <- grep(blacklist_pattern, rownames(ace), ignore.case = T)
  scores[blacklisted.rows, ] <- 0

  rowMaps(ace)[[sprintf("%s_feature_specificity_wilcox", out_name)]] <- scores
  rowMapTypes(ace)[[sprintf("%s_feature_specificity_wilcox", out_name)]] <- "reduction"

  return(ace)
}

computeGeneSpecifity.scran <- function(ace, f, out_name = "cond", pos_only = T, blacklist_pattern = "^MT-|^MT[:.:]|^RPS|^RPL|^MALAT") {
  require(scran)
  print("Running computeGeneSpecifity.scran()")

  if (class(f) != "factor") {
    warning("f must be a factor")
    return(ace)
  }
  f <- droplevels(f)

  out <- scran::findMarkers(logcounts(ace), f)
  metadata(ace)[[sprintf("%s_feature_specificity_scran", out_name)]] <- out

  scores <- do.call(cbind, lapply(out, function(DF) {
    idx <- match(rownames(ace), rownames(DF))
    x <- -log10(DF$p.value[idx]) * sign(DF$summary.logFC[idx])
  }))
  colnames(scores) <- levels(f)
  rownames(scores) <- rownames(ace)
  scores[scores == Inf] <- max(scores[!is.infinite(scores)])
  scores[scores == -Inf] <- min(scores[!is.infinite(scores)])

  if (pos_only == T) {
    scores[scores < 0] <- 0
  } # Only "positive markers" [negative would be markers in other levels]

  blacklisted.rows <- grep(blacklist_pattern, rownames(ace), ignore.case = T)
  scores[blacklisted.rows, ] <- 0

  rowMaps(ace)[[sprintf("%s_feature_specificity_scran", out_name)]] <- scores
  rowMapTypes(ace)[[sprintf("%s_feature_specificity_scran", out_name)]] <- "reduction"

  return(ace)
}

computeGeneSpecifity.limma <- function(ace, f, out_name = "cond", pos_only = T, blacklist_pattern = "^MT-|^MT[:.:]|^RPS|^RPL|^MALAT", weight = NULL) {
  require(limma)
  print("Running computeGeneSpecifity.limma()")

  if (class(f) != "factor") {
    warning("f must be a factor")
    return(ace)
  }
  f <- droplevels(f)

  design <- model.matrix(~ 0. + f)
  if (is.null(weight)) {
    fit0 <- limma::lmFit(logcounts(ace), design)
  } else {
    fit0 <- limma::lmFit(logcounts(ace), design, weights = weight)
  }
  fit <- eBayes(fit0, robust = TRUE, trend = TRUE)

  scores <- matrix(0, nrow = nrow(ace), ncol = length(levels(f)))
  colnames(scores) <- levels(f)
  rownames(scores) <- rownames(ace)
  DFs <- vector("list", length(levels(f)))
  names(DFs) <- levels(f)
  for (i in 1:length(levels(f))) {
    contrast.mat <- matrix(-1 / (ncol(design) - 1), nrow = ncol(design))
    contrast.mat[[i]] <- 1
    rownames(contrast.mat) <- colnames(design)

    fit <- limma::contrasts.fit(fit0, contrasts = contrast.mat)
    fit <- limma::eBayes(fit, robust = TRUE, trend = TRUE)
    DF <- limma::topTable(
      fit = fit,
      number = Inf,
      adjust.method = "BH",
      sort.by = "none"
    )

    DFs[[i]] <- DF
    scores[rownames(DF), i] <- DF$t
  }

  if (pos_only == T) {
    scores[scores < 0] <- 0
  } # Only "positive markers" [negative would be markers in other levels]

  blacklisted.rows <- grep(blacklist_pattern, rownames(ace), ignore.case = T)
  scores[blacklisted.rows, ] <- 0

  metadata(ace)[[sprintf("%s_feature_specificity_limma", out_name)]] <- DFs
  rowMaps(ace)[[sprintf("%s_feature_specificity_limma", out_name)]] <- scores
  rowMapTypes(ace)[[sprintf("%s_feature_specificity_limma", out_name)]] <- "reduction"

  return(ace)
}

computeGeneSpecifity.limma.pb <- function(ace, f, out_name = "cond", pos_only = T, blacklist_pattern = "^MT-|^MT[:.:]|^RPS|^RPL|^MALAT", resolution = 5, min_size = 10) {
  print("Running computeGeneSpecifity.limma.pb()")

  set.seed(0)
  if (resolution == -1) {
    cl <- ace$assigned_archetype
  } else {
    cl <- Leiden.clustering(ace, resolution_parameter = resolution)
  }

  f.prod <- interaction(cl, f)
  IDX <- split(1:ncol(ace), f.prod)
  mask <- sapply(IDX, length) > min_size
  S <- logcounts(ace)
  ll <- lapply(IDX[mask], function(idx) {
    return(fast_row_sums(S[, idx, drop = F]) / length(idx))
  })
  pb <- do.call(cbind, ll)
  weight <- sapply(IDX[mask], function(idx) {
    return(length(idx))
  })
  rownames(pb) <- rownames(ace)
  pb <- as.matrix(pb)

  pb.ace <- ACTIONetExperiment(assays = list(logcounts = pb))
  pb.ace$condition <- factor(sapply(colnames(pb.ace), function(str) stringr::str_split(str, stringr::fixed("."))[[1]][[2]]), levels(f))

  z <- (weight - median(weight)) / mad(weight)
  mask <- z > -1
  pb.ace <- computeGeneSpecifity.limma(pb.ace[, mask], pb.ace$condition[mask], out_name = out_name, pos_only = pos_only, blacklist_pattern = blacklist_pattern, weight = weight[mask])

  metadata(ace)[[sprintf("%s_feature_specificity_limma_pb", out_name)]] <- rowMaps(pb.ace)[[sprintf("%s_feature_specificity_limma", out_name)]]
  rowMaps(ace)[[sprintf("%s_feature_specificity_limma_pb", out_name)]] <- rowMaps(pb.ace)[[sprintf("%s_feature_specificity_limma", out_name)]]
  rowMapTypes(ace)[[sprintf("%s_feature_specificity_limmapb", out_name)]] <- "reduction"

  return(ace)
}

computeGeneSpecifity.ace <- function(ace, var_of_interest, method = "ACTIONet", out_name = "cond", pos_only = T, blacklist_pattern = "^MT-|^MT[:.:]|^RPS|^RPL|^MALAT", resolution = 5, min_size = 10, max_class = 100) {
  f <- process.var.of.interest(ace, var_of_interest, max_class = max_class)

  if (method == "ACTIONet") {
    ace <- computeGeneSpecifity.ACTIONet(ace, f, out_name = out_name, pos_only = pos_only, blacklist_pattern = blacklist_pattern)
  } else if (method == "wilcox") {
    ace <- computeGeneSpecifity.wilcox(ace, f, out_name = out_name, pos_only = pos_only, blacklist_pattern = blacklist_pattern)
  } else if (method == "scran") {
    ace <- computeGeneSpecifity.scran(ace, f, out_name = out_name, pos_only = pos_only, blacklist_pattern = blacklist_pattern)
  } else if (method == "limma") {
    ace <- computeGeneSpecifity.limma(ace, f, out_name = out_name, pos_only = pos_only, blacklist_pattern = blacklist_pattern)
  } else if (method == "limma_pseudobulk") {
    ace <- computeGeneSpecifity.limma.pb(ace, f, out_name = out_name, pos_only = pos_only, blacklist_pattern = blacklist_pattern, resolution = resolution, min_size = min_size)
  }

  return(ace)
}