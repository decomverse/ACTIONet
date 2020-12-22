
.preprocess_annotation_markers <- function(markers){
  if (is.matrix(markers) | is.sparseMatrix(markers)) {
      marker_set = lapply(1:NCOL(markers), function(i) rownames(markers)[markers[,i] > 0])
      names(marker_set) = colnames(markers)
  } else if(is.list(markers)){
      marker_set = markers
  } else{
      err = sprintf("'markers' must be a list or matrix.\n")
      stop(err)
  }

  if(is.null(names(marker_set))) {
      names(marker_set) = sapply(1:length(marker_set), function(i) sprintf("Celltype %s", i))
  }
  return(marker_set)
}




#' Annotate archetypes using prior cell annotations
#' (It uses t-test on the archetype footprint matrix (H))
#'
#' @param ace ACTIONet output object
#' @param labels Annotation of interest (clusters, celltypes, etc.) to test enrichment
#'
#' @return A named list: \itemize{
#' \item Labels: Inferred archetype labels
#' \item Labels.confidence: Confidence of inferred labels
#' \item Enrichment: Full enrichment matrix
#'}
#'
#' @examples
#' arch.annot = annotate.archetypes.using.labels(ace, sce$celltypes)
#' @export
annotate.archetypes.using.labels <- function(ace, labels, archetype.slot = "H_unified") {
    Labels = .preprocess_annotation_labels(labels, ace)

    if (is.matrix(ace) | is.sparseMatrix(ace)) {
        profile = as.matrix(ace)
    } else {
        profile = Matrix::t(colMaps(ace)[[archetype.slot]])
    }
    Annot = names(Labels)[match(sort(unique(Labels)), Labels)]

    # Using t-statistics
    Enrichment.Z = sapply(Annot, function(label) {
        mask = names(Labels) == label
        class.profile = profile[, mask]
        null.profile = profile[, !mask]

        N.class = sum(mask)
        N.null = sum(!mask)

        if ((N.class < 3) | (N.null < 3)) {
            return(rep(0, nrow(profile)))
        }

        mu.class = Matrix::rowMeans(class.profile)
        mu.null = Matrix::rowMeans(null.profile)

        sigma_sq.class = apply(class.profile, 1, var)
        sigma_sq.null = apply(null.profile, 1, var)


        delta.mean = mu.class - mu.null
        t.stat = delta.mean/sqrt((sigma_sq.class/N.class) + (sigma_sq.null/N.null))
        return(t.stat)
    })
	Enrichment.Z[is.na(Enrichment.Z)] = 0

    archetypeLabels = Annot[apply(Enrichment.Z, 1, which.max)]
    Labels.confidence = apply(Enrichment.Z, 1, max)


    rownames(Enrichment.Z) = paste("A", 1:nrow(Enrichment.Z), "-", archetypeLabels,
        sep = "")

    out.list = list(Labels = archetypeLabels, Labels.confidence = Labels.confidence,
        Enrichment = Enrichment.Z)

    return(out.list)
}


#' Annotate clusters using known marker genes
#' (It uses permutation test on cluster specificity scores)
#'
#' @param ace ACTIONet output object
#' @param markers A list of lists (each a set of markers for a given cell type)
#' @param rand_sample_no Number of random permutations (default=1000)
#'
#' @return A named list: \itemize{
#' \item Labels: Inferred archetype labels
#' \item Labels.confidence: Confidence of inferred labels
#' \item Enrichment: Full enrichment matrix
#'}
#'
#' @examples
#' data('curatedMarkers_human') # pre-packaged in ACTIONet
#' markers = curatedMarkers_human$Blood$PBMC$Monaco2019.12celltypes$marker.genes
#' arch.annot = annotate.archetypes.using.markers(ace, markers = markers)
#' @export
annotate.archetypes.using.markers <- function(ace, markers, rand_sample_no = 1000,
    significance.slot = "unified_feature_specificity") {

    marker_set = .preprocess_annotation_markers(markers)
    # features_use = .preprocess_annotation_features(ace, features_use)

    specificity.panel = Matrix::t(as.matrix(log1p(rowMaps(ace)[[significance.slot]])))
    specificity.panel[is.na(specificity.panel)] = 0

    if(is.null(names(marker_set))) {
		names(marker_set) = sapply(1:length(marker_set), function(i) sprintf("Celltype %s", i))
	}

    markers.table = do.call(rbind, lapply(names(marker_set), function(celltype) {
        genes = marker_set[[celltype]]
        if (length(genes) == 0)
            return(data.frame())


        signed.count = sum(sapply(genes, function(gene) grepl("\\+$|-$", gene)))
        is.signed = signed.count > 0

        if (!is.signed) {
            df = data.frame(Gene = (genes), Direction = +1, Celltype = celltype,
                stringsAsFactors = F)
        } else {

            pos.genes = (as.character(sapply(genes[grepl("+", genes, fixed = TRUE)],
                function(gene) stringr::str_replace(gene, stringr::fixed("+"), ""))))
            neg.genes = (as.character(sapply(genes[grepl("-", genes, fixed = TRUE)],
                function(gene) stringr::str_replace(gene, stringr::fixed("-"), ""))))

            df = data.frame(Gene = c(pos.genes, neg.genes), Direction = c(rep(+1,
                length(pos.genes)), rep(-1, length(neg.genes))), Celltype = celltype,
                stringsAsFactors = F)
        }
    }))
    markers.table = markers.table[markers.table$Gene %in% colnames(specificity.panel), ]

    if (dim(markers.table)[1] == 0) {
        print("No markers are left")
        return()
    }
    specificity.panel = specificity.panel[, markers.table$Gene]

    IDX = split(1:dim(markers.table)[1], markers.table$Celltype)

    print("Computing significance scores")
    set.seed(0)
    Z = sapply(IDX, function(idx) {
        markers = (as.character(markers.table$Gene[idx]))
        directions = markers.table$Direction[idx]
        mask = markers %in% colnames(specificity.panel)

        A = as.matrix(specificity.panel[, markers[mask]])
        sgn = as.numeric(directions[mask])
        stat = A %*% sgn

        rand.stats = sapply(1:rand_sample_no, function(i) {
            rand.samples = sample.int(dim(specificity.panel)[2], sum(mask))
            rand.A = as.matrix(specificity.panel[, rand.samples])
            rand.stat = rand.A %*% sgn
        })
        # rand.stats[is.na(rand.stats)] = 0
        cell.zscores = as.numeric((stat - apply(rand.stats, 1, mean))/apply(rand.stats,
            1, sd))

        return(cell.zscores)
    })

    Z[is.na(Z)] = 0
    Labels = colnames(Z)[apply(Z, 1, which.max)]

    # L = names(marker_set) L.levels = L[L %in% Labels] Labels = match(L, L.levels)
    # names(Labels) = L.levels Labels = factor(Labels, levels = L)
    Labels.conf = apply(Z, 1, max)

    names(Labels) = rownames(specificity.panel)
    names(Labels.conf) = rownames(specificity.panel)
    rownames(Z) = rownames(specificity.panel)

    out.list = list(Labels = Labels, Labels.confidence = Labels.conf, Enrichment = Z)

    return(out.list)
}

#' Directly inferring cell annotations from imputed gene expressions using permutation test
#'
#' @param ace ACTIONet output object
#' @param markers A list of lists (each a set of markers for a given cell type)
#' @param rand_sample_no Number of random permutations (default=1000)
#' @param alpha_val Random-walk parameter for gene imputation (if imputation = 'PageRank')
#' @param imputation Gene imputation method. Can be either 'PageRank' (default) or 'archImpute'
#' @param thread_no Number of parallel threads used for gene imputation (if imputation = 'PageRank')
#'
#' @return A named list: \itemize{
#' \item Labels: Inferred archetype labels
#' \item Labels.confidence: Confidence of inferred labels
#' \item Enrichment: Full enrichment matrix
#'}
#'
#' @examples
#' data('curatedMarkers_human') # pre-packaged in ACTIONet
#' markers = curatedMarkers_human$Blood$PBMC$Monaco2019.12celltypes$marker.genes
#' arch.annot = annotate.cells.using.markers(ace, markers = markers)
#' cell.labels = arch.annot$Labels
#' @export
annotate.cells.using.markers <- function(ace, markers, features_use = NULL, rand_sample_no = 100,
    alpha_val = 0.9, thread_no = 8, imputation = "PageRank", assay_name = "logcounts", diffusion_iters = 5) {

    marker_set = .preprocess_annotation_markers(markers)
    features_use = .preprocess_annotation_features(ace, features_use)

    markers.table = do.call(rbind, lapply(names(marker_set), function(celltype) {
        genes = marker_set[[celltype]]
        if (length(genes) == 0)
            return(data.frame())


        signed.count = sum(sapply(genes, function(gene) grepl("\\+$|-$", gene)))
        is.signed = signed.count > 0

        if (!is.signed) {
            df = data.frame(Gene = (genes), Direction = +1, Celltype = celltype)
        } else {

            pos.genes = (as.character(sapply(genes[grepl("+", genes, fixed = TRUE)],
                function(gene) stringr::str_replace(gene, stringr::fixed("+"), ""))))
            neg.genes = (as.character(sapply(genes[grepl("-", genes, fixed = TRUE)],
                function(gene) stringr::str_replace(gene, stringr::fixed("-"), ""))))

            df = data.frame(Gene = c(pos.genes, neg.genes), Direction = c(rep(+1,
                length(pos.genes)), rep(-1, length(neg.genes))), Celltype = celltype)
        }
    }))
    markers.table = markers.table[markers.table$Gene %in% features_use, ]

    if (dim(markers.table)[1] == 0) {
        err = sprintf("No given markers found in feature set.\n")
        stop(err, call. = F)
    }

    rows = match(markers.table$Gene, features_use)
    if (imputation == "PageRank") {
        # PageRank-based imputation
        cat("Imputing marker expression using PageRank ... ")
        expression_profile = impute.genes.using.ACTIONet(ace, markers.table$Gene, features_use = features_use, alpha_val = alpha_val, thread_no = thread_no, assay_name = assay_name, diffusion_iters = diffusion_iters)
        cat(sprintf("done.\n"))
    } else if (imputation == "archImpute") {
        # PCA-based imputation
        cat("Imputing marker expression using archImpute ... ")
        expression_profile = impute.specific.genes.using.archetypes(ace, markers.table$Gene)
        cat(sprintf("done.\n"))
    } else {
        expression_profile = SummarizedExperiment::assays(ace)[[assay_name]]
    }

    IDX = split(1:dim(markers.table)[1], markers.table$Celltype)

    set.seed(0)
    cat("Computing significance scores ... ")
    Z = sapply(IDX, function(idx) {
        markers = (as.character(markers.table$Gene[idx]))
        directions = markers.table$Direction[idx]
        mask = markers %in% colnames(expression_profile)

        A = as.matrix(expression_profile[, markers[mask]])
        sgn = as.numeric(directions[mask])
        stat = A %*% sgn

        rand.stats = sapply(1:rand_sample_no, function(i) {
            rand.samples = sample.int(dim(expression_profile)[2], sum(mask))
            rand.A = as.matrix(expression_profile[, rand.samples])
            rand.stat = rand.A %*% sgn
        })

        cell.zscores = as.numeric((stat - apply(rand.stats, 1, mean))/apply(rand.stats,
            1, sd))

        return(cell.zscores)
    })
    cat(sprintf("done.\n"))

    Z[is.na(Z)] = 0
    Labels = apply(Z, 1, which.max)
    Labels = colnames(Z)[Labels]
    Labels.conf = apply(Z, 1, max)

    res = list(Labels = Labels, Labels.confidence = Labels.conf, Enrichment = Z)

    return(res)
}


#' Annotates cells by interpolating marker-based archetype annotations
#'
#' @param ace ACTIONet output object
#' @param markers A list of lists (each a set of markers for a given cell type)
#' @param rand_sample_no Number of random permutations (default=1000)
#'
#' @return A named list: \itemize{
#' \item Labels: Inferred archetype labels
#' \item Labels.confidence: Confidence of inferred labels
#' \item Enrichment: Full enrichment matrix
#'}
#'
#' @examples
#'
#' data('curatedMarkers_human') # pre-packaged in ACTIONet
#' marker_set = curatedMarkers_human$Blood$PBMC$Monaco2019.12celltypes$marker.genes
#' cell.annotations = annotate.cells.from.archetypes.using.markers(ace, markers)
#' labels = cell.annotations$Labels
#' @export
annotate.cells.from.archetypes.using.markers <- function(ace, markers, rand_sample_no = 1000,
    unified_suffix = "unified") {

    marker_set = markers
    significance.slot = sprintf("%s_feature_specificity", unified_suffix)
    arch.annot = annotate.archetypes.using.markers(ace, markers = marker_set,
        rand_sample_no = rand_sample_no, significance.slot = significance.slot)

    enrichment.mat = arch.annot$Enrichment

    H.slot = sprintf("H_%s", unified_suffix)
    cell.enrichment.mat = map.cell.scores.from.archetype.enrichment(ace, enrichment.mat,
        normalize = T, H.slot = H.slot)
    cell.annotations = colnames(cell.enrichment.mat)[apply(cell.enrichment.mat, 1,
        which.max)]

    Labels = colnames(cell.enrichment.mat)[apply(cell.enrichment.mat, 1, which.max)]
    Labels.confidence = apply(cell.enrichment.mat, 1, max)

    res = list(Labels = Labels, Labels.confidence = Labels.confidence, Enrichment = cell.enrichment.mat)

    return(res)
}

#' Interpolates cell scores from archetype enrichment matrix
#'
#' @param ace ACTIONet output object
#' @param enrichment.matrix Enrichment matrix with rows corresponding to archetypes and columns to an arbitrary annotation
#' @param normalize If TRUE, enrichment matrix will be first doubly-normalized
#'
#' @return Enrichment map of size cell x annotation
#'
#' @examples
#'
#' data('curatedMarkers_human') # pre-packaged in ACTIONet
#' marker_set = curatedMarkers_human$Blood$PBMC$Monaco2019.12celltypes$marker.genes
#' arch.annot = annotate.archetypes.using.markers(ace, markers = markers)
#' enrichment.mat = arch.annot$enrichment
#' cell.enrichment.mat = map.cell.scores.from.archetype.enrichment(ace, enrichment.mat)
#' cell.assignments = colnames(cell.enrichment.mat)[apply(cell.enrichment.mat, 1, which.max)]
#' @export
map.cell.scores.from.archetype.enrichment <- function(ace, enrichment.matrix, normalize = F, H.slot = "H_unified") {

  cell.scores.mat = colMaps(ace)[[H.slot]]

  if (nrow(enrichment.matrix) != ncol(cell.scores.mat)) {
      print("Flipping enrichment matrix")
      enrichment.matrix = t(enrichment.matrix)
  }

  if (normalize == T) {
      enrichment.scaled = doubleNorm(enrichment.matrix)
  } else {
      enrichment.scaled = enrichment.matrix
      enrichment.scaled[enrichment.scaled < 0] = 0
      if (max(enrichment.scaled) > 50) {
          enrichment.scaled = log1p(enrichment.scaled)
      }
  }

  cell.enrichment.mat = cell.scores.mat %*% enrichment.scaled
  colnames(cell.enrichment.mat) = colnames(enrichment.matrix)
  rownames(cell.enrichment.mat) = colnames(ace)

  return(cell.enrichment.mat)
}
