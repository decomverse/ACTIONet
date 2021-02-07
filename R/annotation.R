
.preprocess_annotation_markers <- function(markers) {

    if (is.matrix(markers) || is.sparseMatrix(markers)) {
        marker_set = lapply(1:NCOL(markers), function(i) rownames(markers)[markers[, i] > 0])
        names(marker_set) = colnames(markers)
    } else if (is.data.frame(markers) || is(DataFrame(), "DFrame")) {
        marker_set = as.list(as.data.frame(markers))
    } else if (is.list(markers)) {
        marker_set = markers
    } else {
        err = sprintf("'markers' must be a list or matrix.\n")
        stop(err)
    }

    if (is.null(names(marker_set))) {
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
#' \item Label: Inferred archetype labels
#' \item Confidence: Confidence of inferred labels
#' \item Enrichment: Full enrichment matrix
#'}
#'
#' @examples
#' arch.annot = annotate.archetypes.using.labels(ace, sce$celltypes)
#' @export
annotate.archetypes.using.labels <- function(
  ace,
  labels,
  archetype.slot = "H_unified"
) {

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

        mu.class = fastRowMeans(class.profile)
        mu.null = fastRowMeans(null.profile)

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

    out = list(
      Label = archetypeLabels,
      Confidence = Labels.confidence,
      Enrichment = Enrichment.Z)

    return(out)
}


#' Annotate clusters using known marker genes
#' (It uses permutation test on cluster specificity scores)
#'
#' @param ace ACTIONet output object
#' @param markers A list of lists (each a set of markers for a given cell type)
#' @param rand_sample_no Number of random permutations (default=1000)
#'
#' @return A named list: \itemize{
#' \item Label: Inferred archetype labels
#' \item Confidence: Confidence of inferred labels
#' \item Enrichment: Full enrichment matrix
#'}
#'
#' @examples
#' data('curatedMarkers_human') # pre-packaged in ACTIONet
#' markers = curatedMarkers_human$Blood$PBMC$Monaco2019.12celltypes$marker.genes
#' arch.annot = annotate.archetypes.using.markers(ace, markers = markers)
#' @export
annotate.archetypes.using.markers <- function(
  ace,
  markers,
  rand_sample_no = 1000,
  significance_slot = "unified_feature_specificity"
) {

    marker_set = .preprocess_annotation_markers(markers)

    specificity.panel = Matrix::t(as.matrix(log1p(rowMaps(ace)[[significance_slot]])))
    specificity.panel[is.na(specificity.panel)] = 0

    if (is.null(names(marker_set))) {
        names(marker_set) = sapply(1:length(marker_set), function(i){
            msg = sprintf("Celltype %s", i)
            print(msg)
        })
    }

    markers.table = do.call(rbind, lapply(names(marker_set), function(celltype) {
        genes = marker_set[[celltype]]
        if (length(genes) == 0){
          err = sprintf("No markers left.\n")
          stop(err, call. = F)
        }

        signed.count = sum(sapply(genes, function(gene) grepl("\\+$|-$", gene)))
        is.signed = signed.count > 0

        if (!is.signed) {
            df = data.frame(
              Gene = genes,
              Direction = +1,
              Celltype = celltype,
              stringsAsFactors = FALSE
            )
        } else {

            pos.genes = (as.character(sapply(genes[grepl("+", genes, fixed = TRUE)],
                function(gene) stringr::str_replace(gene, stringr::fixed("+"), ""))))
            neg.genes = (as.character(sapply(genes[grepl("-", genes, fixed = TRUE)],
                function(gene) stringr::str_replace(gene, stringr::fixed("-"), ""))))

            df = data.frame(
              Gene = c(pos.genes, neg.genes),
              Direction = c(rep(+1, length(pos.genes)), rep(-1, length(neg.genes))),
              Celltype = celltype,
              stringsAsFactors = FALSE
            )
        }
    }))
    markers.table = markers.table[markers.table$Gene %in% colnames(specificity.panel),
        ]

    if (dim(markers.table)[1] == 0) {
      err = sprintf("No markers left.\n")
      stop(err, call. = F)
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

    out = list(
      Label = Labels,
      Confidence = Labels.conf,
      Enrichment = Z)

    return(out)
}

#' Infer cell annotations from imputed gene expression for all cells.
#'
#' @param ace ACTIONetExperiment object
#' @param markers A named list of marker genes.
#' @param features_use A vector of features of length NROW(ace) or the name of a column of rowData(ace) containing the genes given in 'markers'.
#' @param alpha_val Random-walk parameter for gene imputation.
#' @param thread_no Number of parallel threads used for gene imputation.
#' @param net_slot Name of slot in colNets(ace) containing the network to use for gene expression imputation (default="ACTIONet").
#' @param assay_name Name of assay for which to impute gene expression (default="logcounts").
#' @return A named list: \itemize{
#' \item Label: Inferred cell type labels
#' \item Confidence: Confidence of inferred labels
#' \item Enrichment: Cell type score matrix.
#'}
#'
#' @examples
#' data('curatedMarkers_human') # pre-packaged in ACTIONet
#' markers = curatedMarkers_human$Blood$PBMC$Monaco2019.12celltypes$marker.genes
#' annots = annotate.cells.using.markers(ace, markers = markers)
#' plot.ACTIONet(ace, annots$Label, annots$Confidence)
#' @export
annotate.cells.using.markers <- function(
  ace,
  markers,
  features_use = NULL,
  alpha_val = 0.9,
  thread_no = 0,
  net_slot = "ACTIONet",
  assay_name = "logcounts",
  max_iter = 5
) {

    marker_set = .preprocess_annotation_markers(markers)
    features_use = .preprocess_annotation_features(ace, features_use)

    marker_mat = as(sapply(marker_set, function(gs) as.numeric(features_use %in% gs) ), "sparseMatrix")
    G = colNets(ace)[[net_slot]]
    S = SummarizedExperiment::assays(ace)[[assay_name]]

    marker_stats = compute_marker_aggregate_stats(
      G = G,
      S = S,
      marker_mat = marker_mat,
      alpha = alpha_val,
      max_it = max_iter,
      thread_no = thread_no
    )
    colnames(marker_stats) = names(marker_set)

    marker_stats[!is.finite(marker_stats)] = 0
    annots = names(marker_set)[apply(marker_stats, 1, which.max)]
    conf = apply(marker_stats, 1, max)

    out = list(
      Label = annots,
      Confidence = conf,
      Enrichment = marker_stats
    )

    return(out)
}


#' Annotates cells by interpolating marker-based archetype annotations
#'
#' @param ace ACTIONet output object
#' @param markers A list of lists (each a set of markers for a given cell type)
#' @param rand_sample_no Number of random permutations (default=1000)
#'
#' @return A named list: \itemize{
#' \item Label: Inferred archetype labels
#' \item Confidence: Confidence of inferred labels
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
annotate.cells.from.archetypes.using.markers <- function(
  ace,
  markers,
  rand_sample_no = 1000,
  unified_suffix = "unified"
) {

    marker_set = markers
    significance_slot = sprintf("%s_feature_specificity", unified_suffix)
    arch.annot = annotate.archetypes.using.markers(
      ace = ace,
      markers = marker_set,
      rand_sample_no = rand_sample_no,
      significance_slot = significance_slot
    )

    enrichment.mat = arch.annot$Enrichment

    H.slot = sprintf("H_%s", unified_suffix)
    cell.enrichment.mat = map.cell.scores.from.archetype.enrichment(
      ace = ace,
      enrichment_mat = enrichment.mat,
      normalize = TRUE,
      H.slot = H.slot
    )
    cell.annotations = colnames(cell.enrichment.mat)[apply(cell.enrichment.mat, 1,
        which.max)]

    Labels = colnames(cell.enrichment.mat)[apply(cell.enrichment.mat, 1, which.max)]
    Labels.confidence = apply(cell.enrichment.mat, 1, max)

    res = list(
      Label = Labels,
      Confidence = Labels.confidence,
      Enrichment = cell.enrichment.mat
    )

    return(res)
}

#' Interpolates cell scores from archetype enrichment matrix
#'
#' @param ace ACTIONet output object
#' @param enrichment_mat Enrichment matrix with rows corresponding to archetypes and columns to an arbitrary annotation
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
map.cell.scores.from.archetype.enrichment <- function(
  ace,
  enrichment_mat,
  normalize = FALSE,
  H.slot = "H_unified"
) {

    cell.scores.mat = colMaps(ace)[[H.slot]]

    if (nrow(enrichment_mat) != ncol(cell.scores.mat)) {
        print("Flipping enrichment matrix")
        enrichment_mat = Matrix::t(enrichment_mat)
    }

    if (normalize == TRUE) {
        enrichment.scaled = doubleNorm(enrichment_mat)
    } else {
        enrichment.scaled = enrichment_mat
        enrichment.scaled[enrichment.scaled < 0] = 0
        if (max(enrichment.scaled) > 50) {
            enrichment.scaled = log1p(enrichment.scaled)
        }
    }

    cell.enrichment.mat = cell.scores.mat %*% enrichment.scaled
    colnames(cell.enrichment.mat) = colnames(enrichment_mat)
    rownames(cell.enrichment.mat) = colnames(ace)

    return(cell.enrichment.mat)
}
