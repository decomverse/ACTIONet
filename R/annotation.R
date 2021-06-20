.preprocess_annotation_markers <- function(markers, feature_set) {
    if (is.matrix(markers) || is.sparseMatrix(markers)) {
		common_features = sort(unique(intersect(feature_set, rownames(markers))))
		row_idx = match(common_features, rownames(markers))
		X = markers[row_idx, ]
		rownames(X) = common_features
    } else if (is.list(markers)) {
		if (is.null(names(markers))) {
			names(marker_set) = sapply(1:length(marker_set), function(i){
				msg = sprintf("Annotation %d", i)
				print(msg)
			})
		}		
		X = do.call(cbind, lapply(markers, function(gs) {
			genes = unlist(strsplit(gs, "[+]|[-]"))
			if(length(grep("[+|-]$", gs)) != 0) {
				x = as.numeric(grepl("[+]", gs)) - as.numeric(grepl("[-]", gs))
				genes = as.character(sapply(gs, function(gene) substr(gene, 1, stringr::str_length(gene)-1)))
			} else {
				x = rep(1, length(genes))
				genes = gs
			}
			common.genes = sort(unique(intersect(genes, feature_set)))
			idx1 = match(common.genes, genes)
			idx2 = match(common.genes, feature_set)

			v = sparseMatrix(i = idx2, j = rep(1, length(idx2)), x = x[idx1], dims = c(length(feature_set), 1))
			return(v)
		}))
		colnames(X) = names(markers)
		rownames(X) = feature_set
	} else if (is.data.frame(markers) || (length(which(is(markers) == "DFrame")) != 0)) { # marker, cell type, [weight]
		if(ncol(markers) == 2) {
			markers$weight = 1
		}
		UL = sort(unique(markers[, 2]))
		X = do.call(cbind, lapply(UL, function(nn) {
			idx = which(markers[, 2]== nn)
			genes = markers[idx, 1]
			x = markers[idx, 3]

			common.genes = sort(unique(intersect(genes, feature_set)))
			idx1 = match(common.genes, genes)
			idx2 = match(common.genes, feature_set)

			v = sparseMatrix(i = idx2, j = rep(1, length(idx2)), x = x[idx1], dims = c(length(feature_set), 1))
			return(v)

		}))
		colnames(X) = UL
		rownames(X) = feature_set
	}


	# cs = fast_column_sums(X)
	# filtered_cols = which(cs < 0)
	# if(length(filtered_cols)) {
	# 	not
	# }
	# col.mask = cs != 0
	# X.sp = as(X, 'dgTMatrix')
	# X.sp@x = X.sp@x / cs[X.sp@j+1]
	# X = as(X.sp, "dgCMatrix")

    return(X)
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
  features_use = NULL,
  significance_slot = "unified_feature_specificity") {


  features_use = .preprocess_annotation_features(ace, features_use)
  marker_mat = .preprocess_annotation_markers(markers, features_use)

  marker_stats = Matrix::t(assess.geneset.enrichment.from.archetypes(ace, marker_mat)$logPvals)
  colnames(marker_stats) = colnames(marker_mat)

  marker_stats[!is.finite(marker_stats)] = 0
  annots = colnames(marker_mat)[apply(marker_stats, 1, which.max)]
  conf = apply(marker_stats, 1, max)

  out = list(
    Label = annots,
    Confidence = conf,
    Enrichment = marker_stats
  )

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

  features_use = .preprocess_annotation_features(ace, features_use)
  marker_mat = .preprocess_annotation_markers(markers, features_use)

  #marker_mat = as(sapply(marker_set, function(gs) as.numeric(features_use %in% gs) ), "sparseMatrix")
  G = colNets(ace)[[net_slot]]
  S = as(SummarizedExperiment::assays(ace)[[assay_name]], "sparseMatrix")

  marker_stats = compute_marker_aggregate_stats(
    G = G,
    S = S,
    marker_mat = marker_mat,
    alpha = alpha_val,
    max_it = max_iter,
    thread_no = thread_no
  )
  colnames(marker_stats) = colnames(marker_mat)

  marker_stats[!is.finite(marker_stats)] = 0
  annots = colnames(marker_mat)[apply(marker_stats, 1, which.max)]
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
  unified_suffix = "unified"
) {

    marker_set = markers
    significance_slot = sprintf("%s_feature_specificity", unified_suffix)
    arch.annot = annotate.archetypes.using.markers(
      ace = ace,
      markers = marker_set,
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
