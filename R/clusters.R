#' Computes feature (i.e. gene) specificity scores for each cluster
#'
#' @param ace ACTIONet output object
#' @param clusters Cluster
#' @param output_slot Name of the output in rowMaps(ace) to store results
#' @param renormalize.logcounts.slot Name of the new assay with updated logcounts adjusted using archetypes
#' Typically it is either 'logcounts' or 'logcounts'

#' @return `ACE` object with specificity scores of each cluster added to rowMaps(ace) as a matrix with name defined by output_slot
#'
#' @examples
#' ace = compute.cluster.feature.specificity(ace, ace$clusters, 'cluster_specificity_scores')
#' @export
compute.cluster.feature.specificity <- function(
  ace,
  clusters,
  output_slot,
  assay_name = "logcounts"
) {

    S = SummarizedExperiment::assays(ace)[[assay_name]]

    if (is.factor(clusters)) {
        UL = levels(clusters)
    } else {
        UL = sort(unique(clusters))
    }
    lables = match(clusters, UL)

    # Compute gene specificity for each cluster
    if (is.matrix(S)) {
        specificity.out = compute_cluster_feature_specificity_full(S, lables)
    } else {
        specificity.out = compute_cluster_feature_specificity(S, lables)
    }

    specificity.out = lapply(specificity.out, function(specificity.scores) {
        rownames(specificity.scores) = rownames(ace)
        colnames(specificity.scores) = paste("A", 1:ncol(specificity.scores))
        return(specificity.scores)
    })

    X = specificity.out[["upper_significance"]]
    colnames(X) = UL

    rowMaps(ace)[[sprintf("%s_feature_specificity", output_slot)]] = X
    rowMapTypes(ace)[[sprintf("%s_feature_specificity", output_slot)]] = "reduction"

    return(ace)
}


#' Annotate clusters using prior cell annotations
#' (It uses Fisher's exact test for computing overlaps -- approximate HGT is used)
#'
#' @param ace ACTIONet output object
#' @param clusters Cluster
#' @param labels Annotation of interest (clusters, celltypes, etc.) to test enrichment
#'
#' @return A named list: \itemize{
#' \item Labels: Inferred archetype labels
#' \item Labels.confidence: Confidence of inferred labels
#' \item Enrichment: Full enrichment matrix
#'}
#'
#' @examples
#' arch.annot = annotate.clusters.using.labels(ace, ace$clusters, sce$celltypes)
#' @export
annotate.clusters.using.labels <- function(
  ace,
  clusters,
  labels
) {

    clusters = .preprocess_annotation_labels(clusters, ace)
    Labels = .preprocess_annotation_labels(labels, ace)

    pop.size = length(Labels)
    pos.size = table(Labels)

    logPvals = sapply(sort(unique(clusters)), function(i) {
        idx = which(clusters == i)
        sample.size = length(idx)
        success.size = sapply(sort(unique(Labels)), function(i) {
            sum(Labels[idx] == i)
        })

        logPval = HGT_tail(
          population.size = pop.size,
          success.count = pos.size,
          sample.size = sample.size,
          observed.success = success.size
        )

        return(logPval)
    })


    cl.Annot = names(clusters)[match(sort(unique(clusters)), clusters)]
    Annot = names(Labels)[match(sort(unique(Labels)), Labels)]

    colnames(logPvals) = cl.Annot
    rownames(logPvals) = Annot

    clusterLabels = Annot[apply(logPvals, 2, which.max)]

    cellLabels = match(clusterLabels[clusters], Annot)
    names(cellLabels) = clusterLabels[clusters]

    res = list(
      Labels = clusterLabels,
      cellLabels = cellLabels,
      Enrichment = logPvals
    )

    return(res)
}


#' Annotate clusters using known marker genes
#' (It uses permutation test on cluster specificity scores)
#'
#' @param ace ACTIONet output object
#' @param marker.genes A list of lists (each a set of markers for a given cell type)
#' @param specificity.slot.name An entry in the rowMaps(ace), precomputed using compute.cluster.feature.specificity() function
#' @param rand.sample.no Number of random permutations (default=1000)
#'
#' @return A named list: \itemize{
#' \item Labels: Inferred archetype labels
#' \item Labels.confidence: Confidence of inferred labels
#' \item Enrichment: Full enrichment matrix
#'}
#'
#' @examples
#' data('curatedMarkers_human') # pre-packaged in ACTIONet
#' marker.genes = curatedMarkers_human$Blood$PBMC$Monaco2019.12celltypes$marker.genes
#' ace = compute.cluster.feature.specificity(ace, ace$clusters, 'cluster_specificity_scores')
#' arch.annot = annotate.clusters.using.markers(ace, marker.genes = marker.genes, specificity.slot.name = 'cluster_specificity_scores')
#' @export
annotate.clusters.using.markers <- function(
  ace,
  marker.genes,
  specificity.slot.name,
  rand.sample.no = 1000
) {

    if (!grepl("_feature_specificity", specificity.slot.name)) {
        specificity.slot.name = paste(specificity.slot.name, "feature_specificity", sep = "_")
    }
    if (!(specificity.slot.name %in% names(rowMaps(ace)))) {
        message(sprintf("%s does not exist in rowMaps(ace)", specificity.slot.name))
    }

    if (is.matrix(marker.genes) | is.sparseMatrix(marker.genes)) {
        marker.genes = apply(marker.genes, 2, function(x) rownames(marker.genes)[x >
            0])
    }

    specificity.panel = Matrix::t(as.matrix(log1p(rowMaps(ace)[[specificity.slot.name]])))

    GS.names = names(marker.genes)
    if (is.null(GS.names)) {
        GS.names = sapply(1:length(GS.names), function(i) sprintf("Celltype %s", i))
    }

    markers.table = do.call(rbind, lapply(names(marker.genes), function(celltype) {
        genes = marker.genes[[celltype]]

        if (length(genes) == 0){
          err = sprintf("No markers left.\n")
          stop(err, call. = FALSE)
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

    markers.table = markers.table[markers.table$Gene %in% colnames(specificity.panel),]

    if (dim(markers.table)[1] == 0) {
      err = sprintf("No markers left.\n")
      stop(err, call. = FALSE)
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

        rand.stats = sapply(1:rand.sample.no, function(i) {
            rand.samples = sample.int(dim(specificity.panel)[2], sum(mask))
            rand.A = as.matrix(specificity.panel[, rand.samples])
            rand.stat = rand.A %*% sgn
        })

        cell.zscores = as.numeric((stat - apply(rand.stats, 1, mean))/apply(rand.stats,
            1, sd))

        return(cell.zscores)
    })

    Z[is.na(Z)] = 0
    Labels = colnames(Z)[apply(Z, 1, which.max)]

    Labels.conf = apply(Z, 1, max)

    names(Labels) = rownames(specificity.panel)
    names(Labels.conf) = rownames(specificity.panel)
    rownames(Z) = rownames(specificity.panel)

    out = list(
      Label = Labels,
      Confidence = Labels.conf,
      Enrichment = Z
    )

    return(out)
}



#' Annotate arbitary feature score matrix using known marker genes
#' (It uses permutation test on cluster specificity scores)
#'
#' @param feature.scores An arbitrary matrix with rows corresponding to features and columns to any given annotation/grouping of cells
#' @param marker.genes A list of lists (each a set of markers for a given cell type)
#' @param rand.sample.no Number of random permutations (default=1000)
#'
#' @return A named list: \itemize{
#' \item Labels: Inferred archetype labels
#' \item Labels.confidence: Confidence of inferred labels
#' \item Enrichment: Full enrichment matrix
#'}
#'
#' @examples
#' data('curatedMarkers_human') # pre-packaged in ACTIONet
#' marker.genes = curatedMarkers_human$Blood$PBMC$Monaco2019.12celltypes$marker.genes
#' arch.annot = annotate.profile.using.markers(my.gene.scores.profile, marker.genes = marker.genes)
#' @export
annotate.profile.using.markers <- function(
  feature.scores,
  marker.genes,
  rand.sample.no = 1000
) {

    if (is.matrix(marker.genes) | is.sparseMatrix(marker.genes)) {
        marker.genes = apply(marker.genes, 2, function(x) rownames(marker.genes)[x > 0])
    }

    specificity.panel = feature.scores

    GS.names = names(marker.genes)
    if (is.null(GS.names)) {
        GS.names = sapply(1:length(GS.names), function(i) sprintf("Celltype %s", i))
    }

    markers.table = do.call(rbind, lapply(names(marker.genes), function(celltype) {
        genes = marker.genes[[celltype]]

        if (length(genes) == 0) {
          err = sprintf("No markers left.\n")
          stop(err, call. = FALSE)
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

    markers.table = markers.table[markers.table$Gene %in% colnames(specificity.panel),]

    if (dim(markers.table)[1] == 0) {
      err = sprintf("No markers left.\n")
      stop(err, call. = FALSE)
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

        rand.stats = sapply(1:rand.sample.no, function(i) {
            rand.samples = sample.int(dim(specificity.panel)[2], sum(mask))
            rand.A = as.matrix(specificity.panel[, rand.samples])
            rand.stat = rand.A %*% sgn
        })

        cell.zscores = as.numeric((stat - apply(rand.stats, 1, mean))/apply(rand.stats, 1, sd))

        return(cell.zscores)
    })

    Z[is.na(Z)] = 0
    Labels = colnames(Z)[apply(Z, 1, which.max)]

    Labels.conf = apply(Z, 1, max)

    names(Labels) = rownames(specificity.panel)
    names(Labels.conf) = rownames(specificity.panel)
    rownames(Z) = rownames(specificity.panel)

    out = list(
      Label = Labels,
      Confidence = Labels.conf,
      Enrichment = Z
    )

    return(out)
}

#' A wrapper function For Leiden algorithm
#'
#' @param G Adjacency matrix of the input graph
#' @param resolution_parameter Resolution of the clustering.
#' The higher the resolution, the more clusters we will get (default=0.5).
#' @param initial.clustering Used as the inital clustering
#' @param seed Random seed
#'
#' @return ace with added annotation
#'
#' @examples
#' clusters = cluster.graph(G, 1.0)
#' @export
cluster.graph <- function(
  G,
  resolution_parameter = 0.5,
  initial.clustering = NULL,
  seed = 0
) {

    if (is.matrix(G)) {
        G = as(G, "sparseMatrix")
    }

    is.signed = FALSE
    if (min(G) < 0) {
        is.signed = TRUE
        print("Graph is signed. Switching to signed graph clustering mode.")
    }

    if (!is.null(initial.clustering)) {
        print("Perform graph clustering with *prior* initialization")

        if (is.signed) {
            clusters = as.numeric(signed_cluster(
              A = G,
              resolution_parameter = resolution_parameter,
              initial_clusters_ = initial.clustering,
              seed = seed
            ))
        } else {
            clusters = as.numeric(unsigned_cluster(
              A = G,
              resolution_parameter = resolution_parameter,
              initial_clusters_ = initial.clustering,
              seed = seed
            ))
        }
    } else {
        print("Perform graph clustering with *uniform* initialization")

        if (is.signed) {
            clusters = as.numeric(signed_cluster(
              A = G,
              resolution_parameter = resolution_parameter,
              initial_clusters_ = NULL,
              seed = seed
            ))
        } else {
            clusters = as.numeric(unsigned_cluster(
              A = G,
              resolution_parameter = resolution_parameter,
              initial_clusters_ = NULL,
              seed = seed
            ))
        }
    }

}


#' A wrapper function For Leiden algorithm applied to an ACE object
#'
#' @param ace Input results to be clustered

#' @param resolution_parameter Resolution of the clustering.
#' The higher the resolution, the more clusters we will get (default=0.5).
#' @param arch.init Whether to use archetype-assignments to initialize clustering (default=TRUE)
#' @param seed Random seed
#'
#' @return clusters
#'
#' @examples
#' clusters = Leiden.clustering(ace)
#' plot.ACTIONet(ace, clusters)
#' @export
Leiden.clustering <- function(
  ace,
  resolution_parameter = 1,
  net.slot = "ACTIONet",
  init.slot = "assigned_archetype",
  seed = 0
) {

    initial.clusters = NULL
    if (!is.null(init.slot)) {
        initial.clusters = ace[[init.slot]]
    }

    G = colNets(ace)[[net.slot]]

    clusters = cluster.graph(G, resolution_parameter, initial.clusters, seed)
    names(clusters) = paste("C", as.character(clusters), sep = "")

    return(clusters)
}


#' A wrapper function For HDBSCAN algorithm applied to an ACE object
#'
#' @param ace Input results to be clustered
#' minPoints, minClusterSize HDBSCAN parameters (default = 30,30)
#' archetype.slot Slot of archeypte to use for clustering (default='H_unified');
#'
#' @return clusters
#'
#' @examples
#' clusters = HDBSCAN.clustering(ace)
#' plot.ACTIONet(ace, clusters)
#' @export
HDBSCAN.clustering <- function(
  ace,
  minPoints = 30,
  minClusterSize = 30,
  archetype.slot = "H_unified"
) {

    X = as.matrix(colMaps(ace)[[archetype.slot]])
    
    out_list = run_HDBSCAN(
      X = X,
      minPoints = minPoints,
      minClusterSize = minClusterSize)

    return(out_list)
}
