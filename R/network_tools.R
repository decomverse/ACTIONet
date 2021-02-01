#' Uses a variant of the label propagation algorithm to infer missing labels
#'
#' @param ace Input results to be clustered
#' (alternatively it can be the ACTIONet igraph object)
#' @param initial_labels Annotations to correct with missing values (NA) in it.
#' It can be either a named annotation (inside ace$annotations) or a label vector.
#' @param double.stochastic Whether to densify adjacency matrix before running label propagation (default=FALSE).
#' @param max_iter How many iterative rounds of correction/inference should be performed (default=3)
#' @param adjust.levels Whether or not re-adjust labels at the end
#'
#' @return ace with updated annotations added to ace$annotations
#'
#' @examples
#' ace = infer.missing.cell.annotations(ace, sce$assigned_archetypes, 'updated_archetype_annotations')
#' @export
infer.missing.cell.annotations <- function(ace, initial_labels, double.stochastic = FALSE, 
    max_iter = 3, adjust.levels = TRUE) {
    
    Adj = colNets(ace)$ACTIONet
    A = as(Adj, "dgTMatrix")
    
    eps = 1e-16
    rs = fastRowSums(A)
    P = Matrix::sparseMatrix(i = A@i + 1, j = A@j + 1, x = A@x/rs[A@i + 1], dims = dim(A))
    
    if (double.stochastic == TRUE) {
        w = sqrt(fastColSums(P) + eps)
        W = P %*% Matrix::Diagonal(x = 1/w, n = length(w))
        P = W %*% Matrix::t(W)
    }
    
    Labels = .preprocess_annotation_labels(initial_labels, ace)
    if (is.null(Labels)) {
        return(ace)
    }
    
    Annot = sort(unique(Labels))
    idx = match(Annot, Labels)
    names(Annot) = names(Labels)[idx]
    
    na.mask = is.na(Labels)
    
    i = 1
    while (sum(na.mask) > 0) {
        new.Labels = assess.label.local.enrichment(P, Labels)
        
        mask = na.mask & (new.Labels$Labels.confidence > 3 + log(length(Labels)))
        Labels[mask] = new.Labels$Labels[mask]
        
        na.mask = is.na(Labels)
        if (i == max_iter) 
            break
        
        i = i + 1
    }
    new.Labels = assess.label.local.enrichment(P, Labels)
    Labels[na.mask] = new.Labels$Labels[na.mask]
    Labels.conf = new.Labels$Labels.confidence
    
    updated.Labels = as.numeric(Labels)
    names(updated.Labels) = names(Annot)[match(Labels, Annot)]
    if (adjust.levels == T) {
        Labels = reannotate.labels(ace, updated.Labels)
    } else {
        Labels = updated.Labels
    }
    
    return(Labels)
}

#' Uses a variant of the label propagation algorithm to correct likely noisy labels
#'
#' @param ace Input results to be clustered
#' (alternatively it can be the ACTIONet igraph object)
#' @param initial_labels Annotations to correct with missing values (NA) in it.
#' It can be either a named annotation (inside ace$annotations) or a label vector.
#' @param LFR.threshold How aggressively to update labels. The smaller the value, the more labels will be changed (default=2)
#' @param double.stochastic Whether to densify adjacency matrix before running label propagation (default=FALSE).
#' @param max_iter How many iterative rounds of correction/inference should be performed (default=3)
#' @param min.cell.fraction Annotations with less that this fraction will be removed
#'
#' @return ace with updated annotations added to ace$annotations
#'
#' @examples
#' ace = add.cell.annotations(ace, cell.labels, 'input_annotations')
#' ace = correct.cell.annotations(ace, 'input_annotations', 'updated_annotations')
#' @export
correct.cell.annotations <- function(ace, initial_labels, LFR.threshold = 2, double.stochastic = FALSE, 
    max_iter = 3, adjust.levels = TRUE, min.cell.fraction = 0.001) {
    
    Adj = colNets(ace)$ACTIONet
    A = as(Adj, "dgTMatrix")
    
    eps = 1e-16
    rs = fastRowSums(A)
    P = Matrix::sparseMatrix(i = A@i + 1, j = A@j + 1, x = A@x/rs[A@i + 1], dims = dim(A))
    
    if (double.stochastic == TRUE) {
        w = sqrt(fastColSums(P) + eps)
        W = P %*% Matrix::Diagonal(x = 1/w, n = length(w))
        P = W %*% Matrix::t(W)
    }
    
    Labels = .preprocess_annotation_labels(initial_labels, ace)
    if (is.null(Labels)) {
        return(ace)
    }
    
    # Prunes 'trivial' annotations and merges them to larger ones
    min.cells = round(min.cell.fraction * length(Labels))
    counts = table(Labels)
    Labels[Labels %in% as.numeric(names(counts)[counts < min.cells])] = NA
    mask = is.na(Labels)
    
    if (sum(mask) > 0) {
        Labels[mask] = NA
        Labels = infer.missing.cell.annotations(ace, initial_labels = Labels)
    }
    
    Annot = sort(unique(Labels))
    idx = match(Annot, Labels)
    names(Annot) = names(Labels)[idx]
    
    for (i in 1:max_iter) {
        new.Labels = assess.label.local.enrichment(P, Labels)
        Enrichment = new.Labels$Enrichment
        curr.enrichment = sapply(1:nrow(Enrichment), function(k) Enrichment[k, names(Labels)[k]])
        
        Diff.LFR = log2((new.Labels$Labels.confidence/curr.enrichment))
        Diff.LFR[is.na(Diff.LFR)] = 0
        
        Labels[Diff.LFR > LFR.threshold] = new.Labels$Labels[Diff.LFR > LFR.threshold]
        names(Labels) = Annot[Labels]
    }
    
    Labels.conf = new.Labels$Labels.confidence
    updated.Labels = as.numeric(Labels)
    names(updated.Labels) = names(Annot)[match(Labels, Annot)]
    
    if (adjust.levels == TRUE) {
        Labels = reannotate.labels(ace, updated.Labels)
    } else {
        Labels = updated.Labels
    }
    
    Labels = names(Labels)
    return(Labels)
}

EnhAdj <- function(Adj) {
    
    Adj[is.na(Adj)] = 0
    Adj[Adj < 0] = 0
    
    A = as(Adj, "dgTMatrix")
    diag(A) = 0
    eps = 1e-16
    rs = fastRowSums(A)
    rs[rs == 0] = 1
    P = Matrix::sparseMatrix(i = A@i + 1, j = A@j + 1, x = A@x/rs[A@i + 1], dims = dim(A))
    
    w = sqrt(fastColSums(P) + eps)
    W = P %*% Matrix::Diagonal(x = 1/w, n = length(w))
    P = W %*% Matrix::t(W)
    P = as.matrix(P)
    diag(P) = 0
    
    return(P)
}


construct.tspanner <- function(backbone, stretch.factor = 10) {
    
    backbone[backbone < 0] = 0
    diag(backbone) = 0
    
    backbone.graph = igraph::graph_from_adjacency_matrix(adjmatrix = backbone, mode = "undirected", 
        weighted = TRUE)
    
    # Construct t-spanner
    t = (2 * stretch.factor - 1)
    
    d = 1 - igraph::E(backbone.graph)$weight
    EL = igraph::get.edgelist(backbone.graph, names = FALSE)
    perm = order(d, decreasing = FALSE)
    
    backbone.graph.sparse = igraph::delete_edges(graph = backbone.graph, edges = igraph::E(backbone.graph))
    
    for (i in 1:length(d)) {
        u = EL[perm[i], 1]
        v = EL[perm[i], 2]
        sp = igraph::distances(backbone.graph.sparse, v = u, to = v)[1, 1]
        
        if (sp > t * d[perm[i]]) {
            backbone.graph.sparse = igraph::add_edges(graph = backbone.graph.sparse, 
                edges = EL[perm[i], ], attr = list(weight = 1 - d[perm[i]]))
        }
    }
    
    G = as(igraph::get.adjacency(backbone.graph.sparse, attr = "weight"), "dgCMatrix")
    
    return(G)
}
