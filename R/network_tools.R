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
infer.missing.cell.annotations <- function(
  ace,
  initial_labels,
  iters = 3,
  lambda = 0, 
  sig_threshold = 3) {
	fixed_labels_ = which(!is.na(initial_labels))	
	initial_labels[is.na(initial_labels)] = -1
	
	Labels = run_LPA(ace$ACTIONet, initial_labels, lambda = lambda, iters = iters, sig_threshold = sig_threshold, fixed_labels_ = fixed_labels_)

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
correct.cell.annotations <- function(
  ace,
  initial_labels,
  iters = 3,
  lambda = 0, 
  sig_threshold = 3,
  min.cell.fraction = 0.001
) {
	min_cells = round(ncol(ace)*min.cell.fraction)
	
	label_type = "numeric"
	if(is.character(initial_labels)) {
		label_type = "char"
		initial_labels.factor = factor(initial_labels)
		initial_labels = as.numeric(initial_labels.factor)
	} else if(is.factor(initial_labels)) {
		label_type = "factor"
		initial_labels.factor = initial_labels		
		initial_labels = as.numeric(initial_labels.factor)
	}
	
	cc = table(initial_labels)
	initial_labels[initial_labels %in% as.numeric(names(cc)[cc < min_cells])] = -1
	initial_labels[is.na(initial_labels)] = -1
	
	Labels = run_LPA(ace$ACTIONet, initial_labels, lambda = lambda, iters = iters, sig_threshold = sig_threshold)

	if(label_type == "char" | label_type == "factor") {
		Labels = levels(initial_labels.factor)[initial_labels.factor]
	}

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
    P = Matrix::sparseMatrix(
      i = A@i + 1,
      j = A@j + 1,
      x = A@x/rs[A@i + 1],
      dims = dim(A)
    )

    w = sqrt(fastColSums(P) + eps)
    W = P %*% Matrix::Diagonal(x = 1/w, n = length(w))
    P = W %*% Matrix::t(W)
    P = as.matrix(P)
    diag(P) = 0

    return(P)
}


construct.tspanner <- function(
  backbone,
  stretch.factor = 10
) {

    backbone[backbone < 0] = 0
    diag(backbone) = 0

    backbone.graph = igraph::graph_from_adjacency_matrix(
      adjmatrix = backbone,
      mode = "undirected",
      weighted = TRUE
    )

    # Construct t-spanner
    t = (2 * stretch.factor - 1)

    d = 1 - igraph::E(backbone.graph)$weight
    EL = igraph::get.edgelist(backbone.graph, names = FALSE)
    perm = order(d, decreasing = FALSE)

    backbone.graph.sparse = igraph::delete_edges(
      graph = backbone.graph,
      edges = igraph::E(backbone.graph)
    )

    for (i in 1:length(d)) {
        u = EL[perm[i], 1]
        v = EL[perm[i], 2]
        sp = igraph::distances(backbone.graph.sparse, v = u, to = v)[1, 1]

        if (sp > t * d[perm[i]]) {
            backbone.graph.sparse = igraph::add_edges(
              graph = backbone.graph.sparse,
              edges = EL[perm[i], ],
              attr = list(weight = 1 - d[perm[i]])
            )
        }
    }

    G = as(igraph::get.adjacency(backbone.graph.sparse, attr = "weight"), "dgCMatrix")

    return(G)
}
