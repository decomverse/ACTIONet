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
#' ace <- infer.missing.cell.annotations(ace, sce$assigned_archetypes, "updated_archetype_annotations")
#' @export
infer.missing.cell.annotations <- function(ace,
                                           initial_labels,
                                           iters = 3,
                                           lambda = 0,
                                           sig_threshold = 3, net_slot = "ACTIONet") {
  # label_type <- "numeric"
  # if (is.character(initial_labels)) {
  #   label_type <- "char"
  #   initial_labels.factor <- factor(initial_labels)
  #   initial_labels <- as.numeric(initial_labels.factor)
  # } else if (is.factor(initial_labels)) {
  #   label_type <- "factor"
  #   initial_labels.factor <- initial_labels
  #   initial_labels <- as.numeric(initial_labels.factor)
  # }

  # fixed_labels_ <- which(!is.na(initial_labels))
  # initial_labels[is.na(initial_labels)] <- -1

  # Labels <- run_LPA(ace$ACTIONet, initial_labels, lambda = lambda, iters = iters, sig_threshold = sig_threshold, fixed_labels_ = fixed_labels_)

  # Labels[Labels == -1] <- NA
  # if (label_type == "char" | label_type == "factor") {
  #   Labels <- levels(initial_labels.factor)[Labels]
  # }

  fixed_samples <- which(!is.na(initial_labels))
  Labels <- networkPropagation(G = colNets(ace)[[net_slot]], initial_labels = initial_labels, iters = iters, lambda = lambda, sig_threshold = sig_threshold, fixed_samples = fixed_samples)

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
#' @param iters How many iterative rounds of correction/inference should be performed (default=3)
#'
#' @return ace with updated annotations added to ace$annotations
#'
#' @examples
#' ace <- add.cell.annotations(ace, cell.labels, "input_annotations")
#' ace <- correct.cell.annotations(ace, "input_annotations", "updated_annotations")
#' @export
correct.cell.annotations <- function(ace,
                                     initial_labels,
                                     algorithm = "lpa",
                                     iters = 3,
                                     lambda = 0,
                                     sig_threshold = 3,
                                     net_slot = "ACTIONet") {
  algorithm <- tolower(algorithm)

  # label_type <- "numeric"
  # if (is.character(initial_labels)) {
  #   label_type <- "char"
  #   initial_labels.factor <- factor(initial_labels)
  #   initial_labels <- as.numeric(initial_labels.factor)
  # } else if (is.factor(initial_labels)) {
  #   label_type <- "factor"
  #   initial_labels.factor <- initial_labels
  #   initial_labels <- as.numeric(initial_labels.factor)
  # }

  # cc <- table(initial_labels)
  # initial_labels[initial_labels %in% as.numeric(names(cc)[cc < min_cells])] <- -1
  # initial_labels[is.na(initial_labels)] <- -1

  # Labels <- run_LPA(ace$ACTIONet, initial_labels, lambda = lambda, iters = iters, sig_threshold = sig_threshold)
  # Labels[Labels == -1] <- NA

  # if (label_type == "char" | label_type == "factor") {
  #   Labels <- levels(initial_labels.factor)[Labels]
  # }
  Labels <- networkPropagation(G = colNets(ace)[[net_slot]], initial_labels = initial_labels, iters = iters, lambda = lambda, sig_threshold = sig_threshold)


  return(Labels)
}

EnhAdj <- function(Adj) {
  Adj[is.na(Adj)] <- 0
  Adj[Adj < 0] <- 0

  A <- as(Adj, "dgTMatrix")
  diag(A) <- 0
  eps <- 1e-16
  rs <- fastRowSums(A)
  rs[rs == 0] <- 1
  P <- Matrix::sparseMatrix(
    i = A@i + 1,
    j = A@j + 1,
    x = A@x / rs[A@i + 1],
    dims = dim(A)
  )

  w <- sqrt(Matrix::colSums(P) + eps)
  W <- P %*% Matrix::Diagonal(x = 1 / w, n = length(w))
  P <- W %*% Matrix::t(W)
  P <- as.matrix(P)
  diag(P) <- 0

  return(P)
}


construct.tspanner <- function(backbone,
                               stretch.factor = 10) {
  backbone[backbone < 0] <- 0
  diag(backbone) <- 0

  backbone.graph <- igraph::graph_from_adjacency_matrix(
    adjmatrix = backbone,
    mode = "undirected",
    weighted = TRUE
  )

  # Construct t-spanner
  t <- (2 * stretch.factor - 1)

  d <- 1 - igraph::E(backbone.graph)$weight
  EL <- igraph::get.edgelist(backbone.graph, names = FALSE)
  perm <- order(d, decreasing = FALSE)

  backbone.graph.sparse <- igraph::delete_edges(
    graph = backbone.graph,
    edges = igraph::E(backbone.graph)
  )

  for (i in 1:length(d)) {
    u <- EL[perm[i], 1]
    v <- EL[perm[i], 2]
    sp <- igraph::distances(backbone.graph.sparse, v = u, to = v)[1, 1]

    if (sp > t * d[perm[i]]) {
      backbone.graph.sparse <- igraph::add_edges(
        graph = backbone.graph.sparse,
        edges = EL[perm[i], ],
        attr = list(weight = 1 - d[perm[i]])
      )
    }
  }

  G <- as(igraph::get.adjacency(backbone.graph.sparse, attr = "weight"), "dgCMatrix")

  return(G)
}

networkDiffusion <- function(
  G,
  scores,
  algorithm = "pagerank",
  alpha = 0.9,
  thread_no = 0,
  max_it = 5,
  res_threshold = 1e-8,
  net_slot = "ACTIONet"
) {

  algorithm <- match.arg(algorithm)

  if(alpha >= 1){
    alpha = 1
  } else if (alpha <= 0) {
    alpha = 0
  }

  if (is(G, "ACTIONetExperiment")) {
    G <- colNets(G)[[net_slot]]
  }

  if (algorithm == "pagerank") {
    x <- compute_network_diffusion_approx(G = G, X0 = as.matrix(scores), thread_no = thread_no, alpha = alpha, max_it = max_it, res_threshold = res_threshold, norm_type = 0)
  } else if (algorithm == "pagerank_sym") {
    x <- compute_network_diffusion_approx(G = G, X0 = as.matrix(scores), thread_no = thread_no, alpha = alpha, max_it = max_it, res_threshold = res_threshold, norm_type = 2)
  }

  return(x)
}

networkPropagation <- function(
  G,
  initial_labels,
  algorithm = "LPA",
  lambda = 0,
  iters = 3,
  sig_threshold = 3,
  net_slot = "ACTIONet",
  fixed_samples = NULL
) {

  algorithm <- match.arg(algorithm)

  if (algorithm == "lpa") {
    label_type <- "numeric"
    if (is.character(initial_labels)) {
      label_type <- "char"
      initial_labels.factor <- factor(initial_labels)
      initial_labels <- as.numeric(initial_labels.factor)
    } else if (is.factor(initial_labels)) {
      label_type <- "factor"
      initial_labels.factor <- initial_labels
      initial_labels <- as.numeric(initial_labels.factor)
    }
    initial_labels[is.na(initial_labels)] <- -1

    if (is(G, "ACTIONetExperiment")) {
      G <- colNets(G)[[net_slot]]
    }
    updated_labels <- run_LPA(G, initial_labels, lambda = lambda, iters = iters, sig_threshold = sig_threshold)
    updated_labels[updated_labels == -1] <- NA

    if (label_type == "char" | label_type == "factor") {
      updated_labels <- levels(initial_labels.factor)[updated_labels]
    }
  }

  return(updated_labels)
}


networkCentrality <- function(
  ace = NULL,
  G = NULL,
  label_attr = NULL,
  algorithm = c("coreness", "pagerank", "localized_coreness", "localized_pagerank"),
  alpha = 0.9,
  net_slot = "ACTIONet",
  centrality_attr = "node_centrality",
  return_raw = FALSE
) {

  algorithm <- match.arg(algorithm)

  if( is.null(ace) && is.null(G) ){
    err = sprintf("Either 'ace' or 'G' must be given.\n")
    stop(err)
  }

  if (is.null(G)) {
    if (!(net_slot %in% names(colNets(ace)))) {
      err <- sprintf("Attribute '%s' is not in 'colNets'.\n", net_slot)
      stop(err)
    }
    G <- colNets(ace)[[net_slot]]
  } else {
    G =  as(G, "dgCMatrix")
  }

  if(alpha >= 1){
    alpha = 0.99
    warning("alpha=0.99")
  } else if (alpha < 0) {
    alpha = 0
    warning("alpha=0")
  }

  if( algorithm %in% c("pagerank", "localized_pagerank") ) {
    if(is.null(label_attr)){
      err = sprintf("'label_attr' cannot be 'NULL' if 'algorithm=%s'.\n", algorithm)
      stop(err)
    }
    # if (alpha == 0) {
    #   err = sprintf("'alpha' must be non-zero positive when 'algorithm=%s'.\n", algorithm)
    #   stop(err)
    # }
  }


  if (!is.null(label_attr)) {
    if(!is.null(ace)){
        assignments = ACTIONetExperiment::get.data.or.split(ace, attr = label_attr, to_return = "levels")$index
    } else {
        assignments = as.numeric(factor(label_attr))
        if ( length(assignments) != NROW(G) ){
          err = sprintf("'label_attr' does not match 'NROW(G)'.\n")
          stop(err)
        }
    }
  }

  if (algorithm == "coreness") {
    centrality <- compute_core_number(G)
  } else if (algorithm == "pagerank") {
    centrality <- networkDiffusion(G = G, scores = as.matrix(rep(1 / NROW(G), NROW(G))), alpha = alpha)
  } else if (algorithm == "localized_coreness") {
    centrality <- compute_archetype_core_centrality(G, assignments)
  } else if (algorithm == "localized_pagerank") {
    design.mat <- model.matrix(~ 0 + as.factor(assignments))
    design.mat <- scale(design.mat, center = FALSE, scale = colSums(design.mat))
    scores <- networkDiffusion(G = G, scores = as.matrix(design.mat), alpha = alpha)
    scores <- apply(scores, 2, function(x) x / max(x))
    centrality <- apply(scores, 1, max)
  }

  centrality = c(centrality)

  if (return_raw == TRUE || is.null(ace)){
      return(centrality)
  } else {
    SummarizedExperiment::colData(ace)[[centrality_attr]] = centrality
    return(ace)
  }

}


networkAutocorrelation <- function(G, scores = NULL, algorithm = "Geary", score_normalization_method = 1L, perm_no = 0, thread_no = 0L, net_slot = "ACTIONet") {
  algorithm <- tolower(algorithm)

  if (is(G, "ACTIONetExperiment")) {
    G <- colNets(G)[[net_slot]]
  }

  if (algorithm == "geary") {
    if (is.sparseMatrix(G)) {
      out <- ACTIONet::autocorrelation_Geary(G = G, scores = scores, normalization_method = score_normalization_method, perm_no = perm_no, thread_no = thread_no)
    } else {
      out <- ACTIONet::autocorrelation_Geary_full(G = G, scores = scores, normalization_method = score_normalization_method, perm_no = perm_no, thread_no = thread_no)
    }
  } else if (algorithm == "moran") {
    if (is.sparseMatrix(G)) {
      out <- ACTIONet::autocorrelation_Moran(G = G, scores = scores, normalization_method = score_normalization_method, perm_no = perm_no, thread_no = thread_no)
    } else {
      out <- ACTIONet::autocorrelation_Moran_full(G = G, scores = scores, normalization_method = score_normalization_method, perm_no = perm_no, thread_no = thread_no)
    }
  } else if (algorithm == "categorical") {
    if (is.character(scores)) {
      label_type <- "char"
      annotations.factor <- factor(scores)
      annotations <- as.numeric(annotations.factor)
    } else if (is.factor(scores)) {
      label_type <- "factor"
      annotations.factor <- scores
      annotations <- as.numeric(annotations.factor)
    } else {
      annotations <- scores
    }
    out <- assess.categorical.autocorrelation(A = G, labels = annotations, perm.no = perm_no)
  }

  return(out)
}



literallyAnyOtherName <- function(ace, z_threshold = 1, global = FALSE, net_slot = "ACTIONet") {
  G <- colNets(ace)[[net_slot]]
  if (global == TRUE) {
    cn <- compute_core_number(G)
  } else {
    cn <- ace$node_centrality
  }
  x <- networkDiffusion(G, scores = as.matrix(cn))
  z <- -(x - median(x)) / mad(x)

  mask <- z_threshold < z

  out <- list(score = exp(z), is_filtered = mask, z.score = z)
  return(out)
}
