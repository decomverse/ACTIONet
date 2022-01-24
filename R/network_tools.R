#' @export
networkDiffusion <- function(
  G = NULL,
  scores = NULL,
  algorithm = c("pagerank", "pagerank_sym"),
  alpha = 0.9,
  thread_no = 0,
  max_it = 5,
  res_threshold = 1e-8,
  net_slot = "ACTIONet"
) {

  algorithm <- match.arg(algorithm)

  if (is(G, "ACTIONetExperiment")) {
    if (!(net_slot %in% names(colNets(G)))) {
      err <- sprintf("Attribute '%s' is not in 'colNets'.\n", net_slot)
      stop(err)
    }
    G <- colNets(G)[[net_slot]]
  }
  if( is.null(G) ){
    err = sprintf("'G' must be given.\n")
    stop(err)
  }
  G =  as(G, "dgCMatrix")

  if( is.null(scores) ){
    err = sprintf("'scores' cannot be 'NULL'.\n")
    stop(err)
  }
  scores = Matrix::as.matrix(scores)
  if ( nrow(scores) != NROW(G) ){
    err = sprintf("'length(scores)' must equal 'NROW(G)'.\n")
    stop(err)
  }

  if(alpha >= 1){
    stop("alpha=>1")
  } else if (alpha < 0) {
    stop("alpha=<0")
  }


  if (algorithm == "pagerank") {
    x <- compute_network_diffusion_approx(
      G = G,
      X0 = scores,
      thread_no = thread_no,
      alpha = alpha,
      max_it = max_it,
      res_threshold = res_threshold,
      norm_type = 0
    )
  } else if (algorithm == "pagerank_sym") {
    x <- compute_network_diffusion_approx(
      G = G,
      X0 = scores,
      thread_no = thread_no,
      alpha = alpha,
      max_it = max_it,
      res_threshold = res_threshold,
      norm_type = 2
    )
  }

  return(x)
}


#' @export
networkCentrality <- function(
  ace = NULL,
  G = NULL,
  labels = NULL,
  algorithm = c("coreness", "pagerank", "localized_coreness", "localized_pagerank"),
  alpha = 0.9,
  net_slot = "ACTIONet",
  diffusion_it = 5,
  res_threshold = 1e-8,
  thread_no = 0
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

  if( algorithm %in% c("pagerank", "localized_pagerank") ) {
    if( is.null(labels) ){
      err = sprintf("'labels' cannot be 'NULL' if 'algorithm=%s'.\n", algorithm)
      stop(err)
    }

    if(alpha >= 1){
      alpha = 0.99
      warning("alpha=0.99")
    } else if (alpha < 0) {
      alpha = 0
      warning("alpha=0")
    }
  }

  if (!is.null(labels)) {
    if(!is.null(ace)){
        assignments = ACTIONetExperiment::get.data.or.split(ace, attr = labels, to_return = "levels")$index
    } else {
        assignments = as.numeric(factor(labels))
        if ( length(assignments) != NROW(G) ){
          err = sprintf("'length(labels)' must equal 'NROW(G)'.\n")
          stop(err)
        }
    }
  }

  if (algorithm == "coreness") {

    centrality <- compute_core_number(G)

  } else if (algorithm == "pagerank") {

    centrality <- networkDiffusion(
      G = G,
      scores = rep(1 / NROW(G), NROW(G)),
      algorithm = "pagerank",
      alpha = alpha,
      thread_no = thread_no,
      max_it = diffusion_it,
      res_threshold = res_threshold,
      net_slot = NULL
    )

  } else if (algorithm == "personalized_coreness") {

    centrality <- compute_archetype_core_centrality(G, assignments)

  } else if (algorithm == "personalized_pagerank") {

    design.mat <- model.matrix(~ 0 + as.factor(assignments))
    design.mat <- scale(design.mat, center = FALSE, scale = colSums(design.mat))
    scores <- networkDiffusion(
      G = G,
      scores = design.mat,
      algorithm = "pagerank",
      alpha = alpha,
      thread_no = thread_no,
      max_it = diffusion_it,
      res_threshold = res_threshold,
      net_slot = NULL
    )
    scores <- apply(scores, 2, function(x) x / max(x))
    centrality <- apply(scores, 1, max)

  }

  return(centrality)

}


#' @export
networkPropagation <- function(
  G = NULL,
  labels = NULL,
  fixed_samples = NULL,
  algorithm = c("LPA"),
  lambda = 0,
  iters = 3,
  sig_th = 3,
  net_slot = "ACTIONet"
) {
 algorithm <- match.arg(algorithm)

  if (is(G, "ACTIONetExperiment")) {
    if (!(net_slot %in% names(colNets(G)))) {
      err <- sprintf("Attribute '%s' is not in 'colNets'.\n", net_slot)
      stop(err)
    }
    G <- colNets(G)[[net_slot]]
  }
  if( is.null(G) ){
    err = sprintf("'G' must be given.\n")
    stop(err)
  }
  G =  as(G, "dgCMatrix")

  if( is.null(labels) ){
    err = sprintf("'scores' cannot be 'NULL'.\n")
    stop(err)
  }
  
  if ( length(labels) != NROW(G) ){
    err = sprintf("'length(labels)' must equal 'NROW(G)'.\n")
    stop(err)
  }
  lf = factor(labels)
  labels = as.numeric(lf)
  keys = levels(lf)

  labels[is.na(labels)] <- -1

  if (algorithm == "LPA") {
    new_labels <- run_LPA(G = G, labels = labels, lambda = lambda, iters = iters, sig_threshold = sig_th, fixed_labels_ = fixed_samples)
    new_labels[new_labels == -1] <- NA
    new_labels <- keys[new_labels]
  } else {
    new_labels <- keys[labels]
  }

  return(new_labels)
}


#' @export
infer.missing.cell.labels <- function(
  ace,
  labels,
  iters = 3,
  lambda = 0,
  sig_th = 3,
  net_slot = "ACTIONet"
) {


  initial_labels = ACTIONetExperiment::get.data.or.split(ace, attr = labels, to_return = "data")
  fixed_samples <- which(!is.na(initial_labels))

  labels <- networkPropagation <- function(
    G = NULL,
    labels = initial_labels,
    fixed_samples = fixed_samples,
    algorithm = "LPA",
    lambda = lambda,
    iters = iters,
    sig_th = sig_th,
    net_slot = net_slot
  )

  return(labels)
}


#' @export
correct.cell.labels <- function(ace,
                                     labels,
                                     algorithm = "LPA",
                                     iters = 3,
                                     lambda = 0,
                                     sig_th = 3,
                                     net_slot = "ACTIONet") {

  # Labels <- networkPropagation(G = colNets(ace)[[net_slot]], initial_labels = initial_labels, iters = iters, lambda = lambda, sig_th = sig_th)

  initial_labels = ACTIONetExperiment::get.data.or.split(ace, attr = labels, to_return = "data")
  labels <- networkPropagation <- function(
    G = NULL,
    labels = initial_labels,
    fixed_samples = fixed_samples,
    algorithm = algorithm,
    lambda = lambda,
    iters = iters,
    sig_th = sig_th,
    net_slot = net_slot
  )

  return(labels)
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


networkCentrality <- function(G, annotations = NULL, algorithm = "personalized_coreness", alpha_val = 0.85, net_slot = "ACTIONet") {
  algorithm <- tolower(algorithm)

  if (is(G, "ACTIONetExperiment")) {
    G <- colNets(G)[[net_slot]]
  }

  if (!is.null(annotations)) {
    label_type <- "numeric"
    if (is.character(annotations)) {
      label_type <- "char"
      annotations.factor <- factor(annotations)
      annotations <- as.numeric(annotations.factor)
    } else if (is.factor(annotations)) {
      label_type <- "factor"
      annotations.factor <- annotations
      annotations <- as.numeric(annotations.factor)
    }
  }

  if (algorithm == "personalized_coreness") {
    centrality <- compute_archetype_core_centrality(G, annotations)
  } else if (algorithm == "personalized_pagerank") {
    design <- model.matrix(~ .0 + as.factor(annotations))
    design <- scale(design, center = F, scale = fastColSums(design))
    scores <- networkDiffusion(G = G, scores = as.matrix(design), alpha = alpha_val)
    scores <- apply(scores, 2, function(x) x / max(x))
    centrality <- apply(scores, 1, max)
  } else if (algorithm == "pagerank") {
    centrality <- networkDiffusion(G = G, scores = as.matrix(rep(1 / nrow(G), nrow(G))), alpha = alpha_val)
  } else if (algorithm == "coreness") {
    centrality <- compute_core_number(G)
  }

  return(centrality)
}


networkAutocorrelation <- function(G, scores = NULL, algorithm = "Geary", score_normalization_method = 1L, perm_no = 0, thread_no = 0L, network_slot = "ACTIONet") {
  algorithm <- tolower(algorithm)

  if (is(G, "ACTIONetExperiment")) {
    G <- colNets(G)[[network_slot]]
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



literallyAnyOtherName <- function(ace, z_threshold = 1, global = FALSE, network_slot = "ACTIONet") {
  G <- colNets(ace)[[network_slot]]
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
