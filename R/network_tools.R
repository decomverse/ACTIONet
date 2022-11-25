#' @export
networkDiffusion <- function(
  obj,
  scores, ## `scores` can be a matrix
  algorithm = c("pagerank", "pagerank_sym"),
  alpha = 0.9,
  thread_no = 0,
  max_it = 5,
  res_threshold = 1e-8,
  net_slot = "ACTIONet"
) {

  algorithm <- tolower(algorithm)
  algorithm <- match.arg(algorithm, several.ok = FALSE)

  scores = Matrix::as.matrix(scores)
  if ( NROW(scores) != NCOL(obj) ){
    err = sprintf("`length(scores)` must equal `NCOL(obj)`.\n")
    stop(err)
  }

  G <- .ace_or_net(
    obj = obj,
    net_slot = net_slot,
    matrix_type = "sparse",
    force_type = TRUE,
    obj_name = "obj"
  )

  if(alpha >= 1){
    stop("`alpha` => 1")
  } else if (alpha < 0) {
    stop("`alpha` < 0")
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
  obj,
  label_attr = NULL,
  algorithm = c("coreness", "pagerank", "local_coreness", "local_pagerank"),
  alpha = 0.9,
  net_slot = "ACTIONet",
  diffusion_it = 5,
  res_threshold = 1e-8,
  thread_no = 0
) {

  algorithm <- tolower(algorithm)
  algorithm <- match.arg(algorithm, several.ok = FALSE)

  G <- .ace_or_net(
    obj = obj,
    net_slot = net_slot,
    matrix_type = "sparse",
    force_type = TRUE,
    obj_name = "obj"
  )

  if( algorithm %in% c("pagerank", "local_pagerank") ) {
    if( is.null(label_attr) ){
      err = sprintf("'label_attr' cannot be 'NULL' if 'algorithm=%s'.\n", algorithm)
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

  if (!is.null(label_attr)) {
    label_attr <- .validate_attr(obj, attr = label_attr, obj_name = "obj", attr_name = "label_attr")
    assignments = as.numeric(factor(label_attr))
  }

  if (algorithm == "coreness") {

    centrality <- compute_core_number(G)

  } else if (algorithm == "pagerank") {

    centrality <- networkDiffusion(
      obj = G,
      scores = rep(1 / NCOL(G), NCOL(G)),
      algorithm = "pagerank",
      alpha = alpha,
      thread_no = thread_no,
      max_it = diffusion_it,
      res_threshold = res_threshold,
      net_slot = NULL
    )

  } else if (algorithm == "local_coreness") {

    centrality <- compute_archetype_core_centrality(G, assignments)

  } else if (algorithm == "local_pagerank") {

    design.mat <- model.matrix(~ 0 + as.factor(assignments))
    design.mat <- scale(design.mat, center = FALSE, scale = colSums(design.mat))
    scores <- networkDiffusion(
      obj = G,
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

  centrality = c(centrality)
  return(centrality)

}


#' @export
networkPropagation <- function(
  obj,
  label_attr,
  fixed_samples = NULL,
  algorithm = c("lpa"),
  lambda = 0,
  iters = 3,
  sig_th = 3,
  net_slot = "ACTIONet",
  thread_no = 0
) {

  algorithm <- tolower(algorithm)
  algorithm <- match.arg(algorithm, several.ok = FALSE)

  G <- .ace_or_net(
    obj = obj,
    net_slot = net_slot,
    matrix_type = "sparse",
    force_type = TRUE,
    obj_name = "obj"
  )

  lf <- factor(.validate_attr(obj, attr = label_attr, obj_name = "obj", attr_name = "label_attr"))
  labels = as.numeric(lf)
  keys = levels(lf)
  labels[is.na(labels)] <- -1

  if (algorithm == "lpa") {
    new_labels <- run_LPA(G = G, labels = labels, lambda = lambda, iters = iters, sig_threshold = sig_th, fixed_labels_ = fixed_samples, thread_no = thread_no)
    new_labels[new_labels == -1] <- NA
    new_labels <- keys[new_labels]
  } else {
    warning("No propagation was performed");
    new_labels[new_labels == -1] <- NA
    new_labels <- keys[labels]
  }

  return(new_labels)
}


#' @export
correct.cell.labels  <- function(
  ace,
  label_attr,
  algorithm = "LPA",
  iters = 3,
  lambda = 0,
  sig_th = 3,
  net_slot = "ACTIONet",
  thread_no = 0
) {

  initial_labels = ACTIONetExperiment::get.data.or.split(ace, attr = label_attr, to_return = "data")
  labels <- networkPropagation(
    obj = ace,
    label_attr = initial_labels,
    fixed_samples = NULL,
    algorithm = algorithm,
    lambda = lambda,
    iters = iters,
    sig_th = sig_th,
    net_slot = net_slot,
    thread_no = thread_no
  )

  return(labels)
}

#' @export
infer.missing.cell.labels  <- function(
  ace,
  label_attr,
  algorithm = "LPA",
  iters = 3,
  lambda = 0,
  sig_th = 3,
  net_slot = "ACTIONet",
   thread_no = 0
) {

  initial_labels = ACTIONetExperiment::get.data.or.split(ace, attr = label_attr, to_return = "data")
  fixed_samples <- which(!is.na(initial_labels))

  labels <- networkPropagation(
    obj = ace,
    label_attr = initial_labels,
    fixed_samples = fixed_samples,
    algorithm = algorithm,
    lambda = lambda,
    iters = iters,
    sig_th = sig_th,
    net_slot = net_slot,
    thread_no = thread_no
  )

  return(labels)
}


networkAutocorrelation <- function(
  obj,
  scores = NULL,
  algorithm = c("geary", "moran", "categorical"),
  score_normalization_method = 1L,
  perm_no = 0,
  thread_no = 0L,
  net_slot = "ACTIONet"
) {

  algorithm <- tolower(algorithm)
  algorithm <- match.arg(algorithm, several.ok = FALSE)

  G <- .ace_or_net(
    obj = obj,
    net_slot = net_slot,
    obj_name = "obj"
  )

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

  G <- as(igraph::get.adjacency(backbone.graph.sparse, attr = "weight"), "dMatrix")

  return(G)
}
