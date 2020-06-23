#' Computes a list of archetype-specific interactomes
#'
#' @param ace ACTIONet output object
#' @param G Baseline network as an adjacency matrix or igraph object. 
#' If it is NULL, PCNet is loaded as the baseline from SCINET package.
#' @param min.edge.weight Used in post-processing to remove insignificant edges
#' @param compute.topo.specificity Whether to compute topological specificity scores (will be stores as V(G)$specificity)
#' @param spec.sample_no Number of random samples for computing topological specificity scores
#' @param thread_no Number of parallel threads 
#'  
#' @return A list of igraph objects
#' 
#' @examples
#' network.list = run.SCINET.archetype(ace)
#' G = network.list[[1]]
#' V(G)$name[order(V(G)$specificity, decreasing = T)[1:10]]
run.SCINET.archetype <- function(ace, G = NULL, core = T, min.edge.weight = 2, spec.sample_no = 1000, 
    thread_no = 4, compute.topo.specificity = T) {
    library(SCINET)
    
    print("Preprocessing the baseline interactome")
    if (is.null(G)) {
        if (!exists("PCNet")) {
            data("PCNet")
        }
        Adj = PCNet
    } else if (is.matrix(G) | is.sparseMatrix(G)) {
        Adj = as(G, "sparseMatrix")
        Adj@x = rep(1, length(Adj@x))
    } else if (is.igraph(G)) {
        Adj = as(get.adjacency(G), "sparseMatrix")
    }
    
    gene.scores = rowMaps(ace)[["unified_feature_specificity"]]
    
    
    common.genes = intersect(rownames(gene.scores), rownames(PCNet))
    if (length(common.genes) == 0) {
        print("No common genes found. Check rownames (or vertex names) for the input graph")
        return(ace)
    }
    A = gene.scores[common.genes, ]
    G = Adj[common.genes, common.genes]
    
    
    print("Constructing networks")
    gene.activity.scores = SCINET::compute_gene_activities_full(A = A, thread_no = thread_no)
    cellstate.nets = SCINET::construct_cell_networks(net = G, gene_activities = gene.activity.scores, 
        thread_no = thread_no)
    cellstate.nets.list = as.list(cellstate.nets)
    
    print("Post-processing networks\n")
    cellstate.nets.list.igraph = lapply(cellstate.nets.list, function(G.Adj) {
        G.Adj@x[G.Adj@x < min.edge.weight] = 0
        filter.mask = Matrix::colSums(G.Adj) == 0
        G = igraph::graph_from_adjacency_matrix(G.Adj[!filter.mask, !filter.mask], 
            mode = "undirected", weighted = T)
        V(G)$name = common.genes[!filter.mask]
        if (compute.topo.specificity == TRUE) {
            z.scores = topo.spec(G, spec.sample_no)
            V(G)$specificity = 1/(1 + exp(-z.scores))
        }
        
        return(G)
    })
    
    if (is.null(colnames(gene.scores))) {
        names(cellstate.nets.list.igraph) = 1:ncol(gene.scores)
    } else {
        names(cellstate.nets.list.igraph) = colnames(gene.scores)
    }
    
    return(cellstate.nets.list.igraph)
}

#' Computes a list of cluster-specific interactomes
#'
#' @param ace ACTIONet output object
#' @param specificity.slot.name An entry in the rowMaps(ace), precomputed using compute.cluster.feature.specificity() function
#' @param G Baseline network as an adjacency matrix or igraph object. 
#' If it is NULL, PCNet is loaded as the baseline from SCINET package.
#' @param min.edge.weight Used in post-processing to remove insignificant edges
#' @param compute.topo.specificity Whether to compute topological specificity scores (will be stores as V(G)$specificity)
#' @param spec.sample_no Number of random samples for computing topological specificity scores
#' @param thread_no Number of parallel threads 
#'  
#' @return A list of igraph objects
#' 
#' @examples
#' data('curatedMarkers_human') # pre-packaged in ACTIONet
#' marker.genes = curatedMarkers_human$Blood$PBMC$Monaco2019.12celltypes$marker.genes
#' ace = compute.cluster.feature.specificity(ace, ace$clusters, 'cluster_specificity_scores')
#' network.list = run.SCINET.clusters(ace, 'cluster_specificity_scores')
#' G = network.list[[1]]
#' V(G)$name[order(V(G)$specificity, decreasing = T)[1:10]]
run.SCINET.clusters <- function(ace, specificity.slot.name, G = NULL, min.edge.weight = 2, 
    spec.sample_no = 1000, thread_no = 8, compute.topo.specificity = T) {
    library(SCINET)
    
    print("Preprocessing the baseline interactome")
    if (is.null(G)) {
        if (!exists("PCNet")) {
            data("PCNet")
        }
        Adj = PCNet
    } else if (is.matrix(G) | is.sparseMatrix(G)) {
        Adj = as(G, "sparseMatrix")
        Adj@x = rep(1, length(Adj@x))
    } else if (is.igraph(G)) {
        Adj = as(get.adjacency(G), "sparseMatrix")
    }
    
    
    if (!(specificity.slot.name %in% names(rowMaps(ace)))) {
        message(sprintf("%s does not exist in rowMaps(ace)", specificity.slot.name))
    }
    
    gene.scores = Matrix::t(as.matrix(log1p(t(rowMaps(ace)[[specificity.slot.name]]))))
    
    
    common.genes = intersect(rownames(gene.scores), rownames(PCNet))
    if (length(common.genes) == 0) {
        print("No common genes found. Check rownames (or vertex names) for the input graph")
        return(ace)
    }
    A = gene.scores[common.genes, ]
    G = Adj[common.genes, common.genes]
    
    
    print("Constructing networks")
    gene.activity.scores = SCINET::compute_gene_activities_full(A = A, thread_no = thread_no)
    cellstate.nets = SCINET::construct_cell_networks(net = G, gene_activities = gene.activity.scores, 
        thread_no = thread_no)
    cellstate.nets.list = as.list(cellstate.nets)
    
    print("Post-processing networks\n")
    cellstate.nets.list.igraph = lapply(cellstate.nets.list, function(G.Adj) {
        G.Adj@x[G.Adj@x < min.edge.weight] = 0
        filter.mask = Matrix::colSums(G.Adj) == 0
        G = igraph::graph_from_adjacency_matrix(G.Adj[!filter.mask, !filter.mask], 
            mode = "undirected", weighted = T)
        V(G)$name = common.genes[!filter.mask]
        if (compute.topo.specificity == TRUE) {
            z.scores = topo.spec(G, spec.sample_no)
            V(G)$specificity = 1/(1 + exp(-z.scores))
        }
        
        return(G)
    })
    
    if (is.null(colnames(gene.scores))) {
        names(cellstate.nets.list.igraph) = 1:ncol(gene.scores)
    } else {
        names(cellstate.nets.list.igraph) = colnames(gene.scores)
    }
    return(cellstate.nets.list.igraph)
}


#' Computes a list of specific interactomes for an arbitrary gene score matrix
#'
#' @param gene.scores gene.scores
#' @param specificity.slot.name An entry in the rowMaps(ace), precomputed using compute.cluster.feature.specificity() function
#' @param G Baseline network as an adjacency matrix or igraph object. 
#' If it is NULL, PCNet is loaded as the baseline from SCINET package.
#' @param min.edge.weight Used in post-processing to remove insignificant edges
#' @param compute.topo.specificity Whether to compute topological specificity scores (will be stores as V(G)$specificity)
#' @param spec.sample_no Number of random samples for computing topological specificity scores
#' @param thread_no Number of parallel threads 
#'  
#' @return A list of igraph objects
#' 
#' @examples
#' network.list = run.SCINET.gene.scores(gene.scores)
#' G = network.list[[1]]
#' V(G)$name[order(V(G)$specificity, decreasing = T)[1:10]]
run.SCINET.gene.scores <- function(gene.scores, G = NULL, min.edge.weight = 2, spec.sample_no = 1000, 
    thread_no = 8, compute.topo.specificity = T) {
    require(SCINET)
    
    print("Preprocessing the baseline interactome")
    if (is.null(G)) {
        if (!exists("PCNet")) {
            data("PCNet")
        }
        Adj = PCNet
    } else if (is.matrix(G) | is.sparseMatrix(G)) {
        Adj = as(G, "sparseMatrix")
        Adj@x = rep(1, length(Adj@x))
    } else if (is.igraph(G)) {
        Adj = as(get.adjacency(G), "sparseMatrix")
    }
    
    
    common.genes = intersect(rownames(gene.scores), rownames(Adj))
    if (length(common.genes) == 0) {
        print("No common genes found. Check rownames (or vertex names) for the input graph")
        return()
    }
    A = gene.scores[common.genes, ]
    G = Adj[common.genes, common.genes]
    
    
    
    print("Constructing networks")
    gene.activity.scores = SCINET::RIN_transform(A = A, thread_no = thread_no)
    cellstate.nets = SCINET::construct_cell_networks(net = G, gene_activities = gene.activity.scores, 
        thread_no = thread_no)
    cellstate.nets.list = as.list(cellstate.nets)
    
    print("Post-processing networks\n")
    cellstate.nets.list.igraph = lapply(cellstate.nets.list, function(G.Adj) {
        G.Adj@x[G.Adj@x < min.edge.weight] = 0
        filter.mask = Matrix::colSums(G.Adj) == 0
        G = igraph::graph_from_adjacency_matrix(G.Adj[!filter.mask, !filter.mask], 
            mode = "undirected", weighted = T)
        V(G)$name = common.genes[!filter.mask]
        if (compute.topo.specificity == TRUE) {
            z.scores = topo.spec(G, spec.sample_no)
            V(G)$specificity.z = z.scores
            V(G)$specificity = 1/(1 + exp(-z.scores))
        }
        
        return(G)
    })
    
    if (is.null(colnames(gene.scores))) {
        names(cellstate.nets.list.igraph) = 1:ncol(gene.scores)
    } else {
        names(cellstate.nets.list.igraph) = colnames(gene.scores)
    }
    
    return(cellstate.nets.list.igraph)
}
