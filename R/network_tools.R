#' A wrapper function For Leiden algorithm
#'
#' @param ACTIONet.out Input results to be clustered
#' @param annotation.name Arbitrary name to be given to the final clustering for future reference.
#' @param resolution_parameter Resolution of the clustering.
#' The higher the resolution, the more clusters we will get (default=0.5).
#' @param arch.init Whether to use archetype-assignments to initialize clustering (default=TRUE)
#' 
#' @return ACTIONet.out with added annotation
#' 
#' @examples
#' ACTIONet.out = cluster.ACTIONet(ACTIONet.out, "Leiden_clusters")
#' plot.ACTIONet(ACTIONet.out, "Leiden_clusters")
#' clusters = ACTIONet.out$annotations$Leiden_clusters$Labels
cluster.ACTIONet <- function(ACTIONet.out, annotation.name = NULL, resolution_parameter = 0.5, arch.init = TRUE) {
	if( !("NetLibR" %in% rownames(installed.packages())) ) {
		message("You need to install NetLibR (https://github.com/shmohammadi86/NetLibR) first to use graph-based clustering.")
		return
	} else {
		library(NetLibR)
	}
	
    if (arch.init == TRUE) {
		if(! ('unification.out' %in% names(ACTIONet.out)) ) {
			print("unification.out is not in ACTIONet.out. Please run unify.cell.states() first.")
			return(ACTIONet.out)
		}
        print("Perform archetype-based initialization")
		initial.clusters = ACTIONet.out$unification.out$assignments.core
        clusters = as.numeric(unsigned_cluster(ACTIONet.out$build.out$ACTIONet, resolution_parameter, 0, initial.clusters))
    } else {
        print("Perform default initialization")
        clusters = as.numeric(unsigned_cluster(ACTIONet.out$build.out$ACTIONet, resolution_parameter, 0))
    }    
    names(clusters) = paste("Cluster", as.character(clusters), sep = " ")
	
	if(! ('annotations' %in% names(ACTIONet.out)) ) {
		ACTIONet.out$annotations = list()
	}

	time.stamp = as.character(Sys.time())
	if(is.null(annotation.name)) {
		annotation.name = sprintf('Clustering_%s', time.stamp)
	}

	res = list(Labels = clusters, Labels.confidence = NULL, DE.profile = NULL, highlight = NULL, cells = ACTIONet.out$log$cells, time.stamp = time.stamp, type = "cluster_ACTIONet")	

	ACTIONet.out$annotations[[annotation.name]] = res	
		
    return(ACTIONet.out)
}


#' Gene expression imputation using network diffusion.
#'
#' @param ACTIONet.out Input results to be clustered
#' (alternatively it can be the ACTIONet igraph object)
#' @param sce Input single-cell profile witht the normalized counts
#' @param genes The list of genes to perform imputation for.
#' @param alpha_val Depth of diffusion between (0, 1).
#' The larger it is, the deeper the diffusion, which results in less nonzeros (default = 0.85).
#' @param prune Whether to post-process diffusion results to remove false-positives (default = FALSE)
#' @param thread_no Number of parallel threads
#' @param diffusion_iters Number of diffusion iterations (default = 5)
#' @param sce.data.attr Slot in the sce object with normalized counts.
#' 
#' @return Imputed gene expression matrix. Column names are set with imputed genes names and rows are cells.
#' 
#' @examples
#' imputed.genes = impute.genes.using.ACTIONet(ACTIONet.out, sce, c("MOBP", "C3", "GAD1", "VCAN"))
#' plot.ACTIONet.gradient(ACTIONet.out, imputed.genes[, 1])
impute.genes.using.ACTIONet <- function(ACTIONet.out, sce, genes, alpha_val = 0.85, thread_no = 8, prune = FALSE, diffusion_iters = 5, sce.data.attr = "logcounts") {
    genes = unique(genes)
    
    if (igraph::is.igraph(ACTIONet.out)) 
        ACTIONet = ACTIONet.out else ACTIONet = ACTIONet.out$ACTIONet
    
    if (length(igraph::V(ACTIONet)) != ncol(sce)) {
        R.utils::printf("Number of cells in the input sce (%d) doesn't match the number of vertices in the ACTIONet (%d)\n", dim(sce)[2], 
            length(V(ACTIONet)))
    }
    G = as(igraph::get.adjacency(ACTIONet, attr = "weight"), "dgTMatrix")
    
    
    matched.genes = intersect(genes, rownames(sce))
    matched.idx = match(matched.genes, rownames(sce))
    
    # Smooth/impute gene expressions
    if (!(sce.data.attr %in% names(SummarizedExperiment::assays(sce)))) {
        R.utils::printf("%s is not in assays of sce\n", sce.data.attr)
    }
    
    if (length(matched.idx) > 1) {
        raw.gene.expression = Matrix::t(as(SummarizedExperiment::assays(sce)[[sce.data.attr]][matched.idx, ], "dgTMatrix"))
        U = raw.gene.expression
        U[U < 0] = 0
        cs = Matrix::colSums(U)
        U = Matrix::sparseMatrix(i = U@i + 1, j = U@j + 1, x = U@x/cs[U@j + 1], dims = dim(U))
        U = U[, cs > 0]
        gg = matched.genes[cs > 0]
    } else {
        raw.gene.expression = SummarizedExperiment::assays(sce)[[sce.data.attr]][matched.idx, ]
        U = raw.gene.expression/sum(raw.gene.expression)
        gg = matched.genes
    }
    
	# Perform network-diffusion
    G = ACTIONet.out$build.out$ACTIONet
	#imputed.gene.expression = batchPR(G, as.matrix(U), alpha_val, thread_no)
	imputed.gene.expression = PageRank_iter(G, as(U, 'sparseMatrix'), alpha_val, max_it)
	    
    imputed.gene.expression[is.na(imputed.gene.expression)] = 0
    
    # Prune values
    if (prune == TRUE) {
		imputed.gene.expression = prune.scores.using.ACTIONet(ACTIONet.out, imputed.gene.expression);
    }
    
    # Rescale the baseline expression of each gene
	imputed.gene.expression = sapply(1:dim(imputed.gene.expression)[2], function(col) {
		x = raw.gene.expression[, col]
		y = imputed.gene.expression[, col]
		
		x.Q = quantile(x, 1)
		y.Q = quantile(y, 1)
		
		if (y.Q == 0) {
			return(array(0, length(x)))
		}
		
		y = y * x.Q/y.Q
		
		y[y > max(x)] = max(x)
		
		return(y)
	})
    
    colnames(imputed.gene.expression) = gg
    
    return(imputed.gene.expression)
}


#' Uses sweepcup algorithm on the ACTIONet to prune/filter a given score per node.
#'
#' @param ACTIONet.out Input results to be clustered
#' (alternatively it can be the ACTIONet igraph object)
#' @param scores A vector with one value per ACTIONet node
#' 
#' @return Pruned scores
#' 
#' @examples
#' pruned.scores = prune.scores.using.ACTIONet(ACTIONet.out, scores)
prune.scores.using.ACTIONet <- function(ACTIONet.out, scores) {
    if (igraph::is.igraph(ACTIONet.out)) 
        ACTIONet = ACTIONet.out else ACTIONet = ACTIONet.out$ACTIONet
    
    if (length(igraph::V(ACTIONet)) != ncol(sce)) {
        R.utils::printf("Number of cells in the input sce (%d) doesn't match the number of vertices in the ACTIONet (%d)\n", dim(sce)[2], 
            length(V(ACTIONet)))
    }
    G = as(igraph::get.adjacency(ACTIONet, attr = "weight"), "dgTMatrix")
    	
    # To be re-implemented
	# pruned.scores = prune_scores(G, scores)
	pruned.scores = scores
}

#' Uses a variant of the label propagation algorithm to infer missing labels
#'
#' @param ACTIONet.out Input results to be clustered
#' (alternatively it can be the ACTIONet igraph object)
#' @param annotation.in Annotations to correct with missing values (NA) in it.
#' It can be either a named annotation (inside ACTIONet.out$annotations) or a label vector.
#' @param annotation.out Name of the updated annotations to be stored.
#' @param double.stochastic Whether to densify adjacency matrix before running label propagation (default=FALSE).
#' @param max_iter How many iterative rounds of correction/inference should be performed (default=3)
#' 
#' @return ACTIONet.out with updated annotations added to ACTIONet.out$annotations
#' 
#' @examples
#' ACTIONet.out = add.cell.annotations(ACTIONet.out, cell.labels, "input_annotations")
#' ACTIONet.out = infer.missing.cell.annotations(ACTIONet.out, "input_annotations", "updated_annotations")
infer.missing.cell.annotations <- function(ACTIONet.out, annotation.in, annotation.out, double.stochastic = FALSE, max_iter = 3) {
    if (igraph::is.igraph(ACTIONet.out)) 
        ACTIONet = ACTIONet.out else ACTIONet = ACTIONet.out$ACTIONet
    
    if (length(igraph::V(ACTIONet)) != ncol(sce)) {
        R.utils::printf("Number of cells in the input sce (%d) doesn't match the number of vertices in the ACTIONet (%d)\n", dim(sce)[2], 
            length(V(ACTIONet)))
    }
    Adj = as(igraph::get.adjacency(ACTIONet, attr = "weight"), "dgTMatrix")
    
       
    
    eps = 1e-16
    rs = Matrix::rowSums(A)
    P = sparseMatrix(i = A@i + 1, j = A@j + 1, x = A@x/rs[A@i + 1], dims = dim(A))
    if (double.stochastic == TRUE) {
        w = sqrt(Matrix::colSums(P) + eps)
        W = P %*% Matrix::Diagonal(x = 1/w, n = length(w))
        P = W %*% Matrix::t(W)
    }
    
    
	Labels = preprocess.labels(ACTIONet.out, annotation.in)			
	if(is.null(Labels)) {
		return(ACTIONet.out)
	}
	
	Annot = sort(unique(Labels))
	idx = match(Annot, Labels)
	names(Annot) = names(Labels)[idx]
    
    na.mask = is.na(Labels)
    
    i = 1
    while (sum(na.mask) > 0) {
        R.utils::printf("iter %d\n", i)
    	
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
    if(adjust.levels == T) {
		Labels = reannotate.labels(ACTIONet.out, updated.Labels)
	} else {
		Labels = updated.Labels
	}
    

   	if(! ('annotations' %in% names(ACTIONet.out)) ) {
		ACTIONet.out$annotations = list()
	}

	time.stamp = as.character(Sys.time())
	if(is.null(annotation.out)) {
		annotation.out = sprintf('InferredMissingLabels_%s', time.stamp)
	}
	
	res = list(Labels = Labels, Labels.confidence = Labels.conf, DE.profile = NULL, highlight = NULL, cells = ACTIONet.out$log$cells, time.stamp = time.stamp, type = "Missing.label.inference")
	
	cmd = ACTIONet.out$annotations[[annotation.out]] = res
	
	return(ACTIONet.out)
}

#' Uses a variant of the label propagation algorithm to correct likely noisy labels
#'
#' @param ACTIONet.out Input results to be clustered
#' (alternatively it can be the ACTIONet igraph object)
#' @param annotation.in Annotations to correct with missing values (NA) in it.
#' It can be either a named annotation (inside ACTIONet.out$annotations) or a label vector.
#' @param annotation.out Name of the updated annotations to be stored.
#' @param LFR.threshold How aggressively to update labels. The smaller the value, the more labels will be changed (default=2)
#' @param double.stochastic Whether to densify adjacency matrix before running label propagation (default=FALSE).
#' @param max_iter How many iterative rounds of correction/inference should be performed (default=3)
#' @param min.cell.fraction Annotations with less that this fraction will be removed
#' 
#' @return ACTIONet.out with updated annotations added to ACTIONet.out$annotations
#' 
#' @examples
#' ACTIONet.out = add.cell.annotations(ACTIONet.out, cell.labels, "input_annotations")
#' ACTIONet.out = correct.cell.annotations(ACTIONet.out, "input_annotations", "updated_annotations")
correct.cell.annotations <- function(ACTIONet.out, annotation.in, annotation.out, LFR.threshold = 2, double.stochastic = FALSE, max_iter = 3, adjust.levels = T, min.cell.fraction = 0.001) {
    if (igraph::is.igraph(ACTIONet.out)) 
        ACTIONet = ACTIONet.out else ACTIONet = ACTIONet.out$ACTIONet
    
    if (length(igraph::V(ACTIONet)) != ncol(sce)) {
        R.utils::printf("Number of cells in the input sce (%d) doesn't match the number of vertices in the ACTIONet (%d)\n", dim(sce)[2], 
            length(V(ACTIONet)))
    }
    Adj = as(igraph::get.adjacency(ACTIONet, attr = "weight"), "dgTMatrix")
    
	eps = 1e-16
    rs = Matrix::rowSums(A)
    P = sparseMatrix(i = A@i + 1, j = A@j + 1, x = A@x/rs[A@i + 1], dims = dim(A))
    if (double.stochastic == TRUE) {
        w = sqrt(Matrix::colSums(P) + eps)
        W = P %*% Matrix::Diagonal(x = 1/w, n = length(w))
        P = W %*% Matrix::t(W)
    }
    
    
	Labels = preprocess.labels(ACTIONet.out, annotation.in)			
	if(is.null(Labels)) {
		return(ACTIONet.out)
	}
	
	# Prunes "trivial" annotations and merges them to larger ones
	min.cells = round(min.cell.fraction*length(Labels))
	counts = table(Labels)
	Labels[Labels %in% as.numeric(names(counts)[counts < min.cells])] = NA
	mask = is.na(Labels)
	if(sum(mask) > 0) {
		Labels[mask] = NA
		ACTIONet.out = infer.missing.cell.annotations(ACTIONet.out, annotation.in = Labels, annotation.out = annotation.out)
		#Labels = preprocess.labels(ACTIONet.out, tmp.label)			
		#ACTIONet.out$annotations = ACTIONet.out$annotations[-which(names(ACTIONet.out$annotations) == tmp.label)]	
	} else {
		ACTIONet.out = add.cell.annotations(ACTIONet.out, Labels, annotation.name = annotation.out, highlight = F)
	}
	
	Labels = ACTIONet.out$annotations[[annotation.out]]$Labels
	
	Annot = sort(unique(Labels))
	idx = match(Annot, Labels)
	names(Annot) = names(Labels)[idx]
        
    
	for(i in 1:max_iter) {    
        R.utils::printf("iter %d\n", i)
    	
        new.Labels = assess.label.local.enrichment(P, Labels)
        Enrichment = new.Labels$Enrichment
        curr.enrichment = sapply(1:nrow(Enrichment), function(k) Enrichment[k, names(Labels)[k]])
        
        Diff.LFR = log2((new.Labels$Labels.confidence / curr.enrichment))
        Diff.LFR[is.na(Diff.LFR)] = 0

        Labels[Diff.LFR > LFR.threshold] = new.Labels$Labels[Diff.LFR > LFR.threshold]
        names(Labels) = Annot[Labels]
    }
	Labels.conf = new.Labels$Labels.confidence
    updated.Labels = as.numeric(Labels)
    names(updated.Labels) = names(Annot)[match(Labels, Annot)]
    
    if(adjust.levels == T) {
		Labels = reannotate.labels(ACTIONet.out, updated.Labels)
	} else {
		Labels = updated.Labels
	}
    
   	if(! ('annotations' %in% names(ACTIONet.out)) ) {
		ACTIONet.out$annotations = list()
	}

	time.stamp = as.character(Sys.time())
	if(is.null(annotation.out)) {
		annotation.out = sprintf('InferredMissingLabels_%s', time.stamp)
	}
	
	res = list(Labels = Labels, Labels.confidence = Labels.conf, DE.profile = NULL, highlight = NULL, cells = ACTIONet.out$log$cells, time.stamp = time.stamp, type = "Label.correction")
	
	cmd = ACTIONet.out$annotations[[annotation.out]] = res

	return(ACTIONet.out)
}
