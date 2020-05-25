#' Computes feature (i.e. gene) specificity scores for each cluster
#'
#' @param ace ACTIONet output object
#' @param clusters Cluster
#' @param output.slot.name Name of the output in rowFactors(ace) to store results
#' @param renormalize.logcounts.slot Name of the new assay with updated logcounts adjusted using archetypes
#' Typically it is either "logcounts" or "logcounts"
 
#' @return `ACE` object with specificity scores of each cluster added to rowFactors(ace) as a matrix with name defined by output.slot.name
#' 
#' @examples
#' ace = compute.cluster.feature.specificity(ace, ace$clusters, "cluster_specificity_scores")
compute.cluster.feature.specificity <- function(ace, clusters, output.slot.name, data.slot = "logcounts") {			
	S = assays(ace)[[data.slot]]
	
	if(is.factor(clusters)) {
		UL = levels(clusters)
	} else {
		UL = sort(unique(clusters))
	}
	lables = match(clusters, UL)
	
	# Compute gene specificity for each cluster
	specificity.out = compute_cluster_feature_specificity(S, lables)
	specificity.out = lapply(specificity.out, function(specificity.scores) {
		rownames(specificity.scores) = rownames(ace)
		colnames(specificity.scores) = paste("A", 1:ncol(specificity.scores))
		return(specificity.scores)
	})
	
	X = specificity.out[["upper_significance"]]
	colnames(X) = UL
	
	rowFactors(ace)[[output.slot.name]] = X
		
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
annotate.clusters.using.labels <- function(ace, clusters, labels) {
	
	clusters = preprocess.labels(clusters, ace)
	Labels = preprocess.labels(labels, ace)


    pop.size = length(Labels)
    pos.size = table(Labels)
    
    logPvals = sapply(sort(unique(clusters)), function(i) {
        idx = which(clusters == i)
        sample.size = length(idx)
        success.size = sapply(sort(unique(Labels)), function(i) {
        	sum(Labels[idx] == i)
        })

        logPval = HGT_tail(pop.size, pos.size, sample.size, success.size)
        
        return(logPval)
    })
    

    cl.Annot = names(clusters)[match(sort(unique(clusters)), clusters)]
    Annot = names(Labels)[match(sort(unique(Labels)), Labels)]

	colnames(logPvals) = cl.Annot
	rownames(logPvals) = Annot
    
    clusterLabels = Annot[apply(logPvals, 2, which.max)]
    
    cellLabels = match(clusterLabels[clusters], Annot)
    names(cellLabels) = clusterLabels[clusters]
    

    res = list(Labels = clusterLabels, cellLabels = cellLabels, Enrichment = logPvals)
    
    #ace$annotations[[cl.idx]]$labelEnrichment = res
    #return(ace)
    return(res)
}


#' Annotate clusters using known marker genes 
#' (It uses permutation test on cluster specificity scores)
#'
#' @param ace ACTIONet output object
#' @param marker.genes A list of lists (each a set of markers for a given cell type)
#' @param specificity.slot.name An entry in the rowFactors(ace), precomputed using compute.cluster.feature.specificity() function
#' @param rand.sample.no Number of random permutations (default=1000)
#' 
#' @return A named list: \itemize{
#' \item Labels: Inferred archetype labels
#' \item Labels.confidence: Confidence of inferred labels
#' \item Enrichment: Full enrichment matrix
#'}
#' 
#' @examples
#' data("curatedMarkers_human") # pre-packaged in ACTIONet
#' marker.genes = curatedMarkers_human$Blood$PBMC$Monaco2019.12celltypes$marker.genes
#' ace = compute.cluster.feature.specificity(ace, ace$clusters, "cluster_specificity_scores")
#' arch.annot = annotate.clusters.using.markers(ace, marker.genes = marker.genes, specificity.slot.name = "cluster_specificity_scores")
annotate.clusters.using.markers <- function(ace, marker.genes, specificity.slot.name, rand.sample.no = 1000) {
	if(is.matrix(marker.genes) | is.sparseMatrix(marker.genes)) {
			marker.genes = apply(marker.genes, 2, function(x) rownames(marker.genes)[x > 0])
	}

	if(! (specificity.slot.name %in% names(rowFactors(ace))) ) {
		message(sprintf("%s does not exist in rowFactors(ace)", specificity.slot.name))
	}

	specificity.panel = as.matrix(log1p(t(rowFactors(ace)[[specificity.slot.name]])))

    GS.names = names(marker.genes)
    if (is.null(GS.names)) {
        GS.names = sapply(1:length(GS.names), function(i) sprintf("Celltype %s", i))
    }

    markers.table = do.call(rbind, lapply(names(marker.genes), function(celltype) {
        genes = marker.genes[[celltype]]
        if (length(genes) == 0)
            return(data.frame())


        signed.count = sum(sapply(genes, function(gene) grepl("\\+$|-$", gene)))
        is.signed = signed.count > 0

        if (!is.signed) {
            df = data.frame(Gene = (genes), Direction = +1, Celltype = celltype, stringsAsFactors = F)
        } else {

            pos.genes = (as.character(sapply(genes[grepl("+", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene,
                stringr::fixed("+"), ""))))
            neg.genes = (as.character(sapply(genes[grepl("-", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene,
                stringr::fixed("-"), ""))))

            df = data.frame(Gene = c(pos.genes, neg.genes), Direction = c(rep(+1, length(pos.genes)), rep(-1, length(neg.genes))),
                Celltype = celltype, stringsAsFactors = F)
        }
    }))
    markers.table = markers.table[markers.table$Gene %in% colnames(specificity.panel), ]

    if (dim(markers.table)[1] == 0) {
        print("No markers are left")
        return()
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

    #L = names(marker.genes)
    #L.levels = L[L %in% Labels]
    #Labels = match(L, L.levels)
    #names(Labels) = L.levels
    #Labels = factor(Labels, levels = L)
    Labels.conf = apply(Z, 1, max)

    names(Labels) = rownames(specificity.panel)
    names(Labels.conf) = rownames(specificity.panel)
    rownames(Z) = rownames(specificity.panel)

    out.list = list(Labels = Labels, Labels.confidence = Labels.conf, Enrichment = Z)

    return(out.list)
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
#' data("curatedMarkers_human") # pre-packaged in ACTIONet
#' marker.genes = curatedMarkers_human$Blood$PBMC$Monaco2019.12celltypes$marker.genes
#' arch.annot = annotate.profile.using.markers(my.gene.scores.profile, marker.genes = marker.genes)
annotate.profile.using.markers <- function(feature.scores, marker.genes, rand.sample.no = 1000) {
    require(ACTIONet)
    require(igraph)
    require(Matrix)
    require(stringr)

	if(is.matrix(marker.genes) | is.sparseMatrix(marker.genes)) {
			marker.genes = apply(marker.genes, 2, function(x) rownames(marker.genes)[x > 0])
	}

	specificity.panel = feature.scores


    GS.names = names(marker.genes)
    if (is.null(GS.names)) {
        GS.names = sapply(1:length(GS.names), function(i) sprintf("Celltype %s", i))
    }

    markers.table = do.call(rbind, lapply(names(marker.genes), function(celltype) {
        genes = marker.genes[[celltype]]
        if (length(genes) == 0)
            return(data.frame())


        signed.count = sum(sapply(genes, function(gene) grepl("\\+$|-$", gene)))
        is.signed = signed.count > 0

        if (!is.signed) {
            df = data.frame(Gene = (genes), Direction = +1, Celltype = celltype, stringsAsFactors = F)
        } else {

            pos.genes = (as.character(sapply(genes[grepl("+", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene,
                stringr::fixed("+"), ""))))
            neg.genes = (as.character(sapply(genes[grepl("-", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene,
                stringr::fixed("-"), ""))))

            df = data.frame(Gene = c(pos.genes, neg.genes), Direction = c(rep(+1, length(pos.genes)), rep(-1, length(neg.genes))),
                Celltype = celltype, stringsAsFactors = F)
        }
    }))
    markers.table = markers.table[markers.table$Gene %in% colnames(specificity.panel), ]

    if (dim(markers.table)[1] == 0) {
        print("No markers are left")
        return()
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

    #L = names(marker.genes)
    #L.levels = L[L %in% Labels]
    #Labels = match(L, L.levels)
    #names(Labels) = L.levels
    #Labels = factor(Labels, levels = L)
    Labels.conf = apply(Z, 1, max)

    names(Labels) = rownames(specificity.panel)
    names(Labels.conf) = rownames(specificity.panel)
    rownames(Z) = rownames(specificity.panel)

    out.list = list(Labels = Labels, Labels.confidence = Labels.conf, Enrichment = Z)

    return(out.list)
}

