add.cell.annotations <- function(ACTIONet.out, cell.annotations, annotation.name = NULL, highlight = F) {	
	if(is.null(annotation.name)) {
		print("Error: You need to provide a name for the annotation")
		return(ACTIONet.out)
	}
	
	cell.annotations = preprocess.labels(ACTIONet.out, cell.annotations)
	
	if(! ('annotations' %in% names(ACTIONet.out)) ) {
		ACTIONet.out$annotations = list()
	}


	time.stamp = as.character(Sys.time())
	h = hashid_settings(salt = time.stamp, min_length = 8)
	annotation.hashtag = ENCODE(length(ACTIONet.out$annotations)+1, h)
	

	res = list(Labels = cell.annotations, Labels.confidence = NULL, DE.profile = NULL, highlight = NULL, cells = ACTIONet.out$log$cells, time.stamp = time.stamp, annotation.name = annotation.hashtag, type = "add.cell.annotations")
	
	cmd = sprintf("ACTIONet.out$annotations$\"%s\" = res", annotation.name)	
	eval(parse(text=cmd))

	R.utils::printf("Annotation %s has been added to ACTIONet.out$annotations list\n", annotation.name)

	if( highlight == T ) {
		print("Adding annotation highlights")
		ACTIONet.out = highlight.annotations(ACTIONet.out, annotation.name = annotation.name)			
	}
	
    return(ACTIONet.out)	
}

add.batch.cell.annotations <- function(ACTIONet.out, annotations.df) {	
	for(annotation.name in colnames(annotations.df)) {
		labels = as.numeric(annotations.df[, annotation.name])
		labels.char = as.character(annotations.df[, annotation.name])
		
		Annot = sort(unique(labels.char))
		Annot.levels = labels[match(Annot, labels.char)]
		perm = order(Annot.levels)
		Annot.levels = Annot.levels[perm]
		names(Annot.levels) = Annot[perm]
		
		names(labels) = names(Annot.levels)[labels]
		
		labels = preprocess.labels(ACTIONet.out, labels)
		ACTIONet.out = add.cell.annotations(ACTIONet.out, cell.annotations = labels, annotation.name = annotation.name)
	}
	
	return(ACTIONet.out)
}

extract.all.annotations <- function(ACTIONet.out) {
	if(length(ACTIONet.out$annotations) == 0) {
		return(DataFrame(row.names = ACTIONet.out$log$cells))
	}
	annotations.df = DataFrame(as.data.frame(sapply(ACTIONet.out$annotations, function(annotation) {
		labels = annotation$Labels

		labels = preprocess.labels(ACTIONet.out, labels)
		Annot = sort(unique(names(labels)))
		Annot = Annot[order(labels[match(Annot, names(labels))])]
		Labels = factor(names(labels), Annot)
	})))
	
	rownames(annotations.df) = ACTIONet.out$log$cells
	return(annotations.df)
}



annotate.archetypes.using.labels <- function(ACTIONet.out, annotation.known, core = T) {
    if(core == T) {
		profile = ACTIONet.out$unification.out$H.core		
	} else {
		profile = ACTIONet.out$reconstruct.out$H_stacked
	}


	if(length(annotation.known) > 1) {
		Labels = annotation.known
	} else {	
		idx = which(names(ACTIONet.out$annotations) == annotation.known)
		if(length(idx) == 0) {
			R.utils::printf('Error in correct.cell.labels: annotation.known "%s" not found\n', annotation.known)
			return(ACTIONet.out)
		}		
		Labels = ACTIONet.out$annotations[[idx]]$Labels    
	}
	Labels = preprocess.labels(ACTIONet.out, Labels)

    Annot = names(Labels)[match(sort(unique(Labels)), Labels)]
	
	# Using t-statistics
    Enrichment.Z = sapply(Annot, function(label) {
        mask = names(Labels) == label
        class.profile = profile[, mask]
        null.profile = profile[, !mask]
        
        N.class = sum(mask)
        N.null = sum(!mask)
        
        if( (N.class < 3) | (N.null < 3) ) {
			return(rep(0, nrow(profile)))
		}
		
        mu.class = Matrix::rowMeans(class.profile)
        mu.null = Matrix::rowMeans(null.profile)
        
        sigma_sq.class = apply(class.profile, 1, var)
        sigma_sq.null = apply(null.profile, 1, var)
        
        
        delta.mean = mu.class - mu.null
        t.stat = delta.mean/sqrt((sigma_sq.class/N.class) + (sigma_sq.null/N.null))
        return(t.stat)
    })

    
    archetypeLabels = Annot[apply(Enrichment.Z, 1, which.max)]
    Labels.confidence = apply(Enrichment.Z, 1, max)
    
    
    rownames(Enrichment.Z) = paste("A", 1:nrow(Enrichment.Z), "-", archetypeLabels, sep = "")
    
    out.list = list(Labels = archetypeLabels, Labels.confidence = Labels.confidence, Enrichment = Enrichment.Z)
    
    return(out.list)
}



annotate.clusters.using.labels <- function(ACTIONet.out, annotation.cluster, annotation.known) {
	
	cl.idx = which(names(ACTIONet.out$annotations) == annotation.cluster)
	if(length(cl.idx) == 0) {
		R.utils::printf('Error in correct.cell.labels: annotation.cluster "%s" not found\n', annotation.cluster)
		return(NULL)
	}		
	clusters = ACTIONet.out$annotations[[cl.idx]]$Labels    
	clusters = preprocess.labels(ACTIONet.out, clusters)

	if(length(annotation.known) > 1) {
		Labels = annotation.known
	} else {	
		idx = which(names(ACTIONet.out$annotations) == annotation.known)
		if(length(idx) == 0) {
			R.utils::printf('Error in correct.cell.labels: annotation.known "%s" not found\n', annotation.known)
			return(NULL)
		}		
		Labels = ACTIONet.out$annotations[[idx]]$Labels    
	}
	Labels = preprocess.labels(ACTIONet.out, Labels)


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
    
    fullLabels = sapply(sort(unique(clusters)), function(i) {
        return(sprintf("Cluster %d (%s)", i, clusterLabels[[i]]))
    })

    cellFullLabels = match(fullLabels[clusters], fullLabels)
    names(cellFullLabels) = fullLabels[clusters]
    
    res = list(Labels = clusterLabels, cellLabels = cellLabels, fullLabels = fullLabels, cellFullLabels = cellFullLabels, Enrichment = logPvals)
    
    #ACTIONet.out$annotations[[cl.idx]]$labelEnrichment = res
    #return(ACTIONet.out)
    return(res)
}


annotate.archetypes.using.markers <- function(ACTIONet.out, marker.genes, rand.sample.no = 1000, core = T) {
    require(ACTIONet)
    require(igraph)
    require(Matrix)
    require(stringr)

        if(is.matrix(marker.genes) | is.sparseMatrix(marker.genes)) {
                marker.genes = apply(marker.genes, 2, function(x) rownames(marker.genes)[x > 0])
        }

        if(core == T) {
                if (("unification.out" %in% names(ACTIONet.out))) {
                        print("Using unification.out$DE.core (merged archetypes)")
                        archetype.panel = as.matrix(log1p(t(SummarizedExperiment::assays(ACTIONet.out$unification.out$DE.core)[["significance"]])))
                } else {
                        print("unification.out is not in ACTIONet.out. Please run unify.cell.states() first.")
                        return()
                }
        } else {
                if (("archetype.differential.signature" %in% names(ACTIONet.out))) {
                        print("Using archetype.differential.signature (all archetypes)")
                        archetype.panel = as.matrix(log1p(t(SummarizedExperiment::assays(ACTIONet.out$archetype.differential.signature)[["significance"]])))
                } else {
                        print("archetype.differential.signature is not in ACTIONet.out. Please run compute.archetype.feature.specificity() first.")
                        return()
                }
        }     


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
    markers.table = markers.table[markers.table$Gene %in% colnames(archetype.panel), ]

    if (dim(markers.table)[1] == 0) {
        print("No markers are left")
        return()
    }   
    archetype.panel = archetype.panel[, markers.table$Gene]

    IDX = split(1:dim(markers.table)[1], markers.table$Celltype)

    print("Computing significance scores")
    set.seed(0)
    Z = sapply(IDX, function(idx) {
        markers = (as.character(markers.table$Gene[idx]))
        directions = markers.table$Direction[idx]
        mask = markers %in% colnames(archetype.panel)

        A = as.matrix(archetype.panel[, markers[mask]])
        sgn = as.numeric(directions[mask])
        stat = A %*% sgn

        rand.stats = sapply(1:rand.sample.no, function(i) {
            rand.samples = sample.int(dim(archetype.panel)[2], sum(mask))
            rand.A = as.matrix(archetype.panel[, rand.samples])
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

    names(Labels) = rownames(archetype.panel)
    names(Labels.conf) = rownames(archetype.panel)
    rownames(Z) = rownames(archetype.panel)

    out.list = list(Labels = Labels, Labels.confidence = Labels.conf, Enrichment = Z, archetype.panel = archetype.panel)

    return(out.list)
}

annotate.profile.using.markers <- function(profile, marker.genes, rand.sample.no = 1000, core = T) {
    require(ACTIONet)
    require(igraph)
    require(Matrix)
    require(stringr)
    
	if(is.matrix(marker.genes) | is.sparseMatrix(marker.genes)) {
		marker.genes = apply(marker.genes, 2, function(x) rownames(marker.genes)[x > 0])
	}
	    
    GS.names = names(marker.genes)
    if (is.null(GS.names)) {
        GS.names = sapply(1:length(GS.names), function(i) sprintf("Celltype %s", i))
    }
    
    profile = Matrix::t(profile)
    
    markers.table = do.call(rbind, lapply(names(marker.genes), function(celltype) {
        genes = marker.genes[[celltype]]
        if (length(genes) == 0) 
            return(data.frame())
        
        
        signed.count = sum(sapply(genes, function(gene) grepl("\\+$|-$", gene)))
        is.signed = signed.count > 0
        
        if (!is.signed) {
            df = data.frame(Gene = (genes), Direction = +1, Celltype = celltype)
        } else {
            
            pos.genes = (as.character(sapply(genes[grepl("+", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene, 
                stringr::fixed("+"), ""))))
            neg.genes = (as.character(sapply(genes[grepl("-", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene, 
                stringr::fixed("-"), ""))))
            
            df = data.frame(Gene = c(pos.genes, neg.genes), Direction = c(rep(+1, length(pos.genes)), rep(-1, length(neg.genes))), 
                Celltype = celltype)
        }
    }))
    markers.table = markers.table[markers.table$Gene %in% colnames(profile), ]
    
    if (dim(markers.table)[1] == 0) {
        print("No markers are left")
        return()
    }    
               
    IDX = split(1:dim(markers.table)[1], markers.table$Celltype)
    
    print("Computing significance scores")
    set.seed(0)
    Z = sapply(IDX, function(idx) {
        markers = (as.character(markers.table$Gene[idx]))
        directions = markers.table$Direction[idx]
        mask = markers %in% colnames(profile)
        
        A = as.matrix(profile[, markers[mask]])
        sgn = as.numeric(directions[mask])
        stat = A %*% sgn
        
        rand.stats = sapply(1:rand.sample.no, function(i) {
            rand.samples = sample.int(dim(profile)[2], sum(mask))
            rand.A = as.matrix(profile[, rand.samples])
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
    
    names(Labels) = rownames(profile)
    names(Labels.conf) = rownames(profile)
    rownames(Z) = rownames(profile)
    
    out.list = list(Labels = Labels, Labels.confidence = Labels.conf, Enrichment = Z, profile = profile)
    
    return(out.list)
}



annotate.clusters.using.markers <- function(ACTIONet.out, sce, annotation.cluster, marker.genes, rand.sample.no = 1000) {
	cl.idx = which(names(ACTIONet.out$annotations) == annotation.cluster)
	if(length(cl.idx) == 0) {
		R.utils::printf('Error in correct.cell.labels: annotation.cluster "%s" not found\n', annotation.cluster)
		return(ACTIONet.out)
	}		
	clusters = ACTIONet.out$annotations[[cl.idx]]$Labels    
	clusters = preprocess.labels(ACTIONet.out, clusters)
	
	if( is.null(ACTIONet.out$annotations[[cl.idx]]$DE.profile) ) {
		ACTIONet.out = compute.annotations.feature.specificity(ACTIONet.out, sce, annotation.cluster)
	}
	X = log1p(as.matrix(SummarizedExperiment::assays(ACTIONet.out$annotations[[cl.idx]]$DE.profile)$significance))
	# colnames(X) = ACTIONet.out$annotations$Leiden$labelEnrichment$fullLabels

    GS.names = names(marker.genes)
    if (is.null(GS.names)) {
        GS.names = sapply(1:length(GS), function(i) sprintf("Celltype %s", i))
    }
    
    markers.table = do.call(rbind, lapply(names(marker.genes), function(celltype) {
        genes = marker.genes[[celltype]]
        if (length(genes) == 0) 
            return(data.frame())
        
        
        signed.count = sum(sapply(genes, function(gene) grepl("\\+$|-$", gene)))
        is.signed = signed.count > 0
        
        if (!is.signed) {
            df = data.frame(Gene = (genes), Direction = +1, Celltype = celltype)
        } else {
            
            pos.genes = (as.character(sapply(genes[grepl("+", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene, 
                stringr::fixed("+"), ""))))
            neg.genes = (as.character(sapply(genes[grepl("-", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene, 
                stringr::fixed("-"), ""))))
            
            df = data.frame(Gene = c(pos.genes, neg.genes), Direction = c(rep(+1, length(pos.genes)), rep(-1, length(neg.genes))), 
                Celltype = celltype)
        }
    }))
    markers.table = markers.table[markers.table$Gene %in% rownames(X), ]
    if (dim(markers.table)[1] == 0) {
        print("No markers are left")
        return()
    }
    
    expression.panel = Matrix::t(X)
    
    IDX = split(1:dim(markers.table)[1], markers.table$Celltype)
    
    print("Computing significance scores")
    set.seed(0)
    Z = sapply(IDX, function(idx) {
        markers = (as.character(markers.table$Gene[idx]))
        directions = markers.table$Direction[idx]
        mask = markers %in% colnames(expression.panel)
        
        A = as.matrix(expression.panel[, markers[mask]])
        sgn = as.numeric(directions[mask])
        stat = A %*% sgn
        
        rand.stats = sapply(1:rand.sample.no, function(i) {
            rand.samples = sample.int(dim(expression.panel)[2], sum(mask))
            rand.A = as.matrix(expression.panel[, rand.samples])
            rand.stat = rand.A %*% sgn
        })
        
        cell.zscores = as.numeric((stat - apply(rand.stats, 1, mean))/apply(rand.stats, 1, sd))
        
        return(cell.zscores)
    })
    
    
    Annot = names(marker.genes)
    
    clusterLabels = Annot[apply(Z, 1, which.max)]
   
    cellLabels = match(clusterLabels[clusters], Annot)
    names(cellLabels) = clusterLabels[clusters]
    
    fullLabels = sapply(sort(unique(clusters)), function(i) {
        return(sprintf("Cluster %d (%s)", i, clusterLabels[[i]]))
    })

    cellFullLabels = match(fullLabels[clusters], fullLabels)
    names(cellFullLabels) = fullLabels[clusters]
    
    res = list(Labels = clusterLabels, cellLabels = cellLabels, fullLabels = fullLabels, cellFullLabels = cellFullLabels, Enrichment = Z)
    
    #ACTIONet.out$annotations[[cl.idx]]$markerEnrichment = res
    #return(ACTIONet.out)
    
    return(res)
}



annotate.cells.using.markers <- function(ACTIONet.out, sce, marker.genes, annotation.name = NULL, alpha_val = 0.9, rand.sample.no = 100, thread_no = 8, imputation = "PageRank") {
    require(ACTIONet)
    require(igraph)
    require(Matrix)
    require(stringr)
    
    rownames(sce) = (rownames(sce))
    
    GS.names = names(marker.genes)
    if (is.null(GS.names)) {
        GS.names = sapply(1:length(GS), function(i) sprintf("Celltype %s", i))
    }
    
    markers.table = do.call(rbind, lapply(names(marker.genes), function(celltype) {
        genes = marker.genes[[celltype]]
        if (length(genes) == 0) 
            return(data.frame())
        
        
        signed.count = sum(sapply(genes, function(gene) grepl("\\+$|-$", gene)))
        is.signed = signed.count > 0
        
        if (!is.signed) {
            df = data.frame(Gene = (genes), Direction = +1, Celltype = celltype)
        } else {
            
            pos.genes = (as.character(sapply(genes[grepl("+", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene, 
                stringr::fixed("+"), ""))))
            neg.genes = (as.character(sapply(genes[grepl("-", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene, 
                stringr::fixed("-"), ""))))
            
            df = data.frame(Gene = c(pos.genes, neg.genes), Direction = c(rep(+1, length(pos.genes)), rep(-1, length(neg.genes))), 
                Celltype = celltype)
        }
    }))
    markers.table = markers.table[markers.table$Gene %in% rownames(sce), ]
    if (dim(markers.table)[1] == 0) {
        print("No markers are left")
        return()
    }
    
    rows = match(markers.table$Gene, rownames(sce))
    if (imputation == "PageRank") {
        # PageRank-based imputation
        print("Using PageRank for imptation of marker genes")
        imputed.marker.expression = impute.genes.using.ACTIONet(ACTIONet.out, sce, markers.table$Gene, alpha_val, thread_no, prune = TRUE)
    } else {
        # PCA-based imputation
        print("Using archImpute for imptation of marker genes")
        imputed.marker.expression = t(ACTIONet.out$signature.profile[rows, ACTIONet.out$core.out$core.archs] %*% ACTIONet.out$core.out$H)
    }
    colnames(imputed.marker.expression) = (colnames(imputed.marker.expression))
    
    
    IDX = split(1:dim(markers.table)[1], markers.table$Celltype)
    
    print("Computing significance scores")
    set.seed(0)
    Z = sapply(IDX, function(idx) {
        markers = (as.character(markers.table$Gene[idx]))
        directions = markers.table$Direction[idx]
        mask = markers %in% colnames(imputed.marker.expression)
        
        A = as.matrix(imputed.marker.expression[, markers[mask]])
        sgn = as.numeric(directions[mask])
        stat = A %*% sgn
        
        rand.stats = sapply(1:rand.sample.no, function(i) {
            rand.samples = sample.int(dim(imputed.marker.expression)[2], sum(mask))
            rand.A = as.matrix(imputed.marker.expression[, rand.samples])
            rand.stat = rand.A %*% sgn
        })
        
        cell.zscores = as.numeric((stat - apply(rand.stats, 1, mean))/apply(rand.stats, 1, sd))
        
        return(cell.zscores)
    })
    
    Z[is.na(Z)] = 0
    Labels = apply(Z, 1, which.max)
    names(Labels) = colnames(Z)[Labels]
    
    Labels.conf = apply(Z, 1, max)
    
    
	if(! ('annotations' %in% names(ACTIONet.out)) ) {
		ACTIONet.out$annotations = list()
	}

	time.stamp = as.character(Sys.time())
	if(is.null(annotation.name)) {
		annotation.name = sprintf('InferredCelltypes_%s', time.stamp)
	}
	h = hashid_settings(salt = time.stamp, min_length = 8)
	annotation.hashtag = ENCODE(length(ACTIONet.out$annotations)+1, h)
	
	Labels = reannotate.labels(ACTIONet.out, Labels)
	res = list(Labels = Labels, Labels.confidence = Labels.conf, DE.profile = NULL, highlight = NULL, cells = ACTIONet.out$log$cells, time.stamp = time.stamp, annotation.name = annotation.hashtag, type = "annotate.cells.using.markers", Enrichment = Z)
	
	cmd = sprintf("ACTIONet.out$annotations$\"%s\" = res", annotation.name)	
	eval(parse(text=cmd))
	

    return(ACTIONet.out)
}


map.cell.scores.from.archetype.enrichment <- function(ACTIONet.out, Enrichment, core = T, scale = T) {
    if( core == T ) {
		cell.scores.mat = Matrix::t(ACTIONet.out$unification.out$H.core)
		
	} else {
		cell.scores.mat = Matrix::t(ACTIONet.out$reconstruct.out$H_stacked)
	}    
	
    if (nrow(Enrichment) != ncol(cell.scores.mat)) {
		print("Flipping enrichment matrix")
        Enrichment = t(Enrichment)
    }

	if(scale == T) {
		Z = (Enrichment - mean(as.numeric(Enrichment))) / sd(as.numeric(Enrichment))
		Z[Z > 3] = 3
		Z[Z < -3] = -3
		Enrichment.scaled = exp(Z)
	} else {
		Enrichment.scaled = Enrichment
	}
    
    cell.Enrichment.mat = cell.scores.mat %*% Enrichment.scaled
    colnames(cell.Enrichment.mat) = colnames(Enrichment)
    rownames(cell.Enrichment.mat) = ACTIONet.out$log$cells
    
    return(cell.Enrichment.mat)
}
	
annotate.cells.from.archetype.enrichment <- function(ACTIONet.out, Enrichment, annotation.name = NULL, scale = T, core = T) {
    cell.Enrichment.mat = map.cell.scores.from.archetype.enrichment(ACTIONet.out, Enrichment, scale = scale)
    
    cell.Labels = apply(cell.Enrichment.mat, 1, which.max)
    names(cell.Labels) = colnames(cell.Enrichment.mat)[cell.Labels]
    cell.Labels.conf = apply(cell.Enrichment.mat, 1, max)

    

	if(! ('annotations' %in% names(ACTIONet.out)) ) {
		ACTIONet.out$annotations = list()
	}

	time.stamp = as.character(Sys.time())
	if(is.null(annotation.name)) {
		annotation.name = sprintf('InferredCelltypes_%s', time.stamp)
	}
	h = hashid_settings(salt = time.stamp, min_length = 8)
	annotation.hashtag = ENCODE(length(ACTIONet.out$annotations)+1, h)
	
	cell.Labels = reannotate.labels(ACTIONet.out, cell.Labels)
	res = list(Labels = cell.Labels, Labels.confidence = cell.Labels.conf, DE.profile = NULL, highlight = NULL, cells = ACTIONet.out$log$cells, time.stamp = time.stamp, annotation.name = annotation.hashtag, type = "annotate.cells.from.archetype.enrichment", Enrichment = cell.Enrichment.mat)
	
	cmd = sprintf("ACTIONet.out$annotations$\"%s\" = res", annotation.name)	
	eval(parse(text=cmd))
	

    return(ACTIONet.out)
}

annotate.cells.from.archetypes.using.markers <- function(ACTIONet.out, marker.genes, annotation.name, rand.sample.no = 1000, confidence.threshold = 3) {
	
	arch.annot.out = annotate.archetypes.using.markers(ACTIONet.out, marker.genes, rand.sample.no = rand.sample.no, core = T)
	Enrichment = arch.annot.out$Enrichment

	confidence = apply(Enrichment, 1, max)
	archLabels = as.character(colnames(Enrichment)[apply(Enrichment, 1, which.max)])
	archLabels[confidence < confidence.threshold] = "?"

	cell.annotations = archLabels[ACTIONet.out$unification.out$assignments.core]
	
	ACTIONet.out = add.cell.annotations(ACTIONet.out, cell.annotations = cell.annotations, annotation.name = annotation.name, highlight = F)

    return(ACTIONet.out)
}


prioritize.celltypes <- function(ACTIONet.out, species = "Human", plot = T) {
    if( !('unification.out' %in% names(ACTIONet.out)) ) {
		print('Error in plot.ACTIONet.gene.view: "unification.out" is not in ACTIONet.out. Please run unify.cell.states() first.')
		return()
	}
	
	DE.profile = as.matrix(SummarizedExperiment::assays(ACTIONet.out$unification.out$DE.core)[["significance"]])
	if(tolower(species) == "human") {
		data("CellMarkerDB_human")
		marker.genes = apply(CellMarkerDB_human, 2, function(x) intersect(rownames(DE.profile), rownames(CellMarkerDB_human)[x > 0]))
	} else if(tolower(species) == "mouse"){
		data("CellMarkerDB_mouse")
		marker.genes = apply(CellMarkerDB_mouse, 2, function(x) intersect(rownames(DE.profile), rownames(CellMarkerDB_mouse)[x > 0]))
	} else {
		R.utils::printf("Unknown species: %s\n", species)
		return()
	}
	CellMarker.annot.out = annotate.archetypes.using.markers(ACTIONet.out, marker.genes, core = T)
	
	
	CellMarker.Enrichment = t(CellMarker.annot.out$Enrichment)
	
	idx = sort(unique(apply(CellMarker.Enrichment, 2, which.max)))	
	selected.celltypes = rownames(CellMarker.Enrichment)[idx]
	genesets = marker.genes[selected.celltypes]
		

	if(plot == T) {
		CellMarker.Enrichment[CellMarker.Enrichment < 0] = 0
		x = sort(CellMarker.Enrichment[CellMarker.Enrichment > 0], decreasing = T)
		nnz = round(sum(x^2)^2 / sum(x^4))
		x.threshold = x[nnz]
		mask = (apply(CellMarker.Enrichment, 1, max) > x.threshold)


		subCellMarker.Enrichment = t(CellMarker.annot.out$Enrichment[, mask])
		rows = sort(unique(apply(subCellMarker.Enrichment, 2, which.max)))
		require(RColorBrewer)
		require(ComplexHeatmap)
		gradPal = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
		ht = Heatmap(subCellMarker.Enrichment[rows, ], row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8), name = "Enrichment", row_title = "Putative Cell-types", column_title = "Cell States", col=gradPal)	
		show(ht)
	}
	
	res = list(Enrichment = CellMarker.Enrichment, suggested.markers = genesets)
	return(res)
}

highlight.annotations <- function(ACTIONet.out, annotation.name, z.threshold = -1) {
	annot.idx = which(names(ACTIONet.out$annotations) == annotation.name)
	
	if(length(annot.idx) == 0) {
		R.utils::printf('Annotation %s not found\n', annotation.name)
		return(ACTIONet.out)
	}
	
	R.utils::printf('Annotation found: name = %s, tag = %s\n', names(ACTIONet.out$annotations)[[annot.idx]], ACTIONet.out$annotations[[annot.idx]]$annotation.name)
	
	clusters = ACTIONet.out$annotations[[annot.idx]]$Labels
    
    IDX = split(1:length(clusters), clusters)
    cluster.cell.connectivity = vector("list", length(IDX))
    cluster.cell.connectivity.smoothed = vector("list", length(IDX))
    cluster.pruned.cells = vector("list", length(IDX))
    
    for (i in 1:length(IDX)) {
		print(i)
        idx = IDX[[i]]
        
        sub.ACTIONet = igraph::induced.subgraph(ACTIONet.out$ACTIONet, V(ACTIONet.out$ACTIONet)[idx])
        
        sub.cn = coreness(sub.ACTIONet)
        if(length(sub.cn) == 1) {
			cluster.cell.connectivity[[i]] = -3
			next
		}
        
		G = as(get.adjacency(sub.ACTIONet, attr = "weight"), 'sparseMatrix')
		U = matrix(sub.cn/sum(sub.cn), nrow = length(V(sub.ACTIONet)))
		pr.out = PageRank_iter(G, as(U, 'sparseMatrix'), alpha = 0.85, max_it = 5)
		
		sub.cn = pr.out*length(pr.out)
		
        
        if (mad(sub.cn) > 0) {
            z = (sub.cn - median(sub.cn))/mad(sub.cn)
        } else if (sd(sub.cn) > 0) {
            z = as.numeric(scale(sub.cn))
        } else {
            z = (as.numeric(rep(0, length(idx))))
        }
        
        cluster.pruned.cells[[i]] = idx[z < z.threshold]
        
        cluster.cell.connectivity[[i]] = z
    }
    all.cell.connectivity.scores = as.numeric(sparseVector(unlist(cluster.cell.connectivity), unlist(IDX), length(clusters)))
    all.pruned.cells = sort(unique(unlist(cluster.pruned.cells)))
    
    if (length(all.pruned.cells) > 0) {
        is.pruned = as.numeric(sparseVector(1, all.pruned.cells, length(clusters)))
    } else {
        is.pruned = rep(0, ncol(sce))
    }
    
    
    out = list(connectivity.scores = all.cell.connectivity.scores, pruned.cells = all.pruned.cells, cluster.connectivity.scores = cluster.cell.connectivity, 
        cluster.pruned.cells = cluster.pruned.cells, is.pruned = is.pruned)
    
    ACTIONet.out$annotations[[annot.idx]]$highlight = out
    
    return(ACTIONet.out)
}


suggest.markers <- function(ACTIONet.out, min.enrichment = 2, DB = "PanglaoDB", species = "human", plot = T) {
	if(DB == "PanglaoDB") {
		if(tolower(species) == "human") {
			data("PanglaoDB_human_markers")
			GS = CellMarkerDB_human
		} else if(tolower(species) == "mouse") {
			data("PanglaoDB_mouse_markers")
			GS = CellMarkerDB_mouse		
		} else {
			message("Unknown species")
			return
		}
	} else if(DB == "CellMarkerDB") {
		if(tolower(species) == "human") {
			data("CellMarkerDB_human_markers")
			GS = CellMarkerDB_human_markers
		} else if(tolower(species) == "mouse") {
			data("CellMarkerDB_mouse_markers")
			GS = CellMarkerDB_mouse_markers		
		} else {
			message("Unknown species")
			return
		}
	} else {
		message("Unknon database selected")
		return
	}
	

	En = geneset.enrichment.archetype(ACTIONet.out, GS)
	En[En < min.enrichment] = 0
	colnames(En) = paste("A", 1:ncol(En), sep="")
	W = as(MWM(En), 'dgTMatrix')
	M = En[W@i+1, W@j+1]
	if(plot == T) {
		Heatmap(M, rect_gp = gpar(col = NA), cluster_rows = F, cluster_columns = F)
	}
	out = list(markers = GS[W@i+1], match.scores = M)
	return(out)
}
