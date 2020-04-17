assess.continuous.autocorrelation <- function(ace, variables) {
	N = ncol(ace$build.out$ACTIONet)
	A = ace$build.out$ACTIONet
	degs = Matrix::colSums(A)
	L = -A
	diag(L) = degs

	W = sum(A@x)
	Cs = apply(variables, 2, function(x) {
		mu = mean(x)
		delta = x - mu
		denom = as.numeric(2*W*(t(delta) %*% delta))
		
		C = as.numeric((N-1) * t(x) %*% L %*% x / denom)
	})
	
	return(1 - Cs/2)
}



assess.categorical.autocorrelation <- function(ace, labels) {
    A = as(ace$build.out$ACTIONet, "dgTMatrix")
    
    n = length(labels)
    counts = table(labels)
    categoties = as.numeric(names(counts))
    
    pvec = as.vector(counts)/n
    
    k = length(pvec)
    
    w = A@x
    s0 = sum(w)
    s1 = sum(4 * w^2)/2
    s2 = sum((colSums(A) + rowSums(A))^2)
    
    m1.rawphi = (s0/(n * (n - 1))) * (n^2 * k * (2 - k) - n * sum(1/pvec))
    
    Q1 = sum(1/pvec)
    Q2 = sum(1/pvec^2)
    Q3 = sum(1/pvec^3)
    Q22 = sum((1/pvec) %*% t(1/pvec))
    E1 = (n^2 * Q22 - n * Q3)/(n * (n - 1))
    E2 = 4 * n^3 * Q1 - 4 * n^3 * k * Q1 + n^3 * k^2 * Q1 - 2 * (2 * n^2 * Q2 - n^2 * k * Q2) + 2 * n * Q3 - n^2 * Q22
    E2 = E2/(n * (n - 1) * (n - 2))
    
    A1 = 4 * n^4 * k^2 - 4 * n^4 * k^3 + n^4 * k^4 - (2 * n^3 * k * Q1 - n^3 * k^2 * Q1)
    A2 = 4 * n^3 * Q1 - 4 * n^3 * k * Q1 + n^3 * k^2 * Q1 - (2 * n^2 * Q2 - n^2 * k * Q2)
    Apart = A1 - 2 * A2
    
    B1 = 4 * n^3 * Q1 - 4 * n^3 * k * Q1 + n^3 * k^2 * Q1 - (2 * n^2 * Q2 - n^2 * k * Q2)
    B2 = 2 * n^2 * Q2 - n^2 * k * Q2 - n * Q3
    B3 = n^2 * Q22 - n * Q3
    Bpart = B1 - B2 - B3
    
    C1 = 2 * n^3 * k * Q1 - n^3 * k^2 * Q1 - n^2 * Q22
    C2 = 2 * n^2 * Q2 - n^2 * k * Q2 - n * Q3
    Cpart = C1 - 2 * C2
    
    E3 = (Apart - 2 * Bpart - Cpart)/(n * (n - 1) * (n - 2) * (n - 3))
    
    m2.rawphi = s1 * E1 + (s2 - 2 * s1) * E2 + (s0^2 - s2 + s1) * E3
    
    v_i = v[A@i + 1]
    v_j = v[A@j + 1]
    
    p_i = pvec[match(v_i, categoties)]
    p_j = pvec[match(v_j, categoties)]
    
    rawphi = sum(w * (2 * (v_i == v_j) - 1)/(p_i * p_j))
    
    mean.rawphi = m1.rawphi
    var.rawphi = m2.rawphi - mean.rawphi^2
    phi.z = (rawphi - mean.rawphi)/sqrt(var.rawphi)
    phi.logPval = -log10(pnorm(phi.z, lower.tail = FALSE))
    
    return(list(z = phi.z, logPval = phi.logPval, phi = rawphi))
}


geneset.enrichment.gProfiler <- function(genes, top.terms = 10, col = "tomato", organism = "hsapiens", category = c("GO:BP", "REAC", "KEGG")) {
    require(gprofiler2)
    require(ggpubr)
    
    gp.out = gprofiler2::gost(genes, ordered_query = FALSE, exclude_iea = FALSE, correction_method = "fdr", sources = category, 
        organism = organism)
    if(is.null(gp.out))
		return()
    
    terms = gp.out$result
		
    terms$logPval = -log10(terms$p_value)
    
    
    too.long = which(sapply(terms$term_name, function(x) stringr::str_length(x)) > 50)
    terms = terms[-too.long, ]
    
    terms = terms[order(terms$logPval, decreasing = TRUE), ]
    sub.terms = terms[1:min(top.terms, sum(terms$logPval > 1)), ]
    
    p = ggbarplot(sub.terms, x = "term_name", y = "logPval", sort.val = "asc", orientation = "horiz", 
        fill = col, xlab = "", ylab = "") + geom_hline(yintercept = -log10(0.05), col = "gray", lty = 2)
        
    return(p)
}


geneset.enrichment.archetype.v2 <- function(ace, genesets, top.length = 1000, min.size = 3, max.size = Inf, min.enrichment = 1, genes.subset = NULL, core = T, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH") {
       
	if(core == T) {
		if (("unification.out" %in% names(ace))) {
			print("Using unification.out$DE.core (merged archetypes)")
			DE.profile = as.matrix(log1p(rowFactors(ace)[["archetype_gene_specificity"]]))
		} else {
			print("unification.out is not in ace. Please run unify.cell.states() first.")
			return()
		}
	} else {
		if (("archetype.differential.signature" %in% names(ace))) {
			print("Using archetype.differential.signature (all archetypes)")
			DE.profile = as.matrix(log1p(SummarizedExperiment::assays(ace$archetype.differential.signature)[["significance"]]))
		} else {
			print("archetype.differential.signature is not in ace. Please run compute.archetype.feature.specificity() first.")
			return()
		}
	}                  
    
    if (is.null(rownames(DE.profile))) {
        print("Rows of the DE profile have to be named with genes.")
    }
    
    require(Matrix)
    if (is.matrix(genesets) | is.sparseMatrix(genesets)) {    
		genesets = apply(genesets, 2, function(x) intersect(rownames(DE.profile), rownames(genesets)[x > 0]))
	}
	ind.mat = as(sapply(genesets, function(gs) as.numeric(rownames(DE.profile) %in% gs)), "sparseMatrix")
	rownames(ind.mat) = rownames(DE.profile)
    
    # prune annotations
    cs = Matrix::colSums(ind.mat)
    mask = (min.size <= cs) & (cs <= max.size)
    ind.mat = as.matrix(ind.mat[, mask])


# logPvals = Chernoff.enrichment.noRowScaling(DE.profile, Ind.mat = ind.mat)
	
	
	A = apply(DE.profile, 2, function(a) sort(a, decreasing = T))[1:top.length, ]
	A.cs = apply(A, 2, function(a) cumsum(a))[1:top.length, ]
	A.cs.sq = apply(A, 2, function(a) cumsum(a^2))[1:top.length, ]
	
	
	p.col = Matrix::colMeans(ind.mat)
	
	logPvals = sapply(1:ncol(A), function(i) {
		print(i)
		perm = order(DE.profile[, i], decreasing = T)
		X = ind.mat[perm[1:top.length], ]
		
		B = as.matrix(A[, i] * X)
		B.cs = apply(B, 2, cumsum)
		
		lp = sapply(1:ncol(X), function(j) {
			x = X[, j]
			p = p.col[j]
			
			Exp = p*A.cs[, i]
			Obs = B.cs[, j]
			Lambda = Obs - Exp
			
			Nu = p*A.cs.sq[, i]
			
			scores = Lambda^2 / (2 * (Nu + max(A[, i])*Lambda/3))
			scores[Lambda < 0] = 0
			scores[is.na(scores)] = 0
			
			v = max(scores)/log(10)
	
			return(v)
		})
	})

    rownames(logPvals) = colnames(ind.mat)
    colnames(logPvals) = colnames(DE.profile)
    
    return(logPvals)
}


geneset.enrichment.annotations.v2 <- function(ace, annotation.name, genesets, top.length = 1000, min.size = 3, max.size = Inf, genes.subset = NULL, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH", core = T) {
	idx = which((names(ace$annotations) == annotation.name) | (sapply(ace$annotations, function(X) X$annotation.name == annotation.name)))
	if(length(idx) == 0) {
		R.utils::printf('Annotation %s not found\n', annotation.name)
		return(ace)
	}
	
	R.utils::printf('Annotation found: name = %s, tag = %s\n', names(ace$annotations)[[idx]], ace$annotations[[idx]]$annotation.name)
	
	if(is.null(ace$annotations[[idx]]$DE.profile)) {
		print("Please run compute.annotations.feature.specificity() first")
		return()		
	}
	DE.profile = as.matrix(SummarizedExperiment::assays(ace$annotations[[idx]]$DE.profile)[["significance"]])
	
    if (is.null(rownames(DE.profile))) {
        print("Rows of the DE profile have to be named with genes.")
    }
    
    require(Matrix)
    if (is.matrix(genesets) | is.sparseMatrix(genesets)) {    
		genesets = apply(genesets, 2, function(x) intersect(rownames(DE.profile), rownames(genesets)[x > 0]))
	}
	ind.mat = as(do.call(cbind, lapply(genesets, function(gs) as.numeric(rownames(DE.profile) %in% gs))), "sparseMatrix")
	rownames(ind.mat) = rownames(DE.profile)
    
    # prune annotations
    cs = Matrix::colSums(ind.mat)
    mask = (min.size <= cs) & (cs <= max.size)
    ind.mat = ind.mat[, mask]
    
# logPvals = Chernoff.enrichment.noRowScaling(DE.profile, Ind.mat = ind.mat)
	
	
	A = apply(DE.profile, 2, function(a) sort(a, decreasing = T))[1:top.length, ]
	A.cs = apply(A, 2, function(a) cumsum(a))[1:top.length, ]
	A.cs.sq = apply(A, 2, function(a) cumsum(a^2))[1:top.length, ]
	
	
	p.col = Matrix::colMeans(ind.mat)
	
	logPvals = sapply(1:ncol(A), function(i) {
		print(i)
		perm = order(DE.profile[, i], decreasing = T)
		X = ind.mat[perm[1:top.length], ]
		
		B = as.matrix(A[, i] * X)
		B.cs = apply(B, 2, cumsum)
		
		lp = sapply(1:ncol(X), function(j) {
			x = X[, j]
			p = p.col[j]
			
			Exp = p*A.cs[, i]
			Obs = B.cs[, j]
			Lambda = Obs - Exp
			
			Nu = p*A.cs.sq[, i]
			
			scores = Lambda^2 / (2 * (Nu + max(A[, i])*Lambda/3))
			scores[Lambda < 0] = 0
			scores[is.na(scores)] = 0
			
			v = max(scores)/log(10)
	
			return(v)
		})
	})
    
    rownames(logPvals) = colnames(ind.mat)
    colnames(logPvals) = colnames(DE.profile)
    
    return(logPvals)
}


geneset.enrichment.gene.scores.v2 <- function(gene.scores, genesets, top.length = 1000, min.size = 3, max.size = Inf, genes.subset = NULL, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH") {
	DE.profile = as.matrix(gene.scores)
	
    if (is.null(rownames(DE.profile))) {
        print("Rows of the DE profile have to be named with genes.")
    }
    
    require(Matrix)
    if (is.matrix(genesets) | is.sparseMatrix(genesets)) {    
		genesets = apply(genesets, 2, function(x) intersect(rownames(DE.profile), rownames(genesets)[x > 0]))
	}
	ind.mat = as(do.call(cbind, as.matrix(lapply(genesets, function(gs) as.numeric(rownames(DE.profile) %in% gs)))), "sparseMatrix")
	colnames(ind.mat) = names(genesets)

    # prune annotations
    cs = Matrix::colSums(ind.mat)
    mask = (min.size <= cs) & (cs <= max.size)
    if(sum(mask) == 1) {
		ind.mat = as.matrix(ind.mat[, mask])
	} else {
		ind.mat = ind.mat[, mask]    
	}
    #rs = Matrix::rowSums(ind.mat)
    #ind.mat = ind.mat[rs > 0, ]
	rownames(ind.mat) = rownames(DE.profile)
	
	
# logPvals = Chernoff.enrichment.noRowScaling(DE.profile, Ind.mat = ind.mat)
	
	
	A = apply(DE.profile, 2, function(a) sort(a, decreasing = T))[1:top.length, ]
	A.cs = apply(A, 2, function(a) cumsum(a))[1:top.length, ]
	A.cs.sq = apply(A, 2, function(a) cumsum(a^2))[1:top.length, ]
	
	
	p.col = Matrix::colMeans(ind.mat)
	
	logPvals = sapply(1:ncol(A), function(i) {
		print(i)
		perm = order(DE.profile[, i], decreasing = T)
		X = ind.mat[perm[1:top.length], ]
		
		B = as.matrix(A[, i] * X)
		B.cs = apply(B, 2, cumsum)
		
		lp = sapply(1:ncol(X), function(j) {
			x = X[, j]
			p = p.col[j]
			
			Exp = p*A.cs[, i]
			Obs = B.cs[, j]
			Lambda = Obs - Exp
			
			Nu = p*A.cs.sq[, i]
			
			scores = Lambda^2 / (2 * (Nu + max(A[, i])*Lambda/3))
			scores[Lambda < 0] = 0
			scores[is.na(scores)] = 0
			
			v = max(scores)/log(10)
	
			return(v)
		})
	})

    
    rownames(logPvals) = colnames(ind.mat)
    colnames(logPvals) = colnames(DE.profile)
        
    return(logPvals)
}


assess.archetypes.TF.activities <- function(ace) {
	if(!exists("ChEA3plusDB")) {
		data("ChEA3plusDB")
	}
	Enrichments = lapply(ChEA3plusDB, function(gs) enrichment = geneset.enrichment.archetype(ace, gs, min.size = 0, max.size = Inf))
	TF.scores = sapply(1:ncol(Enrichments[[1]]), function(j) {
		X = t(sapply(Enrichments, function(enrichment) as.numeric(enrichment[, j])))
		meta.logPval = combine.logPvals(X)
		return(meta.logPval)

	})
	rownames(TF.scores) = names(ChEA3plusDB$Enrichr)
	return(TF.scores)
}

assess.annotations.TF.activities <- function(ace, annotation.name) {
	if(!exists("ChEA3plusDB")) {
		data("ChEA3plusDB")
	}
	
	Enrichments = lapply(ChEA3plusDB, function(gs) enrichment = geneset.enrichment.annotations(ace, annotation.name, gs, min.size = 0, max.size = Inf))
	TF.scores = sapply(1:ncol(Enrichments[[1]]), function(j) {
		X = t(sapply(Enrichments, function(enrichment) as.numeric(enrichment[, j])))
		meta.logPval = combine.logPvals(X)
		return(meta.logPval)

	})
	rownames(TF.scores) = names(ChEA3plusDB$Enrichr)
	colnames(TF.scores) = colnames(Enrichments[[1]])
	
	return(TF.scores)
}



geneset.enrichment.archetype <- function(ace, genesets, min.size = 3, max.size = 500, min.enrichment = 1, genes.subset = NULL, core = T, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH") {
    if (is.matrix(ace) | is.sparseMatrix(ace)) {
        DE.profile = ace
    } else {        
		if(core == T) {
			if (("unification.out" %in% names(ace))) {
				print("Using unification.out$DE.core (merged archetypes)")
				DE.profile = as.matrix(log1p(rowFactors(ace)[["archetype_gene_specificity"]]))
			} else {
				print("unification.out is not in ace. Please run unify.cell.states() first.")
				return()
			}
		} else {
			if (("archetype.differential.signature" %in% names(ace))) {
				print("Using archetype.differential.signature (all archetypes)")
				DE.profile = as.matrix(log1p(SummarizedExperiment::assays(ace$archetype.differential.signature)[["significance"]]))
			} else {
				print("archetype.differential.signature is not in ace. Please run compute.archetype.feature.specificity() first.")
				return()
			}
		}                  
    }
    
    if (is.null(rownames(DE.profile))) {
        print("Rows of the DE profile have to be named with genes.")
    }
    
    require(Matrix)
    if (is.matrix(genesets) | is.sparseMatrix(genesets)) {    
		genesets = apply(genesets, 2, function(x) intersect(rownames(DE.profile), rownames(genesets)[x > 0]))
	}
	ind.mat = as(sapply(genesets, function(gs) as.numeric(rownames(DE.profile) %in% gs)), "sparseMatrix")
	rownames(ind.mat) = rownames(DE.profile)
    
    # prune annotations
    cs = Matrix::colSums(ind.mat)
    mask = (min.size <= cs) & (cs <= max.size)
    ind.mat = as.matrix(ind.mat[, mask])
    #rs = Matrix::rowSums(ind.mat)
    #ind.mat = ind.mat[rs > 0, ]
    
    
    # prune enrichment profile DE.profile[DE.profile < min.enrichment] = 0 rs = Matrix::rowSums(DE.profile) DE.profile = DE.profile[rs
    # > 0, ]
    
    common.genes = intersect(rownames(DE.profile), rownames(ind.mat))
    if (!is.null(genes.subset)) {
        common.genes = intersect(common.genes, genes.subset)
    }    
	filtered.genes = rownames(DE.profile)[grepl(blacklist.pattern, rownames(DE.profile), ignore.case = T)]
	common.genes = setdiff(common.genes, filtered.genes)
	
    
    idx = match(common.genes, rownames(DE.profile))
    DE.profile = DE.profile[idx, ]

    
    idx = match(common.genes, rownames(ind.mat))
    X = ind.mat[idx, ]
    
    # Normalize scores to avoid heavy-tail side-effect(s) pos.scores = DE.profile pos.scores[pos.scores < 0] = 0
    A = DE.profile #/max(DE.profile)
    
    
    p_c = Matrix::colMeans(X)
    Obs = as.matrix(Matrix::t(X) %*% A)
    
    
    Exp = as.matrix(p_c %*% Matrix::t(Matrix::colSums(A)))
    Nu = as.matrix(p_c %*% Matrix::t(Matrix::colSums(A^2)))
    
    Lambda = Obs - Exp
    
    a = apply(A, 2, max)
    ones = array(1, dim = dim(Lambda)[1])
    
    
    logPvals = Lambda^2/(2 * (Nu + (ones %*% Matrix::t(a)) * Lambda/3))
    
    rownames(logPvals) = colnames(ind.mat)
    colnames(logPvals) = colnames(DE.profile)
    
    logPvals[Lambda < 0] = 0
    logPvals[is.na(logPvals)] = 0
    
    logPvals = logPvals/log(10)
    
    
    return(logPvals)
}


geneset.enrichment.annotations <- function(ace, annotation.name, genesets, min.size = 3, max.size = 500, genes.subset = NULL, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH", core = T) {
	idx = which((names(ace$annotations) == annotation.name) | (sapply(ace$annotations, function(X) X$annotation.name == annotation.name)))
	if(length(idx) == 0) {
		R.utils::printf('Annotation %s not found\n', annotation.name)
		return(ace)
	}
	
	R.utils::printf('Annotation found: name = %s, tag = %s\n', names(ace$annotations)[[idx]], ace$annotations[[idx]]$annotation.name)
	
	if(is.null(ace$annotations[[idx]]$DE.profile)) {
		print("Please run compute.annotations.feature.specificity() first")
		return()		
	}
	DE.profile = as.matrix(SummarizedExperiment::assays(ace$annotations[[idx]]$DE.profile)[["significance"]])
	
    if (is.null(rownames(DE.profile))) {
        print("Rows of the DE profile have to be named with genes.")
    }
    
    require(Matrix)
    if (is.matrix(genesets) | is.sparseMatrix(genesets)) {    
		genesets = apply(genesets, 2, function(x) intersect(rownames(DE.profile), rownames(genesets)[x > 0]))
	}
	ind.mat = as(do.call(cbind, lapply(genesets, function(gs) as.numeric(rownames(DE.profile) %in% gs))), "sparseMatrix")
	rownames(ind.mat) = rownames(DE.profile)
    
    # prune annotations
    cs = Matrix::colSums(ind.mat)
    mask = (min.size <= cs) & (cs <= max.size)
    ind.mat = ind.mat[, mask]
    #rs = Matrix::rowSums(ind.mat)
    #ind.mat = ind.mat[rs > 0, ]
    
    
    # prune enrichment profile DE.profile[DE.profile < min.enrichment] = 0 rs = Matrix::rowSums(DE.profile) DE.profile = DE.profile[rs
    # > 0, ]
    
    common.genes = intersect(rownames(DE.profile), rownames(ind.mat))
    if (!is.null(genes.subset)) {
        common.genes = intersect(common.genes, genes.subset)
    }    
	filtered.genes = rownames(DE.profile)[grepl(blacklist.pattern, rownames(DE.profile), ignore.case = T)]
	common.genes = setdiff(common.genes, filtered.genes)
	
    
    idx = match(common.genes, rownames(DE.profile))
    DE.profile = DE.profile[idx, ]
    
    idx = match(common.genes, rownames(ind.mat))
    X = ind.mat[idx, ]
    
    # Normalize scores to avoid heavy-tail side-effect(s) pos.scores = DE.profile pos.scores[pos.scores < 0] = 0
    A = DE.profile 
    
    logPvals = Chernoff.enrichment.noRowScaling(A, X)
    
    # p_c = Matrix::colMeans(X)
    # Obs = as.matrix(Matrix::t(X) %*% A)
    # 
    # 
    # Exp = as.matrix(p_c %*% Matrix::t(Matrix::colSums(A)))
    # Nu = as.matrix(p_c %*% Matrix::t(Matrix::colSums(A^2)))
    # 
    # Lambda = Obs - Exp
    # 
    # a = apply(A, 2, max)
    # ones = array(1, dim = dim(Lambda)[1])
    # 
    # 
    # logPvals = Lambda^2/(2 * (Nu + (ones %*% Matrix::t(a)) * Lambda/3))
    # 
    # rownames(logPvals) = colnames(ind.mat)
    # colnames(logPvals) = colnames(DE.profile)
    # 
    # logPvals[Lambda < 0] = 0
    # logPvals[is.na(logPvals)] = 0
    # 
    # logPvals = logPvals/log(10)
    
    rownames(logPvals) = colnames(ind.mat)
    colnames(logPvals) = colnames(DE.profile)
    
    return(logPvals)
}


geneset.enrichment.gene.scores <- function(gene.scores, genesets, min.size = 3, max.size = 500, genes.subset = NULL, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH") {
	DE.profile = as.matrix(gene.scores)
	
    if (is.null(rownames(DE.profile))) {
        print("Rows of the DE profile have to be named with genes.")
    }
    
    require(Matrix)
    if (is.matrix(genesets) | is.sparseMatrix(genesets)) {    
		genesets = apply(genesets, 2, function(x) intersect(rownames(DE.profile), rownames(genesets)[x > 0]))
	}
	ind.mat = as(do.call(cbind, as.matrix(lapply(genesets, function(gs) as.numeric(rownames(DE.profile) %in% gs)))), "sparseMatrix")
	colnames(ind.mat) = names(genesets)

    # prune annotations
    cs = Matrix::colSums(ind.mat)
    mask = (min.size <= cs) & (cs <= max.size)
    if(sum(mask) == 1) {
		ind.mat = as.matrix(ind.mat[, mask])
	} else {
		ind.mat = ind.mat[, mask]    
	}
    #rs = Matrix::rowSums(ind.mat)
    #ind.mat = ind.mat[rs > 0, ]
	rownames(ind.mat) = rownames(DE.profile)
    
    
    # prune enrichment profile DE.profile[DE.profile < min.enrichment] = 0 rs = Matrix::rowSums(DE.profile) DE.profile = DE.profile[rs
    # > 0, ]
    
    common.genes = intersect(rownames(DE.profile), rownames(ind.mat))
    if (!is.null(genes.subset)) {
        common.genes = intersect(common.genes, genes.subset)
    }    
	filtered.genes = rownames(DE.profile)[grepl(blacklist.pattern, rownames(DE.profile), ignore.case = T)]
	common.genes = setdiff(common.genes, filtered.genes)
	
    
    idx = match(common.genes, rownames(DE.profile))
    DE.profile = DE.profile[idx, ]
    
    idx = match(common.genes, rownames(ind.mat))
    if(ncol(ind.mat) == 1) {
	    X = as.matrix(ind.mat[idx, ])
	} else {
	    X = ind.mat[idx, ]
	}
    
    # Normalize scores to avoid heavy-tail side-effect(s) pos.scores = DE.profile pos.scores[pos.scores < 0] = 0
    A = DE.profile 
    
    logPvals = Chernoff.enrichment.noRowScaling(A, X)
    
    rownames(logPvals) = colnames(ind.mat)
    colnames(logPvals) = colnames(DE.profile)
        
    return(logPvals)
}
