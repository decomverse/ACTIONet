Chernoff.enrichment <- function(A, Ind.mat) {
	print("Running Chernoff.enrichment")
	
	if(nrow(A) != nrow(Ind.mat)) {
		print("Chernoff.Enrichment:: Number of rows do not match");
		return()		
	}
	if(max(A) > 50 & min(A) >= 0) {	
		print("log-transforming")
		A = log1p(A)	
	}
	
	X = as(Ind.mat, 'dgCMatrix')

	p_c = Matrix::colMeans(X)

	p_r = Matrix::rowMeans(X)
	rho = mean(p_r)
	
	Obs = as.matrix(Matrix::t(X) %*% A)
	Exp = as.matrix( p_c %*% ( Matrix::t(p_r) %*% A) ) / rho 
	Nu = as.matrix( (p_c^2) %*% (Matrix::t(p_r) %*% (A^2)) ) / (rho^2) 

	Lambda = Obs - Exp
	
	a = apply(A, 2, max)  
	logPvals = Lambda^2 / (2*(Nu + (p_c %*% t(a))*Lambda/3))
	
	rownames(logPvals) = colnames(X)
	colnames(logPvals) = colnames(A)

	logPvals[Lambda  < 0] = 0
	logPvals[is.na(logPvals)] = 0
	
	logPvals = logPvals / log(10)
			
	return(logPvals)
}

Chernoff.enrichment.noRowScaling <- function(A, Ind.mat) {
	if(nrow(A) != nrow(Ind.mat)) {
		message("Chernoff.enrichment.noRowScaling:: Number of rows do not match");
		return()		
	}
	if(max(A) > 50 & min(A) >= 0) {	
		A = log(1+A)	
	}
		
	X = as(Ind.mat, 'dgCMatrix')

	p_c = Matrix::colMeans(X)


	Obs = as.matrix(Matrix::t(X) %*% A)
	Exp = as.matrix(p_c %*% Matrix::t(Matrix::colSums(A)))
	Nu = as.matrix(p_c %*% Matrix::t(Matrix::colSums(A^2)))

	Lambda = Obs - Exp

	a = apply(A, 2, max)
	ones = array(1, dim = dim(Nu)[1])

	logPvals = Lambda^2 / (2*(Nu + (ones %*% t(a))*Lambda/3))

	
	rownames(logPvals) = colnames(X)
	colnames(logPvals) = colnames(A)

	logPvals[Lambda  < 0] = 0
	logPvals[is.na(logPvals)] = 0
	
	logPvals = logPvals / log(10)
			
	return(logPvals)
}

assess.continuous.autocorrelation.Cpp <- function(ACTIONet.out, variables, rand_perm = 100, num_shuffles = 10000) {
    if (is.igraph(ACTIONet.out)) 
        ACTIONet = ACTIONet.out else ACTIONet = ACTIONet.out$ACTIONet
    
    nV = length(V(ACTIONet))
    G = get.adjacency(ACTIONet, attr = "weight")

	if( (ncol(variables) == nV) & (nrow(variables) != nV))
		variables = t(variables)
		
    Enrichment = assess.computeAutocorrelation_Geary(G, as.matrix(variables), rand_perm, num_shuffles)
    
    return(Enrichment)
}

assess.continuous.autocorrelation <- function(ACTIONet.out, variables) {
	N = ncol(ACTIONet.out$build.out$ACTIONet)
	A = ACTIONet.out$build.out$ACTIONet
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



assess.categorical.autocorrelation <- function(ACTIONet.out, labels) {
    A = as(ACTIONet.out$build.out$ACTIONet, "dgTMatrix")
    
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


geneset.enrichment.archetype.v2 <- function(ACTIONet.out, genesets, top.length = 1000, min.size = 3, max.size = Inf, min.enrichment = 1, genes.subset = NULL, core = T, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH") {
       
	if(core == T) {
		if (("unification.out" %in% names(ACTIONet.out))) {
			print("Using unification.out$DE.core (merged archetypes)")
			DE.profile = as.matrix(log1p(SummarizedExperiment::assays(ACTIONet.out$unification.out$DE.core)[["significance"]]))
		} else {
			print("unification.out is not in ACTIONet.out. Please run unify.cell.states() first.")
			return()
		}
	} else {
		if (("archetype.differential.signature" %in% names(ACTIONet.out))) {
			print("Using archetype.differential.signature (all archetypes)")
			DE.profile = as.matrix(log1p(SummarizedExperiment::assays(ACTIONet.out$archetype.differential.signature)[["significance"]]))
		} else {
			print("archetype.differential.signature is not in ACTIONet.out. Please run compute.archetype.feature.specificity() first.")
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


geneset.enrichment.annotations.v2 <- function(ACTIONet.out, annotation.name, genesets, top.length = 1000, min.size = 3, max.size = Inf, genes.subset = NULL, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH", core = T) {
	idx = which((names(ACTIONet.out$annotations) == annotation.name) | (sapply(ACTIONet.out$annotations, function(X) X$annotation.name == annotation.name)))
	if(length(idx) == 0) {
		R.utils::printf('Annotation %s not found\n', annotation.name)
		return(ACTIONet.out)
	}
	
	R.utils::printf('Annotation found: name = %s, tag = %s\n', names(ACTIONet.out$annotations)[[idx]], ACTIONet.out$annotations[[idx]]$annotation.name)
	
	if(is.null(ACTIONet.out$annotations[[idx]]$DE.profile)) {
		print("Please run compute.annotations.feature.specificity() first")
		return()		
	}
	DE.profile = as.matrix(SummarizedExperiment::assays(ACTIONet.out$annotations[[idx]]$DE.profile)[["significance"]])
	
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


assess.feature.specificity <- function(sce, X, sce.data.attr = "logcounts") {	
	library(Matrix)
	print("Computing feature specificity ... ")
	
    if (is.matrix(sce) | is.sparseMatrix(sce)) {
		A = as(sce, "sparseMatrix")
    } else {        
		A = as(SummarizedExperiment::assays(sce)[[sce.data.attr]], "dgTMatrix")
    }
    
    print("Binarize matrix")
    B = A
    B@x = rep(1, length(B@x))
    
    print("Compute row stats")
    p_r = Matrix::rowMeans(A)
    p_r_bin = Matrix::rowMeans(B)
    alpha = p_r/p_r_bin
    
    print("Compute column stats")
    p_c = Matrix::colMeans(B)
    rho = mean(p_c)
    beta = Matrix::t(p_c)/rho
    
    Gamma = apply(X, 1, function(x) x * beta)
    
    print("Computing observation statistics")
    #Obs = as.matrix(Matrix::t(X %*% Matrix::t(A)))	
    Obs = as.matrix(A %*% Matrix::t(X))	
    

    print("Computing expectation statistics")
    Exp = p_r %*% t(Matrix::colSums(Gamma))
    
    
    print("Computing significance")
    Nu = p_r %*% t(Matrix::colSums(Gamma^2))
    Lambda = Obs - Exp
    
    a = apply(Gamma, 2, max)
    
    logPvals = (Lambda^2)/(2 * (Nu + Lambda %*% diag(a)/3))
    
    logPvals[Lambda < 0] = 0
    logPvals[is.na(logPvals)] = 0
    
    logPvals = logPvals/log(10)
    

    # logPvals.lower = (Lambda^2)/(2 *Nu)
    # 
    # logPvals.lower[Lambda > 0] = 0
    # logPvals.lower[is.na(logPvals.lower)] = 0
    # 
    # logPvals.lower = logPvals.lower/log(10)
    # 
    # logPvals.upper = (Lambda^2)/(2 * (Nu + Lambda %*% diag(a)/3))
    # 
    # logPvals.upper[Lambda < 0] = 0
    # logPvals.upper[is.na(logPvals.upper)] = 0
    # 
    # logPvals.upper = logPvals.upper/log(10)    

    rownames(logPvals) = rownames(Obs) = rownames(sce)
    diff.sce <- SingleCellExperiment(assays = list(significance = logPvals, profile = Obs))
    rownames(diff.sce) = rownames(sce)
    
    if(!(is.matrix(sce) | is.sparseMatrix(sce))) {
		if(class(rowRanges(sce))[[1]] == "GRanges") {
			rowRanges(diff.sce) = rowRanges(sce)
			genome(diff.sce) = genome(sce)
			values(rowRanges(sce)) = c()
		}
	} 
    
    return(diff.sce)
}


compute.archetype.feature.specificity <- function(ACTIONet.out, sce, mode = "sparse", sce.data.attr = "logcounts") {
	R.utils::printf("compute.archetype.feature.specificity:: Using mode %s\n", mode)	
    if (mode == "sparse") {
        X = t(ACTIONet.out$reconstruct.out$C_stacked)
    } else {
        X = ACTIONet.out$reconstruct.out$H_stacked
    }
    
    diff.sce = assess.feature.specificity(sce, X, sce.data.attr = sce.data.attr)

	R.utils::printf("done\n")	
    
    return(diff.sce)
}


compute.annotations.feature.specificity <- function(ACTIONet.out, sce, annotation.name, sce.data.attr = "logcounts") {
	idx = which((names(ACTIONet.out$annotations) == annotation.name))
	if(length(idx) == 0) {
		R.utils::printf('Annotation %s not found\n', annotation.name)
		return(ACTIONet.out)
	}
	
	R.utils::printf('Annotation found: name = %s, tag = %s\n', names(ACTIONet.out$annotations)[[idx]], ACTIONet.out$annotations[[idx]]$annotation.name)

	
	Labels = preprocess.labels(ACTIONet.out, ACTIONet.out$annotations[[idx]]$Labels)	
    labels = as.numeric(Labels)
    if(is.null(names(Labels))) {
    	names(Labels) = as.character(labels)
    	Annot.sorted = as.character(sort(unique(labels)))
    } else {
	    labels.char = names(Labels)
	    Annot = sort(unique(labels.char))
	    Annot.levels = labels[match(Annot, labels.char)]
	    perm = order(Annot.levels)
	    Annot.sorted = Annot[perm]
    }	
    
    X = t(sapply(Annot.sorted, function(l) as.numeric(names(Labels) == l)))
    
    
    diff.sce = assess.feature.specificity(sce, X, sce.data.attr = sce.data.attr)
    
	colnames(SummarizedExperiment::assays(diff.sce)[["significance"]]) = colnames(SummarizedExperiment::assays(diff.sce)[["profile"]]) = colnames(diff.sce) = Annot.sorted
	
	ACTIONet.out$annotations[[idx]]$DE.profile = diff.sce
	
    return(ACTIONet.out)
}


find.cluster.phenotype.associated.genes <- function(sce, phenotypes, individuals = NULL, clusters = NULL, batch = NULL, direction = "up", case.condition = 1) {
    require(scran)
    
    conditions = as.character(unique(phenotypes))
    if (length(conditions) != 2) {
        print("phenotype should be a vector with two conditions (case/control)")
    }
    names(conditions) = conditions
    
    
    if (is.null(clusters)) {
        IDX = list(1:ncol(sce))
        names(IDX) = c("all")
    } else {
        IDX = split(1:ncol(sce), clusters)
    }
    
    DE.stats = vector("list", length(IDX))
    names(DE.stats) = names(IDX)
    
    for (Type in names(IDX)) {
        R.utils::printf("Processing %s\n", Type)
        idx = IDX[[Type]]
        sub.sce = sce[, idx]
        sub.pheno = phenotypes[idx]
        
        if (is.null(individuals)) {
            
            tbl = find.cluster.markers(sub.sce, sub.pheno, direction, batch, test.type = "ttest")
            X = tbl[[conditions[[case.condition]]]][, 2:4]
            colnames(X) = c("pVal", "FDR", "logFC")
            X = X[order(X$pVal), ]
            
            DE.stats[[Type]] = X
        } else {
            sub.individuals = individuals[idx]
            
            sub.case.mask = sub.pheno == conditions[case.condition]
            sub.sce.case = sub.sce[, sub.case.mask]
            sub.sce.ctl = sub.sce[, !sub.case.mask]
            
            sub.individuals.case = sub.individuals[sub.case.mask]
            sub.individuals.ctl = sub.individuals[!sub.case.mask]
            
            IDX.case = split(1:length(sub.individuals.case), sub.individuals.case)
            IDX.case = IDX.case[sapply(IDX.case, length) > 10]
            
            IDX.ctl = split(1:length(sub.individuals.ctl), sub.individuals.ctl)
            IDX.ctl = IDX.ctl[sapply(IDX.ctl, length) > 10]
            
            R.utils::printf("\tCombining %s cases and %s controls\n", length(IDX.case), length(IDX.ctl))
            
            col.id = 1
            logFC.mat = matrix(0, nrow = nrow(sce), ncol = length(IDX.case) * length(IDX.ctl))
            Pval.mat = matrix(1, nrow = nrow(sce), ncol = length(IDX.case) * length(IDX.ctl))
            rownames(Pval.mat) = rownames(logFC.mat) = rownames(sce)
            for (case.id in 1:length(IDX.case)) {
                sub.sce1 = sub.sce.case[, IDX.case[[case.id]]]
                for (ctl.id in 1:length(IDX.ctl)) {
                  sub.sce2 = sub.sce.ctl[, IDX.ctl[[ctl.id]]]
                  
                  joint.labels = c(rep(conditions[[case.condition]], ncol(sub.sce1)), rep(setdiff(conditions, conditions[[case.condition]]), 
                    ncol(sub.sce2)))
                  joint.sce = cbind(sub.sce1, sub.sce2)
                  
                  tbl = find.cluster.markers(joint.sce, joint.labels, direction, batch, test.type = "ttest")
                  
                  X = tbl[[conditions[[case.condition]]]]
                  
                  logFC.mat[, col.id] = X[match(rownames(sce), rownames(X)), grep("logFC", colnames(X))]
                  Pval.mat[, col.id] = X[match(rownames(sce), rownames(X)), grep("val", colnames(X))]
                  
                  col.id = col.id + 1
                }
            }
            combined.logFC = Matrix::rowSums(logFC.mat)
            logFC.std = apply(logFC.mat, 1, sd)
            
            Pval.mat[Pval.mat < 1e-300] = 1e-300
            combined.pVal = ncol(Pval.mat)/Matrix::rowSums(1/Pval.mat)  # The harmonic mean p-value for combining dependent tests
            combined.FDR = p.adjust(combined.pVal, method = "fdr")
            
            X = data.frame(combined.pVal, combined.FDR, combined.logFC, logFC.std)
            colnames(X) = c("pVal", "FDR", "logFC", "logFC.std")
            rownames(X) = rownames(sce)
            
            X = X[order(X$pVal), ]
            
            DE.stats[[Type]] = X
        }
    }
    
    return(DE.stats)
}


find.archtype.binary.phenotype.associated.genes <- function(ACTIONet.out, sce, phenotypes, case.condition = NULL, individuals = NULL, direction = "up", core = T) {
	if(core == T) {
		if (("unification.out" %in% names(ACTIONet.out))) {
			print("Using unification.out$DE.core (merged archetypes)")
			H = as.matrix(ACTIONet.out$unification.out$H.core)
		} else {
			print("unification.out is not in ACTIONet.out. Please run unify.cell.states() first.")
			return()
		}
	} else {
		if (("archetype.differential.signature" %in% names(ACTIONet.out))) {
			print("Using archetype.differential.signature (all archetypes)")
			H = as.matrix(reconstruct.out$H_stacked)
		} else {
			print("archetype.differential.signature is not in ACTIONet.out. Please run compute.archetype.feature.specificity() first.")
			return()
		}
	}   

    
    if(is.null(case.condition) | !is.character(case.condition)) {
		print("You need to specificy the case condition name as a character name")
		return(ACTIONet.out)
	}
    condition.mask = as.numeric(as.character(phenotypes) == case.condition)
    
    
    
    
    if(is.null(individuals)) {
		print("No individual information is provided. Performing baseline DE")
		A = as(SummarizedExperiment::assays(sce)$logcounts, 'dgTMatrix')		

		ctl.cells.mask = (condition.mask == 0)
		case.cells.mask = (condition.mask == 1)

		ctl.mask = ctl.cells.mask[(A@j+1)]
		case.mask = case.cells.mask[(A@j+1)]

		jj = A@j+1
		ctl.j = match(jj[ctl.mask], which(ctl.cells.mask))
		case.j = match(jj[case.mask], which(case.cells.mask))
    
		A.ctl = Matrix::sparseMatrix(i = A@i[ctl.mask]+1, j = ctl.j, x = A@x[ctl.mask], dims = c(nrow(A), sum(ctl.cells.mask)))
		A.case = Matrix::sparseMatrix(i = A@i[case.mask]+1, j = case.j, x = A@x[case.mask], dims = c(nrow(A), sum(case.cells.mask)))
		H.case = H[, case.cells.mask]
		
		DE.genes = assess.feature.specificity.pairwise(A.case, A.ctl, H.case)
    }


	if(core == T) {
		cmd = sprintf("ACTIONet.out$unification.out$\'%s.DE\' = DE.genes", case.condition)
		eval(parse(text = cmd))
	} else {
		cmd = sprintf("ACTIONet.out$\'%s.DE\' = DE.genes", case.condition)
		eval(parse(text = cmd))
	}   

	return(ACTIONet.out)        
}


assess.TF.activities.Viper <- function(ACTIONet.out, inRegulon, core = T) {
    require(viper)
	if(core == T) {
		if (("unification.out" %in% names(ACTIONet.out))) {
			print("Using unification.out$DE.core (merged archetypes)")
			archetype.panel = as.matrix(log1p(SummarizedExperiment::assays(ACTIONet.out$unification.out$DE.core)[["significance"]]))
		} else {
			print("unification.out is not in ACTIONet.out. Please run unify.cell.states() first.")
			return()
		}
	} else {
		if (("archetype.differential.signature" %in% names(ACTIONet.out))) {
			print("Using archetype.differential.signature (all archetypes)")
			archetype.panel = as.matrix(log1p(SummarizedExperiment::assays(ACTIONet.out$archetype.differential.signature)[["significance"]]))
		} else {
			print("archetype.differential.signature is not in ACTIONet.out. Please run compute.archetype.feature.specificity() first.")
			return()
		}
	}   
	
    
    TF_activities <- viper::viper(eset = archetype.panel, regulon = inRegulon, nes = T, method = "none", minsize = 4, eset.filter = F)
    colnames(TF_activities) <- colnames(archetype.panel)
    
    return(TF_activities)
}


assess.archetypes.TF.activities <- function(ACTIONet.out) {
	if(!exists("ChEA3plusDB")) {
		data("ChEA3plusDB")
	}
	Enrichments = lapply(ChEA3plusDB, function(gs) enrichment = geneset.enrichment.archetype(ACTIONet.out, gs, min.size = 0, max.size = Inf))
	TF.scores = sapply(1:ncol(Enrichments[[1]]), function(j) {
		X = t(sapply(Enrichments, function(enrichment) as.numeric(enrichment[, j])))
		meta.logPval = combine.logPvals(X)
		return(meta.logPval)

	})
	rownames(TF.scores) = names(ChEA3plusDB$Enrichr)
	return(TF.scores)
}

assess.annotations.TF.activities <- function(ACTIONet.out, annotation.name) {
	if(!exists("ChEA3plusDB")) {
		data("ChEA3plusDB")
	}
	
	Enrichments = lapply(ChEA3plusDB, function(gs) enrichment = geneset.enrichment.annotations(ACTIONet.out, annotation.name, gs, min.size = 0, max.size = Inf))
	TF.scores = sapply(1:ncol(Enrichments[[1]]), function(j) {
		X = t(sapply(Enrichments, function(enrichment) as.numeric(enrichment[, j])))
		meta.logPval = combine.logPvals(X)
		return(meta.logPval)

	})
	rownames(TF.scores) = names(ChEA3plusDB$Enrichr)
	colnames(TF.scores) = colnames(Enrichments[[1]])
	
	return(TF.scores)
}



geneset.enrichment.archetype <- function(ACTIONet.out, genesets, min.size = 3, max.size = 500, min.enrichment = 1, genes.subset = NULL, core = T, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH") {
    if (is.matrix(ACTIONet.out) | is.sparseMatrix(ACTIONet.out)) {
        DE.profile = ACTIONet.out
    } else {        
		if(core == T) {
			if (("unification.out" %in% names(ACTIONet.out))) {
				print("Using unification.out$DE.core (merged archetypes)")
				DE.profile = as.matrix(log1p(SummarizedExperiment::assays(ACTIONet.out$unification.out$DE.core)[["significance"]]))
			} else {
				print("unification.out is not in ACTIONet.out. Please run unify.cell.states() first.")
				return()
			}
		} else {
			if (("archetype.differential.signature" %in% names(ACTIONet.out))) {
				print("Using archetype.differential.signature (all archetypes)")
				DE.profile = as.matrix(log1p(SummarizedExperiment::assays(ACTIONet.out$archetype.differential.signature)[["significance"]]))
			} else {
				print("archetype.differential.signature is not in ACTIONet.out. Please run compute.archetype.feature.specificity() first.")
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


geneset.enrichment.annotations <- function(ACTIONet.out, annotation.name, genesets, min.size = 3, max.size = 500, genes.subset = NULL, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH", core = T) {
	idx = which((names(ACTIONet.out$annotations) == annotation.name) | (sapply(ACTIONet.out$annotations, function(X) X$annotation.name == annotation.name)))
	if(length(idx) == 0) {
		R.utils::printf('Annotation %s not found\n', annotation.name)
		return(ACTIONet.out)
	}
	
	R.utils::printf('Annotation found: name = %s, tag = %s\n', names(ACTIONet.out$annotations)[[idx]], ACTIONet.out$annotations[[idx]]$annotation.name)
	
	if(is.null(ACTIONet.out$annotations[[idx]]$DE.profile)) {
		print("Please run compute.annotations.feature.specificity() first")
		return()		
	}
	DE.profile = as.matrix(SummarizedExperiment::assays(ACTIONet.out$annotations[[idx]]$DE.profile)[["significance"]])
	
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
