#' Compute the activity score of each TF based on the observed activity of its targets.
#' (Similar to ChEA3; it is based on a meta-analysis using 4 separate datasets)
#'
#' @param scores Matrix of genes x cell type/state
#'
#' @return Matrix of TF x cell type/state indicating inferred TF activity scores
#'
#' @examples
#' scores = rowMaps(ace)$unified_feature_specificity
#' TF.scores = assess.TF.activities.from.scores(scores)
#' @export
assess.TF.activities.from.scores <- function(scores) {
    if (!exists("ChEA3plusDB")) {
        data("ChEA3plusDB")
    }

  Enrichments = lapply(1:length(ChEA3plusDB), function(i) {
    associations = ChEA3plusDB[[i]]
    associations.mat = as(sapply(associations, function(gs) as.numeric(rownames(ace) %in% 
      gs)), "dgCMatrix")
    Enrichment = assess_enrichment(scores, associations.mat)
    Enrichment.mat = t(Enrichment[[1]])
  })


  TF.scores = Matrix::t(sapply(1:ncol(Enrichments[[1]]), function(j) {
    X = t(sapply(Enrichments, function(enrichment) as.numeric(enrichment[, 
      j])))
    meta.logPval = combine.logPvals(X)
    return(meta.logPval)
  }))
  rownames(TF.scores) = names(ChEA3plusDB$Enrichr)

  return(TF.scores)
}

#' Compute the activity score of each TF based on the observed activity of its targets.
#' (Similar to ChEA3; it is based on a meta-analysis using 4 separate datasets)
#'
#' @param ace ACTIONetExperiment (ACE) output object
#'
#' @return Matrix of TF x archetypes indicating inferred TF activity scores
#'
#' @examples
#' TF.scores = assess.TF.activities.from.archetypes(ace)
#' @export
assess.TF.activities.from.archetypes <- function(ace) {
    scores = rowMaps(ace)$unified_feature_specificity

    TF.scores = assess.TF.activities.from.scores(scores)

    return(TF.scores)
}



#' PErforms geneset enrichment analysis on arbitrary gene scores
#'
#' @param scores Gene scores
#' @param associations Either a genes x pathways membership matrix, or a set of genesets
#' @param L Maximum length of the top-ranked genes to consider
#'
#' @return Matrix pathway x cell type/states
#'
#' @examples
#' data('gProfilerDB_human')
#' associations = gProfilerDB_human$SYMBOL$WP
#' scores = rowMaps(ace)$unified_feature_specificity
#' Geneset.enrichments = assess.geneset.enrichment.from.scores(scores, associations)
#' @export
assess.geneset.enrichment.from.scores <- function(scores, associations, L = 1000) {
    if (is.list(associations)) {
        associations = sapply(associations, function(gs) as.numeric(rownames(scores) %in%
            gs))
        rownames(associations) = rownames(scores)
    }
    associations = as.matrix(associations)

    common.features = intersect(rownames(associations), rownames(scores))
    L = min(L, length(common.features))


    Enrichment.mat = t(assess_enrichment(scores[common.features, ], associations[common.features,
        ], L))

    rownames(Enrichment.mat) = colnames(associations)

    return(Enrichment.mat)
}

#' Performs geneset enrichment analysis on archetypes
#'
#' @param ace ACTIONetExperiment (ACE) output object
#' @param associations Either a genes x pathways membership matrix, or a set of genesets
#' @param L Maximum length of the top-ranked genes to consider
#'
#' @return Matrix pathway x cell type/states
#'
#' @examples
#' data('gProfilerDB_human')
#' associations = gProfilerDB_human$SYMBOL$WP
#' Geneset.enrichments = assess.geneset.enrichment.from.archetypes(ace, associations)
#' @export
assess.geneset.enrichment.from.archetypes <- function(ace, associations, min.counts = 0, specificity.slot = "unified_feature_specificity") {
	scores = rowMaps(ace)[[specificity.slot]]
	common.genes = intersect(rownames(ace), rownames(associations))

	scores = as.matrix(scores[common.genes, ])
	if(max(scores) > 100) {
		scores = log1p(scores)
	}
	associations = as(associations[common.genes, ], 'sparseMatrix')
	col.mask = (Matrix::colSums(associations) > min.counts) #& (Matrix::colSums(associations) < nrow(associations)*0.1)
	associations = associations[, col.mask]
	associations = associations[, -1]
	enrichment.out = assess_enrichment(scores, associations)

	rownames(enrichment.out$logPvals) = colnames(associations)
	rownames(enrichment.out$thresholds) = colnames(associations)
	enrichment.out$scores = scores

	return(enrichment.out)
}



#' Performs geneset enrichment analysis on archetypes
#'
#' @param ace ACTIONetExperiment (ACE) output object
#' @param associations Either a genes x pathways membership matrix, or a set of genesets
#' @param L Maximum length of the top-ranked genes to consider
#'
#' @return Matrix pathway x cell type/states
#'
#' @examples
#' data('gProfilerDB_human')
#' associations = gProfilerDB_human$SYMBOL$WP
#' Geneset.enrichments = assess.geneset.enrichment.from.archetypes(ace, associations)
#' @export
assess.peakset.enrichment.from.archetypes <- function(ace, associations, min.counts = 0, specificity.slot = "unified_feature_specificity") {
	scores = rowMaps(ace)[[specificity.slot]]
	common.genes = intersect(rownames(ace), rownames(associations))

	scores = as.matrix(scores[common.genes, ])
	if(max(scores) > 100) {
		scores = log1p(scores)
	}
	associations = as(associations[common.genes, ], 'sparseMatrix')
	col.mask = (Matrix::colSums(associations) > min.counts) #& (Matrix::colSums(associations) < nrow(associations)*0.1)
	associations = associations[, col.mask]
	associations = associations[, -1]
	enrichment.out = assess_enrichment(scores, associations)

	rownames(enrichment.out$logPvals) = colnames(associations)
	rownames(enrichment.out$thresholds) = colnames(associations)
	enrichment.out$scores = scores

	return(enrichment.out)
}


#' Performs geneset enrichment analysis on a given set of genes
#' (It uses permutation test on cluster specificity scores)
#'
#' @param genes Set of genes
#' @param category List of functional categories to analyze (default: c('GO:BP', 'REAC', 'KEGG'))
#' @param organism Species name (default = 'hsapiens')
#' @param top.terms Number of terms to sho
#' @param col Color of the barplot
#'
#' @return Bar plot of the top-ranked terms
#'
#' @examples
#' geneset.enrichment.gProfiler(my_genes)
#' @export
assess.geneset.enrichment.gProfiler <- function(genes, category = c("GO:BP", "REAC",
    "KEGG"), organism = "hsapiens", top.terms = 10, col = "tomato") {
    require(gprofiler2)
    require(ggpubr)

    gp.out = gprofiler2::gost(genes, ordered_query = FALSE, exclude_iea = FALSE,
        correction_method = "fdr", sources = category, organism = organism)
    if (is.null(gp.out))
        return()

    terms = gp.out$result

    terms$logPval = -log10(terms$p_value)


    too.long = which(sapply(terms$term_name, function(x) stringr::str_length(x)) >
        50)
    terms = terms[-too.long, ]

    terms = terms[order(terms$logPval, decreasing = TRUE), ]
    sub.terms = terms[1:min(top.terms, sum(terms$logPval > 1)), ]

    p = ggbarplot(sub.terms, x = "term_name", y = "logPval", sort.val = "asc", orientation = "horiz",
        fill = col, xlab = "", ylab = "") + geom_hline(yintercept = -log10(0.05),
        col = "gray", lty = 2)

    return(p)
}

#' @export
assess.genesets = function(arch.gs, terms.gs, N, min.pval = 1e-100, correct = T) {

    shared = t(sapply(terms.gs, function(gs1) {
        sapply(arch.gs, function(gs2) {
            nn = intersect(gs1, gs2)
        })
    }))
    colnames(shared) = names(arch.gs)

    GS.sizes = sapply(terms.gs, length)
    logPvals.out = sapply(1:ncol(shared), function(i) {
        gg = shared[, i]
        x = as.numeric(sapply(gg, length))

        n.sample = length(arch.gs[[i]])
        n.success = as.numeric(GS.sizes)

        v = rep(0, length(x))

        min.overlap = n.success * n.sample/N
        idx = which(x >= min.overlap)
        if (length(idx) == 0)
            return(v)

        # p = phyper(x[idx]-1, n.success[idx], N - n.success[idx], n.sample, lower.tail =
        # F) if(correct == T) { p = p.adjust(p, 'fdr') } p[p < min.pval] = min.pval
        # v[idx] = -log10(p)

        v[idx] = HGT_tail(population.size = N, success.count = n.success[idx], sample.size = n.sample,
            observed.success = x[idx])

        return(v)
    })
    rownames(logPvals.out) = names(terms.gs)
    colnames(logPvals.out) = names(arch.gs)

    # if(correct == T) { idx = which(logPvals.out > 0) pp = 10^(-logPvals.out[idx])
    # pp.adj = p.adjust(pp, 'BH') logPvals.out[idx] = -log10(pp.adj) }

    return(t(logPvals.out))
}

#' @export
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

compute.Geary.C <- function(L, W, x) {
  N = ncol(L)
	mu = mean(x)
	delta = x - mu
	denom = as.numeric(2*W*(t(delta) %*% delta))

	C = as.numeric((N-1) * t(x) %*% L %*% x / denom)
	
	return(1 - C/2)
}

assess.continuous.autocorrelation <- function(ace, variables, perm.no = 100) {
  set.seed(0)
  
	A = ace$ACTIONet
	degs = Matrix::colSums(A)
	L = -A
	diag(L) = degs
	W = sum(A@x)

	Z = apply(variables, 2, function(x) {
		C = compute.Geary.C(L, W, x)
		
		rand.Cs = sapply(1:perm.no, function(j) {
		  rand.x = sample(x)
  		rand.C = compute.Geary.C(L, W, rand.x)
		  
		})
		
		z = (C - mean(rand.Cs)) / sd(rand.Cs)
	})
	
	return(Z)
}


compute.phi <- function(A, labels, s0, s1, s2) {
    n = length(labels)
    counts = table(labels)
    categoties = as.numeric(names(counts))
    
    w = A@x
    
    pvec = as.vector(counts)/n
    
    k = length(pvec)
   
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
    
    v_i = labels[A@i + 1]
    v_j = labels[A@j + 1]
    
    p_i = pvec[match(v_i, categoties)]
    p_j = pvec[match(v_j, categoties)]
    
    rawphi = sum(w * (2 * (v_i == v_j) - 1)/(p_i * p_j))
    
    mean.rawphi = m1.rawphi
    var.rawphi = m2.rawphi - mean.rawphi^2
    phi.z = (rawphi - mean.rawphi)/sqrt(var.rawphi)
    phi.logPval = -log10(pnorm(phi.z, lower.tail = FALSE))
    
    return(list(z = phi.z, logPval = phi.logPval, phi = rawphi))
}

assess.categorical.autocorrelation <- function(ace, labels, perm.no = 100) {
   set.seed(0)
  
    A = ace$ACTIONet
    labels = as.numeric(as.factor(labels))
    
    w = A@x
    s0 = sum(w)
    s1 = sum(4 * w^2)/2
    # s2 = sum((colSums(A) + rowSums(A))^2)
    s2 = sum((fast_column_sums(A) + fast_row_sums(A))^2)
    
    A = as(A, "dgTMatrix")

    phi = compute.phi(A, labels, s0, s1, s2)$phi

		rand.phis = sapply(1:perm.no, function(j) {
		  rand.labels = sample(labels)
      rand.phi = compute.phi(A, rand.labels, s0, s1, s2)$phi
		})

		z = (phi - mean(rand.phis)) / sd(rand.phis)
		
		return(z)
}
