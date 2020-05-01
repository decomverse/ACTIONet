#' Compute the activity score of each TF based on the observed activity of its targets.
#' (Similar to ChEA3; it is based on a meta-analysis using 4 separate datasets)
#'
#' @param scores Matrix of genes x cell type/state
#' 
#' @return Matrix of TF x cell type/state indicating inferred TF activity scores
#' 
#' @examples
#' scores = rowFactors(ace)$H_unified_upper_significance
#' TF.scores = assess.TF.activities.from.scores(scores)
assess.TF.activities.from.scores <- function(scores) {
	if(!exists("ChEA3plusDB")) {
		data("ChEA3plusDB")
	}
	
	Enrichments = lapply(1:length(ChEA3plusDB), function(i) {
		associations = ChEA3plusDB[[i]]
		associations.mat = sapply(associations, function(gs) as.numeric(rownames(ace) %in% gs))
		Enrichment.mat = t(assess_enrichment(scores, associations.mat, L))	
	})
	
	TF.scores = sapply(1:ncol(Enrichments[[1]]), function(j) {
		X = t(sapply(Enrichments, function(enrichment) as.numeric(enrichment[, j])))
		meta.logPval = combine.logPvals(X)
		return(meta.logPval)
	
	})
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
assess.TF.activities.from.archetypes <- function(ace) {
	scores = rowFactors(ace)$H_unified_upper_significance

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
#' data("gProfilerDB_human")
#' associations = gProfilerDB_human$SYMBOL$WP
#' scores = rowFactors(ace)$H_unified_upper_significance
#' Geneset.enrichments = assess.geneset.enrichment.from.scores(scores, associations)
assess.geneset.enrichment.from.scores <- function(scores, associations, L = 1000) {
	if(is.list(associations)) {
		associations = sapply(associations, function(gs) as.numeric(rownames(scores) %in% gs))
	}
	
	Enrichment.mat = t(assess_enrichment(scores, associations, L))	
	
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
#' data("gProfilerDB_human")
#' associations = gProfilerDB_human$SYMBOL$WP
#' Geneset.enrichments = assess.geneset.enrichment.from.archetypes(ace, associations)
assess.geneset.enrichment.from.archetypes <- function(scores, associations, L = 1000) {
	scores = rowFactors(ace)$H_unified_upper_significance

	Enrichment.mat = assess.geneset.enrichment.from.scores(scores, associations, L)
	
	return(Enrichment.mat)	
}



#' Performs geneset enrichment analysis on a given set of genes
#' (It uses permutation test on cluster specificity scores)
#'
#' @param genes Set of genes
#' @param category List of functional categories to analyze (default: c("GO:BP", "REAC", "KEGG"))
#' @param organism Species name (default = "hsapiens")
#' @param top.terms Number of terms to sho
#' @param col Color of the barplot
#' 
#' @return Bar plot of the top-ranked terms
#' 
#' @examples
#' geneset.enrichment.gProfiler(my_genes)
assess.geneset.enrichment.gProfiler <- function(genes, category = c("GO:BP", "REAC", "KEGG"), organism = "hsapiens", top.terms = 10, col = "tomato") {
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
