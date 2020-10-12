#' Computes quality scores for cells
#'
#' @param ace ACTIONet output
#'
#' @return Updated ace object with dblscore.expression, dblscore.centrality, and dblscore.combined scores
#'
#' @examples
#' ace = compute.ACTIONet.doublet.score(ace)
#' @export
compute.ACTIONet.doublet.score <- function(ace, z.threshold = 1, combination_method = "mean") {
	scores = (ace$unified_feature_specificity)
	scores = apply(scores, 2, function(x)  exp(scale(x)))
	
	associations = as(counts(ace), "dgCMatrix")
	associations@x = rep(1, length(associations@x))
	
	
	enrichment.out = assess_enrichment(scores, associations)

	X = apply(enrichment.out$logPvals, 2, function(x) {
		x / sum(x)
	})
	h = apply(X, 1, ineq::entropy)
	
	arch.doublet_score = scale(-h)
	
	
	ACTIONet.dbl.score = -scale(ace$node_centrality)

	if(combination_method == "min") {
		combined = pmin(arch.doublet_score, ACTIONet.dbl.score)
	} else if(combination_method == "max") {
		combined = pmax(arch.doublet_score, ACTIONet.dbl.score)
	} else {
		combined = Matrix::rowMeans(cbind(arch.doublet_score, ACTIONet.dbl.score))
	}
		
	ace$dblscore.expression = arch.doublet_score
	ace$dblscore.centrality = ACTIONet.dbl.score
	ace$dblscore.combined = combined

	ace$dblscore.expression.bin = z.threshold < arch.doublet_score
	ace$dblscore.centrality.bin = z.threshold < ACTIONet.dbl.score
	ace$dblscore.combined.bin = z.threshold < combined

	
	return(ace)
}	
