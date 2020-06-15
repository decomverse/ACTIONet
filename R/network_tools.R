#' Uses a variant of the label propagation algorithm to infer missing labels
#'
#' @param ace Input results to be clustered
#' (alternatively it can be the ACTIONet igraph object)
#' @param annotation.in Annotations to correct with missing values (NA) in it.
#' It can be either a named annotation (inside ace$annotations) or a label vector.
#' @param annotation.out Name of the updated annotations to be stored.
#' @param double.stochastic Whether to densify adjacency matrix before running label propagation (default=FALSE).
#' @param max_iter How many iterative rounds of correction/inference should be performed (default=3)
#' @param adjust.levels Whether or not re-adjust labels at the end
#'
#' @return ace with updated annotations added to ace$annotations
#'
#' @examples
#' ace = infer.missing.cell.annotations(ace, sce$assigned_archetypes, "updated_archetype_annotations")
infer.missing.cell.annotations <- function(ace, annotation.in, annotation.out, double.stochastic = FALSE, max_iter = 3, adjust.levels = T) {
    Adj = colNets(ace)$ACTIONet
    A = as(Adj, 'dgTMatrix')

    eps = 1e-16
    rs = Matrix::rowSums(A)
    P = sparseMatrix(i = A@i + 1, j = A@j + 1, x = A@x/rs[A@i + 1], dims = dim(A))
    if (double.stochastic == TRUE) {
        w = sqrt(Matrix::colSums(P) + eps)
        W = P %*% Matrix::Diagonal(x = 1/w, n = length(w))
        P = W %*% Matrix::t(W)
    }


	Labels = preprocess.labels(annotation.in, ace)
	if(is.null(Labels)) {
		return(ace)
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
		Labels = reannotate.labels(ace, updated.Labels)
	} else {
		Labels = updated.Labels
	}

	colData(ace)[[annotation.out]] = Labels

	return(ace)
}

#' Uses a variant of the label propagation algorithm to correct likely noisy labels
#'
#' @param ace Input results to be clustered
#' (alternatively it can be the ACTIONet igraph object)
#' @param annotation.in Annotations to correct with missing values (NA) in it.
#' It can be either a named annotation (inside ace$annotations) or a label vector.
#' @param annotation.out Name of the updated annotations to be stored.
#' @param LFR.threshold How aggressively to update labels. The smaller the value, the more labels will be changed (default=2)
#' @param double.stochastic Whether to densify adjacency matrix before running label propagation (default=FALSE).
#' @param max_iter How many iterative rounds of correction/inference should be performed (default=3)
#' @param min.cell.fraction Annotations with less that this fraction will be removed
#'
#' @return ace with updated annotations added to ace$annotations
#'
#' @examples
#' ace = add.cell.annotations(ace, cell.labels, "input_annotations")
#' ace = correct.cell.annotations(ace, "input_annotations", "updated_annotations")
correct.cell.annotations <- function(ace, annotation.in, annotation.out, LFR.threshold = 2, double.stochastic = FALSE, max_iter = 3, adjust.levels = T, min.cell.fraction = 0.001) {
     Adj = colNets(ace)$ACTIONet
    A = as(Adj, 'dgTMatrix')

    eps = 1e-16
    rs = Matrix::rowSums(A)
    P = sparseMatrix(i = A@i + 1, j = A@j + 1, x = A@x/rs[A@i + 1], dims = dim(A))
    if (double.stochastic == TRUE) {
        w = sqrt(Matrix::colSums(P) + eps)
        W = P %*% Matrix::Diagonal(x = 1/w, n = length(w))
        P = W %*% Matrix::t(W)
    }


	Labels = preprocess.labels(annotation.in, ace)
	if(is.null(Labels)) {
		return(ace)
	}


	# Prunes "trivial" annotations and merges them to larger ones
	min.cells = round(min.cell.fraction*length(Labels))
	counts = table(Labels)
	Labels[Labels %in% as.numeric(names(counts)[counts < min.cells])] = NA
	mask = is.na(Labels)
	if(sum(mask) > 0) {
		Labels[mask] = NA
		ace = infer.missing.cell.annotations(ace, annotation.in = Labels, annotation.out = annotation.out)

		Labels = colData(ace)[[annotation.out]]
	}

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
		Labels = reannotate.labels(ace, updated.Labels)
	} else {
		Labels = updated.Labels
	}

	colData(ace)[[annotation.out]] = Labels

	return(ace)
}
