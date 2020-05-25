#' Imports data from a 10X experiment folder and constructs an `SingleCellExeriment` object
#'
#' @param input_path Folder containing input files
#' @param mtx_file Count file in Matrix Market format (default="matrix.mtx.gz")
#' @param feature_annotations Table of the same size as number of rows in the count matrix (default="features.tsv.gz")
#' @param sample_annotations Table of the same size as number of columns in the count matrix (default="barcodes.tsv.gz")
#' @param sep Column-separator used in the row/column annotations files (default='\\t')
#' @param prefilter Whether to prefilter genes/cells based on the counts
#' @param min.cell.frac.per.gene Minimum fraction of cells capturing a gene for it to be retained, if prefilter=TRUE (default=0.005)
#' @param min.genes.per.cell Minimum number of required captured genes per cell, if prefilter=TRUE (default=500)
#' 
#' @return `SingleCellExeriment` object
#' 
#' @examples
#' sce = import.sce.from.10X(input_path, prefilter=TRUE)
import.sce.from.10X <- function(input_path, mtx_file = "matrix.mtx.gz", feature_annotations = "features.tsv.gz", sample_annotations = "barcodes.tsv.gz", 
    sep = "\t", prefilter = FALSE, min.cell.frac.per.gene = 0.005, min.genes.per.cell = 500) {
    require(Matrix)
    require(scran)
    
    count.file = paste(input_path, mtx_file, sep = "/")
	if( !file.exists(count.file) ) {
		R.utils::printf("File %s not found. Consider changing `mtx_file` or `input_path` options.", count.file)
		return()
	}
	
    feature.file = paste(input_path, feature_annotations, sep = "/")
    if( !file.exists(feature.file) ) {
		R.utils::printf("File %s not found. Consider changing `mtx_file` or `input_path` options.", feature.file)
		return()
	}
	
    barcode.file = paste(input_path, sample_annotations, sep = "/")
    if( !file.exists(barcode.file) ) {
		R.utils::printf("File %s not found. Consider changing `mtx_file` or `input_path` options.", barcode.file)
		return()
	}

	print("Reading counts ...")
    counts = Matrix::readMM(count.file)
	
    
    feature_table = read.table(feature.file, header = F, sep = sep, as.is = TRUE)
    if(nrow(feature_table) == (nrow(counts) + 1)) {
		rowAnnot = DataFrame(feature_table[-1, ])
		colnames(rowAnnot) = feature_table[1, ] 
	} else {
		rowAnnot = DataFrame(feature_table)
	}
	if(ncol(rowAnnot) == 1) {
		colnames(rowAnnot) = "Gene"
	} else if(ncol(rowAnnot) == 2) {
		colnames(rowAnnot) = c("ENSEMBL", "Gene")		
	} else if(ncol(rowAnnot) == 3) {
		colnames(rowAnnot) = c("ENSEMBL", "Gene", "Feature")		
	}
    
    sample_annotations = read.table(barcode.file, header = F, sep = sep, as.is = TRUE)
    if(ncol(sample_annotations) == (ncol(counts) + 1)) {
		colAnnot = DataFrame(sample_annotations[-1, ])
		colnames(colAnnot) = sample_annotations[1, ] 
	} else {
		colAnnot = DataFrame(sample_annotations)
	}

	# Feature-barcoding
    if(ncol(rowAnnot) > 2) {
		IDX = split(1:nrow(rowAnnot), rowAnnot$V3)
		expression.counts = counts[IDX$`Gene Expression`, ]
		gene.table = rowAnnot[IDX$`Gene Expression`, 1:2]    
		IDX = IDX[!grepl("Gene Expression", names(IDX))]
	} else {
		expression.counts = counts
		gene.table = rowAnnot
		IDX = list()
	}

    rownames(expression.counts) = rowAnnot[, 1]
    colnames(expression.counts) = colAnnot[, 1]
    
    
    if (ncol(sample_annotations) > 1) {
        sce <- SingleCellExperiment(assays = list(counts = expression.counts), colData = colAnnot, rowData = rowAnnot)
    } else {
        sce <- SingleCellExperiment(assays = list(counts = expression.counts))
    }
    
    # Load additional barcoded features
    for(feature.name in names(IDX)) {
		feature.counts = counts[IDX[[feature.name]], ]
		
		row.annotations = rowAnnot[IDX[[feature.name]], ]
		
		rownames(feature.counts) = rowAnnot[, 1]
		colnames(feature.counts) = colAnnot[, 1]
		
		
		reducedDims(sce)[[feature.name]] = Matrix::t(feature.counts)
	}
	
    
    if (prefilter) {
		min.cells.per.gene = round(ncol(sce) * min.cell.frac.per.gene)
        cell.counts = Matrix::rowSums(SummarizedExperiment::assays(sce)$counts > 0)
        sce = sce[cell.counts > min.cells.per.gene, ]
        
        feature.counts = Matrix::colSums(SummarizedExperiment::assays(sce)$counts > 0)
        sce = sce[, feature.counts > min.genes.per.cell]
    }
    
    return(sce)
}


#' Constructs an `SingleCellExeriment` object from count matrix, gene names, and sample_annotations
#'
#' @param counts Matrix of counts
#' @param gene.names A list of gene names to annotate columns of the count matrix
#' @param sample_annotations Optional table of annotations for samples (columns of the count matrix)
#' @param prefilter Whether to prefilter genes/cells based on the counts
#' @param min.cell.frac.per.gene Minimum fraction of cells capturing a gene for it to be retained, if prefilter=TRUE (default=0.005)
#' @param min.genes.per.cell Minimum number of required captured genes per cell, if prefilter=TRUE (default=500)
#' 
#' @return `SingleCellExeriment` object
#' 
#' @examples
#' sce = import.sce.from.count.matrix(counts.mat, gene_names, prefilter=TRUE)
import.sce.from.count.matrix <- function(counts, gene.names, sample_annotations = NULL, prefilter = FALSE, min.cell.frac.per.gene = 0.005, min.genes.per.cell = 500) {
	if(!is.sparseMatrix(counts)) {
		counts = as(counts, "sparseMatrix")
	}	
	rownames(counts) = gene.names    
    
    if(is.null(sample_annotations)) {
		sce <- SingleCellExperiment(assays = list(counts = counts))
    } else {
        sce <- SingleCellExperiment(assays = list(counts = counts), colData = sample_annotations)
	}
	
    if (prefilter) {
		min.cells.per.gene = round(ncol(sce) * min.cell.frac.per.gene)
        cell.counts = Matrix::rowSums(SummarizedExperiment::assays(sce)$counts > 0)
        sce = sce[cell.counts > min.cells.per.gene, ]
        
        feature.counts = Matrix::colSums(SummarizedExperiment::assays(sce)$counts > 0)
        sce = sce[, feature.counts > min.genes.per.cell]
    }
    
    return(sce)
}

#' Constructs an `SingleCellExeriment` object from a full count matrix file
#'
#' @param fname Full path to the count matrix file
#' @param sep Column-separator used in count matrix file (default='\\t')
#' @param prefilter Whether to prefilter genes/cells based on the counts
#' @param min.cell.frac.per.gene Minimum fraction of cells capturing a gene for it to be retained, if prefilter=TRUE (default=0.005)
#' @param min.genes.per.cell Minimum number of required captured genes per cell, if prefilter=TRUE (default=500)
#' 
#' @return `SingleCellExeriment` object
#' 
#' @examples
#' sce = import.sce.from.table(file_name, prefilter=TRUE)
import.sce.from.table <- function(fname, sep = "\t", prefilter = FALSE, min.cell.frac.per.gene = 0.001, min.genes.per.cell = 300) {
    require(Matrix)
    require(SingleCellExperiment)
    
    counts = read.table(fname, header = TRUE, sep = sep, as.is = TRUE)
    
    if (!is.numeric(counts[1, 1])) {
        row.names = counts[, 1]
        counts = counts[, -1]
        rownames(counts) = row.names
    }
    
    counts = as(as.matrix(counts), "sparseMatrix")
    
    sce <- SingleCellExperiment(assays = list(counts = counts))
    
    if (prefilter) {
		min.cells.per.gene = round(ncol(sce) * min.cell.frac.per.gene)
        cell.counts = Matrix::rowSums(SummarizedExperiment::assays(sce)$counts > 0)
        sce = sce[cell.counts > min.cells.per.gene, ]
        
        feature.counts = Matrix::colSums(SummarizedExperiment::assays(sce)$counts > 0)
        sce = sce[, feature.counts > min.genes.per.cell]
    }
    
    return(sce)
}


#' Constructs an `SingleCellExeriment` object from a Seurat object
#' Please refer to: https://satijalab.org/seurat/v3.0/conversion_vignette.html
#'
#' @param Seurat.obj Seurat object
#' 
#' @return `SingleCellExeriment` object
#' 
#' @examples
#' sce = import.sce.from.Seurat(file_name)
import.sce.from.Seurat <- function(Seurat.obj) {
    
    sce <- Seurat::as.SingleCellExperiment(Seurat.obj)
    
    return(sce)
}

#' Constructs an `SingleCellExeriment` object from AnnData
#' Please refer to: https://satijalab.org/seurat/v3.0/conversion_vignette.html
#'
#' @param fname Path to the AnnData file
#' 
#' @return `SingleCellExeriment` object
#' 
#' @examples
#' sce = import.sce.from.AnnData(file_name)
import.sce.from.AnnData <- function(fname) {
    sce <- Seurat::as.SingleCellExperiment(Seurat::ReadH5AD(file = AnnData.fname))
    
    return(sce)
}


#' Constructs an `SingleCellExeriment` object from AnnData
#' This function depends on sceasy (https://github.com/cellgeni/sceasy) and 
#' LoomExperiment(BioC) for file conversion.
#'
#' @param fname Path to the AnnData file
#' 
#' @return `SingleCellExeriment` object
#' 
#' @examples
#' sce = import.sce.from.loom(file_name)
import.sce.from.loom <- function(fname) {
    sce <- sceasy:::readExchangeableLoom(fname)
    
    return(sce)
}




#' Constructs an `SingleCellExeriment` object from CDS format in Monocle
#'
#' @param monocle_cds CDS object
#' 
#' @return `SingleCellExeriment` object
#' 
#' @examples
#' sce = import.sce.from.CDS(monocle_cds)
import.sce.from.CDS <- function(monocle_cds) {
    
    counts = exprs(monocle_cds)
    gene_annotations = fData(monocle_cds)
    sample_annotations = pData(monocle_cds)
    
    sce <- SingleCellExperiment(assays = list(counts = counts), colData = sample_annotations, rowData = gene_annotations)
    
    return(sce)
}


#' Translates rownames() of the `sce` object for mouse and human datasets.
#' A typical use-case is when input ids are in ENSEMBL format, but user is interested to work with gene symbols.
#'
#' @param sce Input `sce` object
#' @param from Source annotation (default="ENSEMBL")
#' @param to Target annotation (default="SYMBOL")
#' @param species Either "mouse" or "human" (default="human")
#' 
#' @return `SingleCellExeriment` object with renamed rows
#' 
#' @examples
#' sce = import.sce.from.CDS(monocle_cds)
convert.sce.rownames <- function(sce, from = "ENSEMBL", to = "SYMBOL", species = "human") {
    if (species == "human") {
        library(org.Hs.eg.db)
        suppressWarnings(ids <- mapIds(org.Hs.eg.db, keys = row.names(sce), keytype = from, column = to, multiVals = "first"))
        ids[is.na(ids)] = ""
        
        sce$original_rownames = rownames(sce)
        
        rownames(sce) = ids
    } else if (species == "mouse") {
        library(org.Mm.eg.db)
        suppressWarnings(ids <- mapIds(org.Mm.eg.db, keys = row.names(sce), keytype = from, column = to, multiVals = "first"))        
        ids[is.na(ids)] = ""
        
        sce$original_rownames = rownames(sce)
        
        rownames(sce) = ids
    }
    
    return(sce)
}


preprocessDF <- function(df, drop_single_values = TRUE) {
	if(ncol(df) > 0) {
		nn = colnames(df)
		for(n in nn) {
			x = df[, n] 
			if(length(unique(x)) < 50) {
				x = factor(x, sort(unique(x)))
				df[, n] = x
			}
		}
	}
    if (ncol(df) == 0) df[['name']] <- rownames(df)
    if (drop_single_values) {
        k_singular <- sapply(df, function(x) length(unique(x)) == 1)
        if (sum(k_singular) > 0)
            warning(paste('Dropping single category variables:'),
                    paste(colnames(df)[k_singular], collapse=', '))
            df <- df[, !k_singular, drop=F]
        if (ncol(df) == 0) df[['name']] <- rownames(df)
    }
    return(df)
}

import.ace.from.legacy <- function(ACTIONet.out, sce, full.import = T) {
    ace = as(sce, "ACTIONetExperiment")
    
	ACTION.out = ACTIONet.out$ACTION.out
    pruning.out = ACTIONet.out$reconstruct.out
    G = ACTIONet.out$build.out$ACTIONet
    
	colNets(ace)$ACTIONet = G
    vis.out = ACTIONet.out$vis.out
    
    reducedDims(ace)$ACTIONet2D = vis.out$coordinates
    reducedDims(ace)$ACTIONet3D = vis.out$coordinates_3D
    colFactors(ace)$denovo_color = col2rgb(vis.out$colors)
    
    

	if(full.import == T) {
		colFactors(ace)[["H_stacked"]] = as(ACTIONet.out$reconstruct.out$H_stacked, 'sparseMatrix')
		colFactors(ace)[["C_stacked"]] = as(t(ACTIONet.out$reconstruct.out$C_stacked), 'sparseMatrix')
	}
	
	unification.out = ACTIONet.out$unification.out
	colFactors(ace)[["H_unified"]] = as(ACTIONet.out$unification.out$H.core, 'sparseMatrix')
	colFactors(ace)[["C_unified"]] = as(t(ACTIONet.out$unification.out$C.core), 'sparseMatrix')
	
	ace$assigned_archetype = ACTIONet.out$unification.out$assignments.core

	#ace$node_centrality = compute_archetype_core_centrality(G, ace$assigned_archetype)

	specificity.out = ACTIONet.out$unification.out$DE.core
	rowFactors(ace)[["H_unified_profile"]] = specificity.out[["profile"]]
	rowFactors(ace)[["H_unified_upper_significance"]] = specificity.out[["significance"]]
	
	
	# Prepare output
	trace = list(ACTION.out = ACTION.out, pruning.out = pruning.out, vis.out = vis.out, unification.out = unification.out)
    trace$log = list(genes = rownames(ace), cells = colnames(ace), time = Sys.time())
      
    out = list(ace = ace, trace = trace)
    
    return(out)
}
