#' Imports data from a 10X experiment folder and constructs an `SingleCellExeriment` object
#'
#' @param input_path Folder containing input files
#' @param mtx_file Count file in Matrix Market format (default="matrix.mtx.gz")
#' @param feature_annotations Table of the same size as number of rows in the count matrix (default="features.tsv.gz")
#' @param sample_annotations Table of the same size as number of columns in the count matrix (default="barcodes.tsv.gz")
#' @param sep Column-separator used in the row/column annotations files (default='\\t')
#' @param prefilter Whether to prefilter genes/cells based on the counts.mat
#' @param min.cell.frac.per.gene Minimum fraction of cells capturing a gene for it to be retained, if prefilter=TRUE (default=0.005)
#' @param min.genes.per.cell Minimum number of required captured genes per cell, if prefilter=TRUE (default=500)
#' 
#' @return `SingleCellExeriment` object
#' 
#' @examples
#' ace = import.ace.from.10X.generic(input_path, prefilter=TRUE)
import.ace.from.10X.generic <- function(input_path, sub_path = "filtered_gene_bc_matrices", mtx_file = "matrix.mtx.gz", feature_annotations = "features.tsv.gz", sample_annotations = "barcodes.tsv.gz", 
    sep = "\t", prefilter = FALSE, min.cell.frac.per.gene = 0.005, min.genes.per.cell = 500, use.names = TRUE) {
    require(ACTIONet)
    require(S4Vectors)
    
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
    counts.mat = Matrix::readMM(count.file)
	
    
    feature_table = read.table(feature.file, header = F, sep = sep, as.is = TRUE)
    if(nrow(feature_table) == (nrow(counts.mat) + 1)) {
		rowAnnot = S4Vectors::DataFrame(feature_table[-1, ])
		colnames(rowAnnot) = feature_table[1, ] 
	} else {
		rowAnnot = S4Vectors::DataFrame(feature_table)
	}
	if(ncol(rowAnnot) == 1) {
		colnames(rowAnnot) = "Gene"
	} else if(ncol(rowAnnot) == 2) {
		colnames(rowAnnot) = c("ENSEMBL", "Gene")		
	} else if(ncol(rowAnnot) == 3) {
		colnames(rowAnnot) = c("ENSEMBL", "Gene", "Feature")		
	}
    
    sample_annotations = read.table(barcode.file, header = F, sep = sep, as.is = TRUE)
    if(ncol(sample_annotations) == (ncol(counts.mat) + 1)) {
		colAnnot = S4Vectors::DataFrame(sample_annotations[-1, ])
		colnames(colAnnot) = sample_annotations[1, ] 
	} else {
		colAnnot = S4Vectors::DataFrame(sample_annotations)
	}

	# Feature-barcoding
    if(ncol(rowAnnot) > 2) {
		IDX = split(1:nrow(rowAnnot), rowAnnot[, 3])
		expression.counts.mat = counts.mat[IDX$`Gene Expression`, ]
		gene.table = rowAnnot[IDX$`Gene Expression`, 1:2]    
		IDX = IDX[!grepl("Gene Expression", names(IDX))]
	} else {
		expression.counts.mat = counts.mat
		gene.table = rowAnnot
		IDX = list()
	}

    if(use.names == T && (ncol(rowAnnot) >= 2)) {
		rownames(expression.counts.mat) = gene.table[, 2]
	} else {
		rownames(expression.counts.mat) = gene.table[, 1]
	}
    
    colnames(expression.counts.mat) = colAnnot[, 1]
    
    
    if (ncol(sample_annotations) > 1) {
        ace <- ACTIONetExperiment(assays = list(counts = expression.counts.mat), colData = colAnnot, rowData = rowAnnot)
    } else {
        ace <- ACTIONetExperiment(assays = list(counts = expression.counts.mat))
    }
    
    # Load additional barcoded features
    for(feature.name in names(IDX)) {
		feature.counts.mat = counts.mat[IDX[[feature.name]], ]
		
		row.annotations = rowAnnot[IDX[[feature.name]], ]
		
		rownames(feature.counts.mat) = rowAnnot[IDX[[feature.name]], 1]
		colnames(feature.counts.mat) = colAnnot[, 1]
		
		colFactors(ace)[[feature.name]] = feature.counts.mat
	}
	
    
    if (prefilter) {
		min.cells.per.gene = round(ncol(ace) * min.cell.frac.per.gene)
        cell.counts.mat = Matrix::rowSums(SummarizedExperiment::assays(ace)$counts.mat > 0)
        ace = ace[cell.counts.mat > min.cells.per.gene, ]
        
        feature.counts.mat = Matrix::colSums(SummarizedExperiment::assays(ace)$counts.mat > 0)
        ace = ace[, feature.counts.mat > min.genes.per.cell]
    }
    
    return(ace)
}


#' A simple wrapper to import data from a 10X experiment folder (v2 or v3) and constructs an `SingleCellExeriment` object
#'
#' @param input_path Folder containing input files
#' @param version 2 or 3
#' @param prefilter Whether to prefilter genes/cells based on the counts.mat
#' @param min.cell.frac.per.gene Minimum fraction of cells capturing a gene for it to be retained, if prefilter=TRUE (default=0.005)
#' @param min.genes.per.cell Minimum number of required captured genes per cell, if prefilter=TRUE (default=500)
#' 
#' @return `SingleCellExeriment` object
#' 
#' @examples
#' ace = import.ace.from.10X(input_path, prefilter=TRUE)
import.ace.from.10X <- function(input_path, version = 3, prefilter = F, min.cell.frac.per.gene = 0.005, min.genes.per.cell = 500) {
	if (file.exists(paste(input_path, "genes.tsv", sep = "/"))) {
		version = 2
	}
	
	if((2 <= version) & (version < 3)) {
		mtx_file = "matrix.mtx"
		feature_annotations = "genes.tsv"
		sample_annotations = "barcodes.tsv"
	} else 
	if((3 <= version) & (version < 4)) {
		mtx_file = "matrix.mtx.gz"
		feature_annotations = "features.tsv.gz"
		sample_annotations = "barcodes.tsv.gz"
	} else {
		message("Unknown version")
	}
	
	ace = import.ace.from.10X.generic(input_path = input_path, sub_path = sub_path, mtx_file = mtx_file, feature_annotations = feature_annotations, sample_annotations = sample_annotations)		

	return(ace)
}


#' A simple wrapper to import data from a 10X experiment folder (v2 or v3) and constructs an `SingleCellExeriment` object
#'
#' @param fname Input HDF5 file
#' @param version 2 or 3
#' @param genome genome to import. Default is to import the first (if many)
#' @param prefilter Whether to prefilter genes/cells based on the counts.mat
#' @param min.cell.frac.per.gene Minimum fraction of cells capturing a gene for it to be retained, if prefilter=TRUE (default=0.005)
#' @param min.genes.per.cell Minimum number of required captured genes per cell, if prefilter=TRUE (default=500)
#' 
#' @return `SingleCellExeriment` object
#' 
#' @examples
#' ace = import.ace.from.10X.h5(fname = fname, prefilter=TRUE)
import.ace.from.10X.h5 <- function(fname, version = 3, genome = NULL, prefilter = FALSE, min.cell.frac.per.gene = 0.005, min.genes.per.cell = 500, use.names = TRUE) {
	if (!requireNamespace('hdf5r', quietly = TRUE)) {
		stop("Please install hdf5r to read HDF5 files")
	}
	if (!file.exists(fname)) {
		stop("File not found")
	}
	
	h5file <- hdf5r::H5File$new(filename = fname, mode = 'r')
	if(is.null(genome))
		genome <- names(x = h5file)[[1]]
	
	if (h5file$attr_exists("PYTABLES_FORMAT_VERSION")) {
		version = 2
	}
	
	if( (version <= 3) & (version < 4) ) {
		if (use.names) {
		feature_slot <- 'features/name'
		} else {
		feature_slot <- 'features/id'
		}
	} else if( (version <= 2) & (version < 3) ) {
		if (use.names) {
		feature_slot <- 'gene_names'
		} else {
		feature_slot <- 'genes'
		}
	}
	counts.mat <- h5file[[paste0(genome, '/data')]]
	indices <- h5file[[paste0(genome, '/indices')]]
	indptr <- h5file[[paste0(genome, '/indptr')]]
	shp <- h5file[[paste0(genome, '/shape')]]
	features <- h5file[[paste0(genome, '/', feature_slot)]][]
	barcodes <- h5file[[paste0(genome, '/barcodes')]]
	sparse.mat <- sparseMatrix(i = indices[] + 1, p = indptr[], x = as.numeric(x = counts.mat[]), dims = shp[], giveCsparse = FALSE)
	if (length(unique(features)) < length(features)) {
		features <- make.names(names = features, unique = T)
	}
	rownames(x = sparse.mat) <- features
	colnames(x = sparse.mat) <- barcodes[]
	sparse.mat <- as(object = sparse.mat, Class = 'dgCMatrix')
	
	if (h5file$exists(name = paste0(genome, '/features'))) {
		types <- h5file[[paste0(genome, '/features/feature_type')]][]
		types.unique <- unique(x = types)
		if (length(x = types.unique) > 1) {
			message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
			mats <- sapply(
			X = types.unique,
			FUN = function(x) {return(sparse.mat[which(x = types == x), ])}, simplify = FALSE, USE.NAMES = TRUE)
		}
		else {
			mats = list(sparse.mat)
		}		
	} else {
		mats = list(sparse.mat)
	}
	
	h5file$close_all()
	
	ace <- ACTIONetExperiment(assays = list(counts = mats[[1]]))
	
	# Load additional barcoded features
	for(feature.name in names(mats)[-1]) {
		colFactors(ace)[[feature.name]] = mats[[feature.name]]
	}
	
	
	if (prefilter) {
		min.cells.per.gene = round(ncol(ace) * min.cell.frac.per.gene)
		cell.counts.mat = Matrix::rowSums(SummarizedExperiment::assays(ace)$counts.mat > 0)
		ace = ace[cell.counts.mat > min.cells.per.gene, ]
		
		feature.counts.mat = Matrix::colSums(SummarizedExperiment::assays(ace)$counts.mat > 0)
		ace = ace[, feature.counts.mat > min.genes.per.cell]
	}
	
	return(ace)
}



#' Constructs an `SingleCellExeriment` object from count matrix, gene names, and sample_annotations
#'
#' @param counts.mat Matrix of counts.mat
#' @param gene.names A list of gene names to annotate columns of the count matrix
#' @param sample_annotations Optional table of annotations for samples (columns of the count matrix)
#' @param prefilter Whether to prefilter genes/cells based on the counts.mat
#' @param min.cell.frac.per.gene Minimum fraction of cells capturing a gene for it to be retained, if prefilter=TRUE (default=0.005)
#' @param min.genes.per.cell Minimum number of required captured genes per cell, if prefilter=TRUE (default=500)
#' 
#' @return `SingleCellExeriment` object
#' 
#' @examples
#' ace = import.ace.from.count.matrix(counts.mat.mat, gene_names, prefilter=TRUE)
import.ace.from.count.matrix <- function(counts.mat, gene.names, sample_annotations = NULL, prefilter = FALSE, min.cell.frac.per.gene = 0.005, min.genes.per.cell = 500) {
	if(!is.sparseMatrix(counts.mat)) {
		counts.mat = as(counts.mat, "sparseMatrix")
	}	
	rownames(counts.mat) = gene.names    
    
    if(is.null(sample_annotations)) {
		ace <- ACTIONetExperiment(assays = list(counts = counts.mat))
    } else {
        ace <- ACTIONetExperiment(assays = list(counts = counts.mat), colData = sample_annotations)
	}
	
    if (prefilter) {
		min.cells.per.gene = round(ncol(ace) * min.cell.frac.per.gene)
        cell.counts.mat = Matrix::rowSums(SummarizedExperiment::assays(ace)$counts.mat > 0)
        ace = ace[cell.counts.mat > min.cells.per.gene, ]
        
        feature.counts.mat = Matrix::colSums(SummarizedExperiment::assays(ace)$counts.mat > 0)
        ace = ace[, feature.counts.mat > min.genes.per.cell]
    }
    
    return(ace)
}

#' Constructs an `SingleCellExeriment` object from a full count matrix file
#'
#' @param fname Full path to the count matrix file
#' @param sep Column-separator used in count matrix file (default='\\t')
#' @param prefilter Whether to prefilter genes/cells based on the counts.mat
#' @param min.cell.frac.per.gene Minimum fraction of cells capturing a gene for it to be retained, if prefilter=TRUE (default=0.005)
#' @param min.genes.per.cell Minimum number of required captured genes per cell, if prefilter=TRUE (default=500)
#' 
#' @return `SingleCellExeriment` object
#' 
#' @examples
#' ace = import.ace.from.table(file_name, prefilter=TRUE)
import.ace.from.table <- function(fname, sep = "\t", prefilter = FALSE, min.cell.frac.per.gene = 0.001, min.genes.per.cell = 300) {
    require(Matrix)
    require(ACTIONetExperiment)
    
    counts.mat = read.table(fname, header = TRUE, sep = sep, as.is = TRUE)
    
    if (!is.numeric(counts.mat[1, 1])) {
        row.names = counts.mat[, 1]
        counts.mat = counts.mat[, -1]
        rownames(counts.mat) = row.names
    }
    
    counts.mat = as(as.matrix(counts.mat), "sparseMatrix")
    
    ace <- ACTIONetExperiment(assays = list(counts = counts.mat))
    
    if (prefilter) {
		min.cells.per.gene = round(ncol(ace) * min.cell.frac.per.gene)
        cell.counts.mat = Matrix::rowSums(SummarizedExperiment::assays(ace)$counts.mat > 0)
        ace = ace[cell.counts.mat > min.cells.per.gene, ]
        
        feature.counts.mat = Matrix::colSums(SummarizedExperiment::assays(ace)$counts.mat > 0)
        ace = ace[, feature.counts.mat > min.genes.per.cell]
    }
    
    return(ace)
}


#' Constructs an `SingleCellExeriment` object from a Seurat object
#' Please refer to: https://satijalab.org/seurat/v3.0/conversion_vignette.html
#'
#' @param Seurat.obj Seurat object
#' 
#' @return `SingleCellExeriment` object
#' 
#' @examples
#' ace = import.ace.from.Seurat(file_name)
import.ace.from.Seurat <- function(Seurat.obj) {
	if (!requireNamespace('Seurat', quietly = TRUE)) {
		stop("Please install Seurat to read Seurat objects")
	}	    
    ace <- as(Seurat::as.SingleCellExperiment(Seurat.obj), "ACTIONetExperiment")
    
    return(ace)
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
#' ace = import.ace.from.loom(file_name)
import.ace.from.loom <- function(fname) {
	if (!requireNamespace('sceasy', quietly = TRUE)) {
		stop("Please install sceasy to read loom files")
	}	
    ace <- as(sceasy:::readExchangeableLoom(fname), "ACTIONetExperiment")
    
    return(ace)
}




#' Constructs an `SingleCellExeriment` object from CDS format in Monocle
#'
#' @param monocle_cds CDS object
#' 
#' @return `SingleCellExeriment` object
#' 
#' @examples
#' ace = import.ace.from.CDS(monocle_cds)
import.ace.from.CDS <- function(monocle_cds) {
    
    counts.mat = exprs(monocle_cds)
    gene_annotations = fData(monocle_cds)
    sample_annotations = pData(monocle_cds)
    
    ace <- ACTIONetExperiment(assays = list(counts = counts.mat), colData = sample_annotations, rowData = gene_annotations)
    
    return(ace)
}


#' Translates rownames() of the `ace` object for mouse and human datasets.
#' A typical use-case is when input ids are in ENSEMBL format, but user is interested to work with gene symbols.
#'
#' @param ace Input `ace` object
#' @param from Source annotation (default="ENSEMBL")
#' @param to Target annotation (default="SYMBOL")
#' @param species Either "mouse" or "human" (default="human")
#' 
#' @return `SingleCellExeriment` object with renamed rows
#' 
#' @examples
#' ace = import.ace.from.CDS(monocle_cds)
convert.ace.rownames <- function(ace, from = "ENSEMBL", to = "SYMBOL", species = "human") {
    if (species == "human") {
        library(org.Hs.eg.db)
        suppressWarnings(ids <- mapIds(org.Hs.eg.db, keys = row.names(ace), keytype = from, column = to, multiVals = "first"))
        ids[is.na(ids)] = ""
        
        ace$original_rownames = rownames(ace)
        
        rownames(ace) = ids
    } else if (species == "mouse") {
        library(org.Mm.eg.db)
        suppressWarnings(ids <- mapIds(org.Mm.eg.db, keys = row.names(ace), keytype = from, column = to, multiVals = "first"))        
        ids[is.na(ids)] = ""
        
        ace$original_rownames = rownames(ace)
        
        rownames(ace) = ids
    }
    
    return(ace)
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

import.ace.from.legacy <- function(ACTIONet.out, sce, full.import = T, return.all = F) {
    ace = as(ace, "ACTIONetExperiment")
    
    if("S_r" %in% names(reducedDims(ace))) {
		reducedDims(ace)[["ACTION"]] = reducedDims(ace)[["S_r"]]
	}
	
	ACTION.out = ACTIONet.out$ACTION.out
    pruning.out = ACTIONet.out$reconstruct.out
    G = ACTIONet.out$build.out$ACTIONet
    
	colNets(ace)$ACTIONet = G
    vis.out = ACTIONet.out$vis.out
    
    reducedDims(ace)$ACTIONet2D = vis.out$coordinates
    reducedDims(ace)$ACTIONet3D = vis.out$coordinates_3D
    reducedDims(ace)$denovo_color = vis.out$colors
    
    

	if(full.import == T) {
		colFactors(ace)[["H_stacked"]] = as(ACTIONet.out$reconstruct.out$H_stacked, 'sparseMatrix')
		colFactors(ace)[["C_stacked"]] = as(t(ACTIONet.out$reconstruct.out$C_stacked), 'sparseMatrix')
	}
	
	unification.out = ACTIONet.out$unification.out
	colFactors(ace)[["H_unified"]] = as(ACTIONet.out$unification.out$H.core, 'sparseMatrix')
	colFactors(ace)[["C_unified"]] = as(t(ACTIONet.out$unification.out$C.core), 'sparseMatrix')
	
	ace$assigned_archetype = ACTIONet.out$unification.out$assignments.core
	ace$node_centrality = compute_archetype_core_centrality(colNets(ace)$ACTIONet, sample_assignments = ace$assigned_archetype)	
	
	#ace$node_centrality = compute_archetype_core_centrality(G, ace$assigned_archetype)

	specificity.out = ACTIONet.out$unification.out$DE.core
	rowFactors(ace)[["H_unified_profile"]] = specificity.out[["profile"]]
	rowFactors(ace)[["H_unified_upper_significance"]] = specificity.out[["significance"]]
	
	
	# Prepare output
	if(return.all) {
		trace = list(ACTION.out = ACTION.out, pruning.out = pruning.out, vis.out = vis.out, unification.out = unification.out)
		trace$log = list(genes = rownames(ace), cells = colnames(ace), time = Sys.time())
      
		out = list(ace = ace, trace = trace)
    
		return(out)
    } else {
		return(ace)
	}
}
