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

import_ace_from_legacy <- function(ACTIONet.out, sce, full.import = T) {
    ace = as(sce, "ACTIONetExperiment")
    
	ACTION.out = ACTIONet.out$ACTION.out
    pruning.out = ACTIONet.out$reconstruct.out
    G = ACTIONet.out$build.out$ACTIONet
    
	colNets(ace)$ACTIONet = G
    vis.out = ACTIONet.out$vis.out
    
    reducedDims(ace)$ACTIONet2D = vis.out$coordinates
    reducedDims(ace)$ACTIONet3D = vis.out$coordinates_3D
    ace$denovo_color = vis.out$colors

	if(full.import == T) {
		colFactors(ace)[["H_stacked"]] = ACTIONet.out$reconstruct.out$H_stacked
		colFactors(ace)[["C_stacked"]] = t(ACTIONet.out$reconstruct.out$C_stacked)
	}
	
	unification.out = ACTIONet.out$unification.out
	colFactors(ace)[["H_unified"]] = ACTIONet.out$unification.out$H.core
	colFactors(ace)[["C_unified"]] = t(ACTIONet.out$unification.out$C.core)
	ace$archetype_assignment = ACTIONet.out$unification.out$assignments.core

	#ace$node_centrality = compute_archetype_core_centrality(G, ace$archetype_assignment)

	specificity.out = ACTIONet.out$unification.out$DE.core
	rowFactors(ace)[["H_unified_profile"]] = specificity.out[["profile"]]
	rowFactors(ace)[["H_unified_upper_significance"]] = specificity.out[["significance"]]
	
	
	# Prepare output
	trace = list(ACTION.out = ACTION.out, pruning.out = pruning.out, vis.out = vis.out, unification.out = unification.out)
    trace$log = list(genes = rownames(ace), cells = colnames(ace), time = Sys.time())
      
    out = list(ace = ace, trace = trace)
    
    return(out)
}

AnnData2ACE <- function(inFile) {
	require(ACTIONet)
	require(reticulate)
	anndata <- reticulate::import('anndata', convert = FALSE)
	scipy <- reticulate::import('scipy', convert = TRUE)
	
	
	from = anndata$read_h5ad(inFile)
	
	meta.data <- py_to_r(from$obs)
	for (key in colnames(meta.data)) {
		if (from$obs[key]$dtype$name == "category") {
			meta.data[key] = py_to_r(from$obs[key]$astype("str"))
		}
	}

	meta.features <- py_to_r(from$var)
	for (key in colnames(meta.features)) {
		if (from$var[key]$dtype$name == "category") {
			meta.features[key] = py_to_r(from$var[key]$astype("str"))
		}
	}

	if( scipy$sparse$issparse(from$X) ) {
		data.matrix <- Matrix::sparseMatrix(
			i = as.numeric(x = from$X$indices),
			p = as.numeric(x = from$X$indptr),
			x = as.numeric(x = from$X$data),
			index1 = FALSE
		)
	} else {
		data.matrix = Matrix::t(py_to_r(from$X))
	}
	rownames(x = data.matrix) <- rownames(x = meta.features)
	colnames(x = data.matrix) <- rownames(x = meta.data)
	
	ace = ACTIONetExperiment(assays = list(logcounts = data.matrix), rowData = meta.features, colData = meta.data)

		
	obsm_keys <- toString(from$obsm$keys())
	obsm_keys <- gsub("KeysView(AxisArrays with keys: ", "", obsm_keys, fixed = TRUE)
	obsm_keys <- substr(obsm_keys, 1, nchar(obsm_keys) - 1)
	obsm_keys <- strsplit(obsm_keys, split = ", ", fixed = TRUE)[[1]]	
	for (key in obsm_keys) {
		R.utils::printf("Importing obsm: %s ... ", key)
		mat = py_to_r(from$obsm$get(key))
		if (startsWith(key, "X_")) {
			R.utils::printf("as a reducedDim()\n")
			key <- substr(key, 3, nchar(key))
			reducedDims(ace)[[key]] = mat					
		} else {
			R.utils::printf("as a colFactor()\n")
			colFactors(ace)[[key]] = Matrix::t(mat)
		}
	}	
	
	varm_keys <- toString(from$varm$keys())
	varm_keys <- gsub("KeysView(AxisArrays with keys: ", "", varm_keys, fixed = TRUE)
	varm_keys <- substr(varm_keys, 1, nchar(varm_keys) - 1)
	varm_keys <- strsplit(varm_keys, split = ", ", fixed = TRUE)[[1]]	
	for (key in varm_keys) {
		R.utils::printf("Importing varm: %s", key)
		mat = py_to_r(from$varm$get(key))
		rowFactors(ace)[[key]] = mat
	}		
	
	obsp_keys <- toString(from$obsp$keys())
	obsp_keys <- gsub("KeysView(PairwiseArrays with keys: ", "", obsp_keys, fixed = TRUE)
	obsp_keys <- substr(obsp_keys, 1, nchar(obsp_keys) - 1)
	obsp_keys <- strsplit(obsp_keys, split = ", ", fixed = TRUE)[[1]]	
	for (key in obsp_keys) {
		key <- substr(key, 1, nchar(key) - 15)

		R.utils::printf("Importing colNets: %s", key)
		mat = py_to_r(from$obsp$get(key))
		colNets(ace)[[key]] = mat
	}	
	
	
	varp_keys <- toString(from$varp$keys())
	varp_keys <- gsub("KeysView(PairwiseArrays with keys: ", "", varp_keys, fixed = TRUE)
	varp_keys <- substr(varp_keys, 1, nchar(varp_keys) - 1)
	varp_keys <- strsplit(varp_keys, split = ", ", fixed = TRUE)[[1]]	
	for (key in varp_keys) {
		R.utils::printf("Importing colNets: %s", key)
		mat = py_to_r(from$varp$get(key))
		rowNets(ace)[[key]] = mat
	}	
	
	return(ace)
}

ACE2AnnData <- function(ace, outFile = NULL, main_layer = 'logcounts', transfer_layers = c()) {
	
    assay_names <- SummarizedExperiment::assayNames(ace)
    main_layer <- match.arg(main_layer, assay_names)
    transfer_layers <- transfer_layers[transfer_layers %in% assay_names]
    transfer_layers <- transfer_layers[transfer_layers != main_layer]

    X <- SummarizedExperiment::assay(ace, main_layer)

    obs <- preprocessDF(as.data.frame(SummarizedExperiment::colData(ace)))
    rownames(obs) = make.names(colnames(ace), unique = TRUE);
	
	
    var <- preprocessDF(as.data.frame(SummarizedExperiment::rowData(ace)))
	rownames(var) = make.names(rownames(ace), unique = TRUE)
	
    obsm <- NULL
    reductions <- SingleCellExperiment::reducedDimNames(ace)
    if (length(reductions) > 0) {
        obsm <- sapply(
            reductions,
            function(name) as.matrix(
                    SingleCellExperiment::reducedDim(ace, type=name)),
            simplify = FALSE
        )
        names(obsm) <- paste0(
            'X_', tolower(SingleCellExperiment::reducedDimNames(ace)))
    }
	
    Fn.o = names(ACTIONet::colFactors(ace))
    if (length(Fn.o) > 0) {
        obsm.ext <- sapply(Fn.o, function(name) Matrix::t(as.matrix(ACTIONet::colFactors(ace)[[name]])), simplify = FALSE)
        names(obsm.ext) <- Fn.o
        obsm = c(obsm, obsm.ext)
    }

	varm = NULL
	Fn.v = names(ACTIONet::rowFactors(ace))
	if (length(Fn.v) > 0) {
		varm <- sapply(Fn.v, function(name) as.matrix(ACTIONet::rowFactors(ace)[[name]]), simplify = FALSE)
		names(varm) <- Fn.v
	} 	

	
	varp = NULL
	Nn.v = names(ACTIONet::rowNets(ace))
	if (length(Nn.v) > 0) {
		varp <- sapply(Nn.v, function(name) as(ACTIONet::rowNets(ace)[[name]], 'sparseMatrix'), simplify = FALSE)
		names(obsp) <- paste0(Nn.v, '_connectivities');
	} 	
	
	obsp = NULL
	Nn.o = names(ACTIONet::colNets(ace))
	if (length(Nn.o) > 0) {
		obsp <- sapply(Nn.o, function(name) as(ACTIONet::colNets(ace)[[name]], 'sparseMatrix'), simplify = FALSE)
		names(varp) <- paste0(Nn.o, '_connectivities');
	} 		
	
    layers <- list()
    for (layer in transfer_layers) {
        mat <- SummarizedExperiment::assay(ace, layer)
        if (all(dim(mat) == dim(X))) layers[[layer]] <- Matrix::t(mat)
    }

	
    require(reticulate)
    anndata <- reticulate::import('anndata', convert = FALSE)

    adata <- anndata$AnnData(
        X = Matrix::t(X),
        obs = obs,
        obsm = obsm,
        obsp = obsp,
        var = var,
        varm = varm,
        varp = varp,
        layers = layers
    )

    if (!is.null(outFile))
        adata$write(outFile, compression = 'gzip')

    adata
}


ACE2AnnData.minimal <- function(ace, outFile = NULL, main_layer = 'logcounts') {
	
    X <- SummarizedExperiment::assay(ace, main_layer)

    obs <- preprocessDF(as.data.frame(SummarizedExperiment::colData(ace)))
    rownames(obs) = make.names(colnames(ace), unique = TRUE);
	
	
    var <- preprocessDF(as.data.frame(SummarizedExperiment::rowData(ace)))
	rownames(var) = make.names(rownames(ace), unique = TRUE)
	
	obsm = list(X_ACTIONet2D = ACTIONet.out$vis.out$coordinates, X_ACTIONet3D = ACTIONet.out$vis.out$coordinates_3D)
		
    require(reticulate)
    anndata <- reticulate::import('anndata', convert = FALSE)

    adata <- anndata$AnnData(
        X = Matrix::t(X),
        obs = obs,
        var = var,
        obsm = obsm
    )

    if (!is.null(outFile))
        adata$write(outFile, compression = 'gzip')

    adata
}
