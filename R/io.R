#' Imports data from a 10X experiment folder and constructs an `SingleCellExeriment` object
#'
#' @param input_path Folder containing input files.
#' @param mtx_file Count file in Matrix Market format (default='matrix.mtx.gz').
#' @param feature_metadata Table of the same size as number of rows in the count matrix (default='features.tsv.gz').
#' @param sample_metadata Table of the same size as number of columns in the count matrix (default='barcodes.tsv.gz').
#' @param sep Column-separator used in the row/column annotations files (default='\\t').
#' @param prefilter Whether to prefilter rows/columns of input counts matrix. Must specify filtering parameters to pass to filter.ace().
#'
#' @return `SingleCellExeriment` object
#'
#' @examples
#' ace = import.ace.from.10X.generic(input_path, prefilter=TRUE, min_feats_per_cell = 500)
#' @export
import.ace.from.10X.generic <- function(input_path, mtx_file = "matrix.mtx.gz", feature_metadata = "features.tsv.gz", 
    sample_metadata = "barcodes.tsv.gz", sep = "\t", use.names = T, prefilter = FALSE, 
    ...) {
    require(ACTIONet)
    require(S4Vectors)
    
    count.file = paste(input_path, mtx_file, sep = "/")
    if (!file.exists(count.file)) {
        err = sprintf("File %s not found. Consider changing `mtx_file` or `input_path` options.", 
            count.file)
        stop(err)
    }
    
    feature.file = paste(input_path, feature_metadata, sep = "/")
    if (!file.exists(feature.file)) {
        err = sprintf("File %s not found. Consider changing `mtx_file` or `input_path` options.", 
            feature.file)
        stop(err)
    }
    
    barcode.file = paste(input_path, sample_metadata, sep = "/")
    if (!file.exists(barcode.file)) {
        err = sprintf("File %s not found. Consider changing `mtx_file` or `input_path` options.", 
            barcode.file)
        stop(err)
    }
    
    message(sprintf("Reading counts ...\n"))
    counts.mat = Matrix::readMM(count.file)
    
    feature_table = read.table(feature.file, header = F, sep = sep, as.is = TRUE)
    if (nrow(feature_table) == (nrow(counts.mat) + 1)) {
        rowAnnot = S4Vectors::DataFrame(feature_table[-1, ])
        colnames(rowAnnot) = feature_table[1, ]
    } else {
        rowAnnot = S4Vectors::DataFrame(feature_table)
    }
    if (ncol(rowAnnot) == 1) {
        colnames(rowAnnot) = "Gene"
    } else if (ncol(rowAnnot) == 2) {
        colnames(rowAnnot) = c("ENSEMBL", "Gene")
    } else if (ncol(rowAnnot) == 3) {
        colnames(rowAnnot) = c("ENSEMBL", "Gene", "Feature")
    }
    
    sample_metadata = read.table(barcode.file, header = F, sep = sep, as.is = TRUE)
    if (ncol(sample_metadata) == (ncol(counts.mat) + 1)) {
        colAnnot = S4Vectors::DataFrame(sample_metadata[-1, ])
        colnames(colAnnot) = sample_metadata[1, ]
    } else {
        colAnnot = S4Vectors::DataFrame(sample_metadata)
    }
    
    # Feature-barcoding
    if (ncol(rowAnnot) > 2) {
        IDX = split(1:nrow(rowAnnot), rowAnnot[, 3])
        expression.counts.mat = counts.mat[IDX$`Gene Expression`, ]
        gene.table = rowAnnot[IDX$`Gene Expression`, 1:2]
        IDX = IDX[!grepl("Gene Expression", names(IDX))]
    } else {
        expression.counts.mat = counts.mat
        gene.table = rowAnnot
        IDX = list()
    }
    
    if (use.names == T && (ncol(rowAnnot) >= 2)) {
        rownames(expression.counts.mat) = gene.table[, 2]
    } else {
        rownames(expression.counts.mat) = gene.table[, 1]
    }
    
    colnames(expression.counts.mat) = colAnnot[, 1]
    
    
    if (ncol(sample_metadata) > 1) {
        ace <- ACTIONetExperiment(assays = list(counts = expression.counts.mat), 
            colData = colAnnot, rowData = rowAnnot)
    } else {
        ace <- ACTIONetExperiment(assays = list(counts = expression.counts.mat))
    }
    
    # Load additional barcoded features
    for (feature.name in names(IDX)) {
        feature.counts.mat = counts.mat[IDX[[feature.name]], ]
        row.annotations = rowAnnot[IDX[[feature.name]], ]
        rownames(feature.counts.mat) = rowAnnot[IDX[[feature.name]], 1]
        colnames(feature.counts.mat) = colAnnot[, 1]
        colMaps(ace)[[feature.name]] = Matrix::t(feature.counts.mat)
    }
    
    if (prefilter) {
        ace = filter.ace(ace, assay_name = "counts", return_fil_ace = TRUE, ...)
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
#' ace = import.ace.from.10X(input_path, prefilter=TRUE, min_feats_per_cell = 500)
#' @export
import.ace.from.10X <- function(input_path, version = 3, prefilter = FALSE, ...) {
    if (file.exists(paste(input_path, "genes.tsv", sep = "/"))) {
        version = 2
    }
    
    if ((2 <= version) & (version < 3)) {
        mtx_file = "matrix.mtx"
        feature_metadata = "genes.tsv"
        sample_metadata = "barcodes.tsv"
    } else if ((3 <= version) & (version < 4)) {
        mtx_file = "matrix.mtx.gz"
        feature_metadata = "features.tsv.gz"
        sample_metadata = "barcodes.tsv.gz"
    } else {
        message("Unknown version")
    }
    
    ace = import.ace.from.10X.generic(input_path = input_path, mtx_file = mtx_file, 
        feature_metadata = feature_metadata, sample_metadata = sample_metadata, prefilter = prefilter, 
        ...)
    
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
#' @export
import.ace.from.10X.h5 <- function(fname, version = 3, genome = NULL, use.names = TRUE, 
    prefilter = FALSE, ...) {
    if (!requireNamespace("hdf5r", quietly = TRUE)) {
        stop("Please install hdf5r to read HDF5 files")
    }
    if (!file.exists(fname)) {
        stop("File not found")
    }
    
    h5file <- hdf5r::H5File$new(filename = fname, mode = "r")
    if (is.null(genome)) 
        genome <- names(x = h5file)[[1]]
    
    if (h5file$attr_exists("PYTABLES_FORMAT_VERSION")) {
        version = 2
    }
    
    if ((version <= 3) & (version < 4)) {
        if (use.names) {
            feature_slot <- "features/name"
        } else {
            feature_slot <- "features/id"
        }
    } else if ((version <= 2) & (version < 3)) {
        if (use.names) {
            feature_slot <- "gene_names"
        } else {
            feature_slot <- "genes"
        }
    }
    counts.mat <- h5file[[paste0(genome, "/data")]]
    indices <- h5file[[paste0(genome, "/indices")]]
    indptr <- h5file[[paste0(genome, "/indptr")]]
    shp <- h5file[[paste0(genome, "/shape")]]
    features <- h5file[[paste0(genome, "/", feature_slot)]][]
    barcodes <- h5file[[paste0(genome, "/barcodes")]]
    sparse.mat <- sparseMatrix(i = indices[] + 1, p = indptr[], x = as.numeric(x = counts.mat[]), 
        dims = shp[], giveCsparse = FALSE)
    if (length(unique(features)) < length(features)) {
        features <- make.names(names = features, unique = T)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = "dgCMatrix")
    
    if (h5file$exists(name = paste0(genome, "/features"))) {
        types <- h5file[[paste0(genome, "/features/feature_type")]][]
        types.unique <- unique(x = types)
        if (length(x = types.unique) > 1) {
            message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
            mats <- sapply(X = types.unique, FUN = function(x) {
                return(sparse.mat[which(x = types == x), ])
            }, simplify = FALSE, USE.NAMES = TRUE)
        } else {
            mats = list(sparse.mat)
        }
    } else {
        mats = list(sparse.mat)
    }
    
    h5file$close_all()
    
    ace <- ACTIONetExperiment(assays = list(counts = mats[[1]]))
    
    # Load additional barcoded features
    for (feature.name in names(mats)[-1]) {
        colMaps(ace)[[feature.name]] = Matrix::t(mats[[feature.name]])
    }
    
    if (prefilter) {
        ace = filter.ace(ace, assay_name = "counts", return_fil_ace = TRUE, ...)
    }
    
    return(ace)
}

#' Constructs an `SingleCellExeriment` object from count matrix, gene names, and sample_metadata
#'
#' @param counts_mat Matrix of counts.mat
#' @param feature_metadata Vector or data frame  of annotations for features (rows of the count matrix)
#' @param sample_metadata Vector or data frame of annotations for samples (columns of the count matrix)
#' @param prefilter Whether to prefilter genes/cells based on the counts.mat
#' @param min.cell.frac.per.gene Minimum fraction of cells capturing a gene for it to be retained, if prefilter=TRUE (default=0.005)
#' @param min.genes.per.cell Minimum number of required captured genes per cell, if prefilter=TRUE (default=500)
#'
#' @return `SingleCellExeriment` object
#'
#' @examples
#' ace = import.ace.from.count.matrix(counts.mat.mat, gene_names, prefilter=TRUE)
#' @export
import.ace.from.counts <- function(counts_mat, feature_metadata = NULL, sample_metadata = NULL, 
    prefilter = FALSE, ...) {
    # if (!is.sparseMatrix(counts.mat)) { counts.mat = as(counts.mat, 'sparseMatrix')
    # } rownames(counts.mat) = gene.names
    
    if (is.null(feature_metadata)) {
        feature_metadata = .default_rowData(NROW(counts_mat))
    }
    
    if (is.null(sample_metadata)) {
        sample_metadata = .default_colData(NCOL(counts_mat))
    }
    
    ace <- ACTIONetExperiment(assays = list(counts = counts_mat), rowData = feature_metadata, 
        colData = sample_metadata)
    
    if (prefilter) {
        ace = filter.ace(ace, assay_name = "counts", return_fil_ace = TRUE, ...)
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
#' @export
import.ace.from.table <- function(fname, sep = "\t", prefilter = FALSE, ...) {
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
        ace = filter.ace(ace, assay_name = "counts", return_fil_ace = TRUE, ...)
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
#' @export
import.ace.from.Seurat <- function(Seurat.obj) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Please install Seurat to read Seurat objects")
    }
    
    Seurat.obj = Seurat::UpdateSeuratObject(Seurat.obj)
    sce = Seurat::as.SingleCellExperiment(Seurat.obj)
    ace <- as(sce, "ACTIONetExperiment")
    
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
#' @export
import.ace.from.loom <- function(fname, ...) {
    .check_and_load_package("sceasy")
    
    SE = SummarizedExperiment::makeSummarizedExperimentFromLoom(fname, ...)
    ace = as(SE, "ACTIONetExperiment")
    
    # ace <- as(sceasy:::readExchangeableLoom(fname), 'ACTIONetExperiment')
    
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
#' @export
import.ace.from.CDS <- function(monocle_cds) {
    .check_and_load_package("monocle3")
    
    counts.mat = exprs(monocle_cds)
    gene_annotations = fData(monocle_cds)
    sample_metadata = pData(monocle_cds)
    
    ace <- ACTIONetExperiment(assays = list(counts = counts.mat), colData = sample_metadata, 
        rowData = gene_annotations)
    
    # This one is NOT stable. SORRY!  sce = monocle::exportCDS(monocle_cds = CDS,
    # export_to = 'Scater', export_all = T)
    
    return(ace)
}


#' Translates rownames() of the `ace` object for mouse and human datasets.
#' A typical use-case is when input ids are in ENSEMBL format, but user is interested to work with gene symbols.
#'
#' @param ace Input `ace` object
#' @param from Source annotation (default='ENSEMBL')
#' @param to Target annotation (default='SYMBOL')
#' @param species Either 'mouse' or 'human' (default='human')
#'
#' @return `SingleCellExeriment` object with renamed rows
#'
#' @examples
#' ace = import.ace.from.CDS(monocle_cds)
#' @export
convert.ace.rownames <- function(ace, from = "ENSEMBL", to = "SYMBOL", species = "human") {
    if (species == "human") {
        library(org.Hs.eg.db)
        suppressWarnings(ids <- mapIds(org.Hs.eg.db, keys = row.names(ace), keytype = from, 
            column = to, multiVals = "first"))
        ids[is.na(ids)] = ""
        
        rownames(ace) = ids
    } else if (species == "mouse") {
        library(org.Mm.eg.db)
        suppressWarnings(ids <- mapIds(org.Mm.eg.db, keys = row.names(ace), keytype = from, 
            column = to, multiVals = "first"))
        ids[is.na(ids)] = ""
        
        rownames(ace) = ids
    }
    
    return(ace)
}


preprocessDF <- function(df, drop_single_values = TRUE) {
    if (ncol(df) > 0) {
        nn = colnames(df)
        for (n in nn) {
            x = df[, n]
            if (length(unique(x)) < 50) {
                x = factor(x, sort(unique(x)))
                df[, n] = x
            }
        }
    }
    if (ncol(df) == 0) 
        df[["name"]] <- rownames(df)
    if (drop_single_values) {
        k_singular <- sapply(df, function(x) length(unique(x)) == 1)
        if (sum(k_singular) > 0) 
            warning(paste("Dropping single category variables:"), paste(colnames(df)[k_singular], 
                collapse = ", "))
        df <- df[, !k_singular, drop = F]
        if (ncol(df) == 0) 
            df[["name"]] <- rownames(df)
    }
    return(df)
}

import.ace.from.legacy <- function(ACTIONet.out, sce, full.import = T, return.all = F) {
    ace = as(sce, "ACTIONetExperiment")
    
    if ("S_r" %in% names(colMaps(ace))) {
        colMaps(ace)[["ACTION"]] = colMaps(ace)[["S_r"]]
    }
    
    ACTION.out = ACTIONet.out$ACTION.out
    pruning.out = ACTIONet.out$reconstruct.out
    G = ACTIONet.out$build.out$ACTIONet
    
    colNets(ace)$ACTIONet = G
    vis.out = ACTIONet.out$vis.out
    
    colMaps(ace)$ACTIONet2D = vis.out$coordinates
    colMaps(ace)$ACTIONet3D = vis.out$coordinates_3D
    colMaps(ace)$denovo_color = vis.out$colors
    
    
    
    if (full.import == T) {
        colMaps(ace)[["H_stacked"]] = as(Matrix::t(ACTIONet.out$reconstruct.out$H_stacked), 
            "sparseMatrix")
        colMaps(ace)[["C_stacked"]] = as(ACTIONet.out$reconstruct.out$C_stacked, 
            "sparseMatrix")
    }
    
    unification.out = ACTIONet.out$unification.out
    colMaps(ace)[["H_unified"]] = as(Matrix::t(ACTIONet.out$unification.out$H.core), 
        "sparseMatrix")
    colMaps(ace)[["C_unified"]] = as(ACTIONet.out$unification.out$C.core, "sparseMatrix")
    
    ace$assigned_archetype = ACTIONet.out$unification.out$assignments.core
    ace$node_centrality = compute_archetype_core_centrality(colNets(ace)$ACTIONet, 
        sample_assignments = ace$assigned_archetype)
    
    # ace$node_centrality = compute_archetype_core_centrality(G,
    # ace$assigned_archetype)
    
    specificity.out = ACTIONet.out$unification.out$DE.core
    if (is.list(specificity.out)) {
        rowMaps(ace)[["unified_feature_profile"]] = specificity.out[["archetypes"]]
        rowMapTypes(ace)[["unified_feature_profile"]] = "internal"
        
        rowMaps(ace)[["unified_feature_specificity"]] = specificity.out[["upper_significance"]]
        rowMapTypes(ace)[["unified_feature_specificity"]] = "reduction"
    } else {
        rowMaps(ace)[["unified_feature_profile"]] = assays(specificity.out)[["profile"]]
        rowMapTypes(ace)[["unified_feature_profile"]] = "internal"
        rowMaps(ace)[["unified_feature_specificity"]] = assays(specificity.out)[["significance"]]
        rowMapTypes(ace)[["unified_feature_specificity"]] = "reduction"
    }
    
    # Prepare output
    if (return.all) {
        trace = list(ACTION.out = ACTION.out, pruning.out = pruning.out, vis.out = vis.out, 
            unification.out = unification.out)
        trace$log = list(genes = rownames(ace), cells = colnames(ace), time = Sys.time())
        
        out = list(ace = ace, trace = trace)
        
        return(out)
    } else {
        return(ace)
    }
}

export.sparse.format <- function(ace, path) {
    counts.mat = counts(ace)
    cell_metadata = colData(ace)
    features = rownames(ace)
    
    write.table(features, file.path(path, "features.txt"), row.names = F, col.names = F, 
        quote = F)
    write.table(cell_metadata, file.path(path, "metadata.tsv"), row.names = F, col.names = T, 
        quote = F, sep = "\t")
    Matrix::writeMM(counts.mat, file.path(path, "counts.mtx"))
}
