h5addAttr.str <- function(h5group, attr.name, attr.val) {
    dtype = H5T_STRING$new(type="c", size=Inf)
    dtype = dtype$set_cset(cset = "UTF-8")

    space = H5S$new(type="scalar")
	h5group$create_attr(attr_name = attr.name, dtype = dtype, space = space)
	attr = h5group$attr_open_by_name(attr_name = attr.name, ".")
	attr$write(attr.val)
}

h5addAttr.str_array <- function(h5group, attr.name, attr.val) {
    #dtype <- guess_dtype(x=attr.val, scalar=F, string_len=Inf)
    dtype = H5T_STRING$new(type="c", size=Inf)
    dtype = dtype$set_cset(cset = "UTF-8")

    space = H5S$new(type="simple", dims = length(attr.val), maxdims = length(attr.val))
	h5group$create_attr(attr_name = attr.name, dtype = dtype, space = space)
	attr = h5group$attr_open_by_name(attr_name = attr.name, ".")
	attr$write(attr.val)
}

write.HD5DF <- function(h5file, gname, DF, compression.level = 0) {

	string.dtype = H5T_STRING$new(type="c", size=Inf)    
	string.dtype = string.dtype$set_cset(cset = "UTF-8")

	N = nrow(DF)

	cat.vars = which(sapply(1:ncol(DF), function(i) length(unique(DF[, i])) < 256 ))
	noncat.num.vars = which(sapply(1:ncol(DF), function(i) {
		if(is.numeric(DF[, i])) 
			return(sum(round(DF[, i]) != DF[, i]) != 0)
		else 
			return(FALSE)
	}))
	cat.vars = setdiff(cat.vars, noncat.num.vars)
	cn = colnames(DF)[c(cat.vars, noncat.num.vars)]

	h5group = h5file$create_group(gname)

	h5addAttr.str(h5group, "_index", "index")
	h5addAttr.str(h5group, "encoding-version", "0.1.0")
	h5addAttr.str(h5group, "encoding-type", "dataframe")
	

	if(length(cn) == 0) {
		dtype = H5T_STRING$new(type="c", size=Inf)
		dtype = dtype$set_cset(cset = "UTF-8")
		space = H5S$new(type="simple", dims = 0, maxdims = 10)
		
		h5group$create_attr(attr_name = "column-order", dtype = dtype, space = space)
		
		# attr = h5group$attr_open_by_name(attr_name = "column-order", ".")
		# attr$write()			
	} else {
		h5addAttr.str_array(h5group, "column-order", cn)			
	}

	

		
	
	if(length(cat.vars) > 0) {
		cat = h5group$create_group("__categories")

		for(i in 1:length(cat.vars)) {
			x =  DF[, cat.vars[i]]
			if(class(x) == "factor") {
				l = as.character(levels(x))
				v = as.numeric(x) - 1
			} else {
				x = as.character(x)
				l = sort(unique(x))
				v = match(x, l) - 1
			}

			dtype = H5T_STRING$new(type="c", size=Inf)
			dtype = dtype$set_cset(cset = "UTF-8")
			l.enum = cat$create_dataset(colnames(DF)[cat.vars[i]], l, gzip_level = compression.level, dtype = dtype)


		    dtype = H5T_ENUM$new(labels = c("FALSE", "TRUE"), values = 0:1)
		    space = H5S$new(type="scalar")
			res = l.enum$create_attr(attr_name = "ordered", dtype = dtype, space = space)

			attr = l.enum$attr_open_by_name(attr_name = "ordered", ".")
			attr$write(0)

			l.vec = h5group$create_dataset(colnames(DF)[cat.vars[i]], as.integer(v), gzip_level = compression.level, dtype = h5types$H5T_NATIVE_INT8)

			ref = cat$create_reference(name = colnames(DF)[cat.vars[i]])

		    dtype = guess_dtype(ref)
		    space = H5S$new(type="scalar")
			res = l.vec$create_attr(attr_name = "categories", dtype = dtype, space = space)
			attr = l.vec$attr_open_by_name(attr_name = "categories", ".")
			attr$write(ref)

		}
	}
	index = rownames(DF)
	if(length(unique(index)) < length(index)) {
		index = make.names(index, unique = TRUE)
	}


	h5group$create_dataset("index", index, gzip_level = compression.level, dtype = string.dtype)
	if(length(noncat.num.vars) > 0) {
		for(i in noncat.num.vars) {
			x = DF[, i]
			nn = colnames(DF)[i]
			h5group$create_dataset(nn, as.single(x), gzip_level = compression.level, dtype = h5types$H5T_IEEE_F32LE)
		}	
	}
}

write.HD5SpMat <- function(h5file, gname, X, compression.level = compression.level) {
	X = Matrix::t(as(X, 'dgCMatrix'))
	Xgroup = h5file$create_group(gname)


	Xgroup$create_dataset("indices", X@i, gzip_level = compression.level, dtype = h5types$H5T_NATIVE_INT32)
	Xgroup$create_dataset("indptr", X@p, gzip_level = compression.level, dtype = h5types$H5T_NATIVE_INT32)
	Xgroup$create_dataset("data", as.single(X@x), gzip_level = compression.level, dtype = h5types$H5T_IEEE_F32LE)

	h5addAttr.str(Xgroup, "encoding-type", "csc_matrix")
	h5addAttr.str(Xgroup, "encoding-version", "0.1.0")
	h5attr(Xgroup, "shape") = dim(X)
}

read.HD5DF <- function(h5file, gname, compression.level = 0) {

	h5group = h5file[[gname]]

	if(! (h5group$attr_open_by_name("encoding-type", ".")$read() == "dataframe") ) {
		R.utils::printf("%s is not a dataframe. Abort.\n", gname)
		return()
	}
	rn = h5group[[h5group$attr_open_by_name("_index", ".")$read()]]$read()
	
	cn = h5group$attr_open_by_name("column-order", ".")
	if(cn$get_storage_size() == 0) {
		DF = DataFrame(row.names = rn)	
		return(DF)		
	}

	column.names = cn$read()
	vars = vector("list", length(column.names))
	names(vars) = column.names
	for(vn in names(vars)) {
		vars[[vn]] = h5group[[vn]]$read()
	}

	if("__categories" %in% names(h5group)) {
		cat = h5group[["__categories"]]
		for(nn in names(cat)) {
			l = cat[[nn]]$read()
			vars[[nn]] = factor(l[vars[[nn]]+1], l)
		}
	}

	DF = DataFrame(vars)
	rownames(DF) = rn

	return(DF)
}

read.HD5SpMat <- function(h5file, gname, compression.level = compression.level) {

	h5group = h5file[[gname]]
	attr = h5attributes(h5group)
	if(! ( ("encoding-type" %in% names(attr)) & (attr[["encoding-type"]] %in% c("csc_matrix", "csr_matrix"))) ) {
		R.utils::printf("%s is not a sparse matrix. Abort.\n", gname)
		return()
	}

	data = h5group[["data"]]$read()
	indices = h5group[["indices"]]$read()
	indptr = h5group[["indptr"]]$read()

	if(attr[["encoding-type"]] == "csc_matrix") {
		X = Matrix::t(new("dgCMatrix", i = indices, p = indptr, x = data, Dim = attr$shape))
	} else if(attr[["encoding-type"]] == "csr_matrix") {
		X = Matrix::t(new("dgRMatrix", j = indices, p = indptr, x = data, Dim = attr$shape))
	}

	return(X)
}

ACE2AnnData <- function(ace, fname = "ACTIONet.h5ad", main.assay = "logcounts", minimal.export = F, compression.level = 0) {
	if(file.exists(fname)) {
		file.remove(fname)
	}	
	
	h5file = H5File$new(fname, mode = "w")

	## Write X (logcounts() in ACE, in either sparse or dense format)
	X = assays(ace)[[main.assay]]
	if(is.sparseMatrix(X)) {
		write.HD5SpMat(h5file, gname = "X", X, compression.level = compression.level)
	} else {
		 h5file$create_dataset("X", X, gzip_level = compression.level, dtype = h5types$H5T_IEEE_F32LE)
	}

	remaining.assays = setdiff(names(assays(ace)), main.assay)
	if( (minimal.export == F) & (0 < length(remaining.assays)) ) {
		layers = h5file$create_group("layers")

		for(an in remaining.assays) {
			Xr = as(assays(ace)[[an]], 'dgCMatrix')
			write.HD5SpMat(layers, gname = an, Xr, compression.level = compression.level)
		}
	}



	## Write obs (colData() in ACE)
	obs.DF = as.data.frame(colData(ace))
	write.HD5DF(h5file, gname = "obs", obs.DF, compression.level = compression.level)

	## Write var (matching rowData() in ACE)
	var.DF = as.data.frame(rowData(ace))
	write.HD5DF(h5file, "var", var.DF, compression.level = compression.level)

	## Write subset of obsm related to the cell embeddings (reducedDims() with 2 or 3 columns)
	obsm = h5file$create_group("obsm")
	RD = reducedDims(ace)
	embeddings.idx = which(sapply(RD, ncol) %in% c(2, 3))
	if(length(embeddings.idx) > 0) {
		subRD = RD[embeddings.idx]
		names(subRD) = paste("X", names(subRD), sep = "_")
		subRD = lapply(subRD, function(x) as.matrix(x))
		for(i in 1:length(subRD)) {
			obsm$create_dataset(names(subRD)[[i]], Matrix::t(subRD[[i]]), gzip_level = compression.level, dtype = h5types$H5T_IEEE_F32LE)
		}
	}
	if(!minimal.export) {
		# Export additional "obsm" matrices. Anything that doesn't start with "X_" in obsm will be map to colFactors() upon reading.
		nonembeddings.idx = setdiff(1:length(reducedDims(ace)), embeddings.idx)
		if((length(nonembeddings.idx) > 0)) {
			subRD = RD[nonembeddings.idx]
			subRD = lapply(subRD, function(x) as.matrix(x))
			for(i in 1:length(subRD)) {
				obsm$create_dataset(names(subRD)[[i]], Matrix::t(subRD[[i]]), gzip_level = compression.level, dtype = h5types$H5T_IEEE_F32LE)
			}
		}


		# Export additional "obsm"-associated matrices, i.e. colFactors(): obs in AnnData ~ columns in SCE ~ cells => AA results
		CF = colFactors(ace)
		if((length(CF) > 0)) {
			CF = lapply(CF, function(x) Matrix::t(as.matrix(x)))
			for(i in 1:length(CF)) {
				obsm$create_dataset(names(CF)[[i]], Matrix::t(CF[[i]]), gzip_level = compression.level, dtype = h5types$H5T_IEEE_F32LE)
			}
		}

		# Export "varm"-associated matrices, i.e. rowFactors(): variables in AnnData ~ rows in SCE ~ genes => such as DE matrices
		RF = rowFactors(ace)
		if((length(RF) > 0)) {
			varm = h5file$create_group("varm")
			RF = lapply(RF, function(x) as.matrix(x))
			for(i in 1:length(RF)) {
				varm$create_dataset(names(RF)[[i]], Matrix::t(RF[[i]]), gzip_level = compression.level, dtype = h5types$H5T_IEEE_F32LE)
			}
		}


		# Export "obsp"-associated matrices, i.e. colNets(): obs in AnnData ~ cols in SCE ~ cells => cell-cell networks (such as ACTIONet)
		CN = colNets(ace)
		if((length(CN) > 0)) {
			obsp = h5file$create_group("obsp")
			CN = lapply(CN, function(x) as(x, "dgCMatrix"))

			for(i in 1:length(CN)) {
				write.HD5SpMat(obsp, gname = names(CN)[[i]], CN[[i]], compression.level = compression.level)
			}
		}


		# Export "varp"-associated matrices, i.e. rowNets(): var in AnnData ~ rows in SCE ~ genes => gene-gene networks (such as SCINET)
		RN = rowNets(ace)
		if((length(RN) > 0)) {
			varp = h5file$create_group("varp")
			RN = lapply(RN, function(x) as(x, "dgCMatrix"))

			for(i in 1:length(RN)) {
				write.HD5SpMat(varp, gname = names(RN)[[i]], RN[[i]], compression.level = compression.level)
			}
		}
	}

	h5file$close_all()

}

AnnData2ACE <- function(fname = "ACTIONet.h5ad", main.assay = "logcounts", minimal.export = T, compression.level = 0) {

	h5file = H5File$new(fname, mode = "r")

	objs = names(h5file)
	missing.elemenets = setdiff(list("X", "obs", "var"), objs)
	if(0 < length(missing.elemenets) ) {
		R.utils::printf("[%s] missing from the h5ad file. Abort.\n", paste(missing.elemenets, sep = ','))
		return()
	}


	X.attr = h5attributes(h5file[["X"]])
	if(length(X.attr) == 0) { # Full matrix
		X = h5file[["X"]]$read()
	} else {
		X = read.HD5SpMat(h5file = h5file, gname = "X", compression.level = compression.level)
	}
	assays = list(X)
	names(assays) = main.assay

	if("layers" %in% objs) {
		layers = h5file[["layers"]]
		additional_assays = vector("list", length(names(layers)))
		names(additional_assays) = names(layers)

		for(an in names(layers)) {
			attr = h5attributes(layers[[an]])
			if(length(attr) == 0) { # Dense matrix
				additional_assays[[an]] = layers[[an]]$read()
			} else {
				additional_assays[[an]] = read.HD5SpMat(h5file = layers, gname = an, compression.level = compression.level)
			}
		}
		assays = c(assays, additional_assays)
	}


	if("obs" %in% objs) {
		obs.DF = read.HD5DF(h5file = h5file, gname = "obs", compression.level = compression.level)			
	}

	if("var" %in% objs) {
		var.DF = read.HD5DF(h5file = h5file, gname = "var", compression.level = compression.level)
	}
	
	assays = lapply(assays, function(X) {
		rownames(X) = rownames(var.DF)
		colnames(X) = rownames(obs.DF)
		return(X)
	})

	ACE = ACTIONetExperiment(assays = assays, rowData = var.DF, colData = obs.DF)

	if("obsm" %in% objs) {
		obsm = h5file[["obsm"]]
		for(mn in names(obsm)) {
			Xr = obsm[[mn]]$read()
			if(sum(grepl(pattern = "^X_", mn))) { # It is an embedding/dimension reduction
				reducedDims(ACE)[[stringr::str_sub(mn, start = 3)]] = Matrix::t(Xr)
			} else if(nrow(Xr) <= 100 & (sum(grepl(pattern = "^C_", mn)|grepl(pattern = "^H_", mn))==0)) { # Keep small factors as reducedDims() -- i.e., "ACTION" or "PCA" reductions
				reducedDims(ACE)[[mn]] = Matrix::t(Xr)
			} else {
				colFactors(ACE)[[mn]] = Xr
			}
		}
	}

	if("varm" %in% objs) {
		varm = h5file[["varm"]]
		for(mn in names(varm)) {
			Xr = Matrix::t(varm[[mn]]$read())
			rowFactors(ACE)[[mn]] = Xr
		}
	}

	if("obsp" %in% objs) {
		obsp = h5file[["obsp"]]
		for(pn in names(obsp)) {
			Net = read.HD5SpMat(obsp, pn)
			colNets(ACE)[[pn]] = Net
		}
	}

	if("varp" %in% objs) {
		varp = h5file[["varp"]]
		for(pn in names(varp)) {
			Net = read.HD5SpMat(varp, pn)
			rowNets(ACE)[[pn]] = Net
		}
	}

	h5file$close_all()

	return(ACE)
}


AnnData2ACE.python <- function(inFile) {
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

ACE2AnnData.python <- function(ace, outFile = NULL, main_layer = 'logcounts', transfer_layers = c()) {

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


ACE2AnnData.minimal.python <- function(ace, outFile = NULL, main_layer = 'logcounts') {

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
