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
	X = as(X, 'dgCMatrix')
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

	Dims = attr$shape
	if(attr[["encoding-type"]] == "csc_matrix") {
		csc_sort_indices_inplace(indptr, indices, data)		
		Xt = sparseMatrix(i = indices+1, p = indptr+1, x = data, dims = Dims)
		X = Matrix::t(Xt)
	} else if(attr[["encoding-type"]] == "csr_matrix") {
		csr_sort_indices_inplace(indptr, indices, data)
		Xt = new("dgRMatrix", j = indices+1, p = indptr+1, x = data, Dim = Dims)
		X = Matrix::t(Xt)
	}
	
	return(X)
}

ACE2AnnData <- function(ace, fname = "ACTIONet.h5ad", main.assay = "logcounts", minimal.export = F, compression.level = 0) {
	if (!requireNamespace('hdf5r', quietly = TRUE)) {
		stop("Please install hdf5r to wrote HDF5 files")
	}
	
  	if(file.exists(fname)) {
		file.remove(fname)
	}	
	
	h5file = H5File$new(fname, mode = "w")

	## Write X (logcounts() in ace, in either sparse or dense format)
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



	## Write obs (colData() in ace)
	obs.DF = as.data.frame(colData(ace))
	write.HD5DF(h5file, gname = "obs", obs.DF, compression.level = compression.level)

	## Write var (matching rowData() in ace)
	var.DF = as.data.frame(rowData(ace))
	write.HD5DF(h5file, "var", var.DF, compression.level = compression.level)

	## Write subset of obsm related to the cell embeddings (Dim=2 or 3)
	obsm = h5file$create_group("obsm")
	CF = colMaps(ace)
	embeddings.idx = which(sapply(ace@colMapsAnnot, function(x) x$type == "embedding"))
	if(length(embeddings.idx) > 0) {
		subCF = CF[embeddings.idx]
		names(subCF) = paste("X", names(subCF), sep = "_")
		subCF = lapply(subCF, function(x) as.matrix(x))
		for(i in 1:length(subCF)) {
			obsm$create_dataset(names(subCF)[[i]], subCF[[i]], gzip_level = compression.level, dtype = h5types$H5T_IEEE_F32LE)
		}
	}
	if(!minimal.export) {
		# Export additional "obsm" matrices. Anything that doesn't start with "X_" in obsm will be map to colMaps() upon reading.
		nonembeddings.idx = setdiff(1:length(CF), embeddings.idx)
		if((length(nonembeddings.idx) > 0)) {
			subCF = CF[nonembeddings.idx]
			subCF = lapply(subCF, function(x) as.matrix(x))
			for(i in 1:length(subCF)) {
				obsm$create_dataset(names(subCF)[[i]], subCF[[i]], gzip_level = compression.level, dtype = h5types$H5T_IEEE_F32LE)
			}
		}
		
		# Export "varm"-associated matrices, i.e. rowMaps(): variables in AnnData ~ rows in SCE ~ genes => such as DE matrices
		RF = rowMaps(ace)
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
	if (!requireNamespace('hdf5r', quietly = TRUE)) {
		stop("Please install hdf5r to read HDF5 files")
	}
	
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

	ace = ACTIONetExperiment(assays = assays, rowData = var.DF, colData = obs.DF)

	if("obsm" %in% objs) {
		obsm = h5file[["obsm"]]
		for(mn in names(obsm)) {
			Xr = obsm[[mn]]$read()
			if(sum(grepl(pattern = "^X_", mn))) { # # Keep the smallest factors as "embeddings"
				nn = stringr::str_sub(mn, start = 3)
				colMaps(ace)[[nn]] = Xr
				ace@colMapsAnnot[[nn]] = list(type = "embedding")
			} else if(nrow(Xr) <= 100 & (sum(grepl(pattern = "^C_", mn)|grepl(pattern = "^H_", mn))==0)) { # Keep small factors as "reductions"
				colMaps(ace)[[mn]] = Xr
				ace@colMapsAnnot[[mn]] = list(type = "reduction")								
			} else {
				colMaps(ace)[[mn]] = Xr
				ace@colMapsAnnot[[mn]] = list(type = "generic")								
			}
		}
	}

	if("varm" %in% objs) {
		varm = h5file[["varm"]]
		for(mn in names(varm)) {
			Xr = Matrix::t(varm[[mn]]$read())
			rowMaps(ace)[[mn]] = Xr
			ace@rowMapsAnnot[[mn]] = list(type = "generic")
		}
	}

	if("obsp" %in% objs) {
		obsp = h5file[["obsp"]]
		for(pn in names(obsp)) {
			Net = read.HD5SpMat(obsp, pn)
			colNets(ace)[[pn]] = Net
		}
	}

	if("varp" %in% objs) {
		varp = h5file[["varp"]]
		for(pn in names(varp)) {
			Net = read.HD5SpMat(varp, pn)
			rowNets(ace)[[pn]] = Net
		}
	}

	h5file$close_all()

	return(ace)
}
