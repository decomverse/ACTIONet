#' @import hdf5r
h5addAttr.str <- function(h5group, attr.name, attr.val) {

    dtype = H5T_STRING$new(type = "c", size = Inf)
    dtype = dtype$set_cset(cset = "UTF-8")

    space = H5S$new(type = "scalar")
    h5group$create_attr(attr_name = attr.name, dtype = dtype, space = space)
    attr = h5group$attr_open_by_name(attr_name = attr.name, ".")
    attr$write(attr.val)
}

#' @import hdf5r
h5addAttr.str_array <- function(h5group, attr.name, attr.val) {

    # dtype <- guess_dtype(x=attr.val, scalar=F, string_len=Inf)
    dtype = H5T_STRING$new(type = "c", size = Inf)
    dtype = dtype$set_cset(cset = "UTF-8")

    space = H5S$new(type = "simple", dims = length(attr.val), maxdims = length(attr.val))
    h5group$create_attr(attr_name = attr.name, dtype = dtype, space = space)
    attr = h5group$attr_open_by_name(attr_name = attr.name, ".")
    attr$write(attr.val)
}

#' @import hdf5r
write.HD5DF <- function(
  h5file,
  gname,
  DF,
  compression_level = 0
) {

    string.dtype = H5T_STRING$new(type = "c", size = Inf)
    string.dtype = string.dtype$set_cset(cset = "UTF-8")

    DF = as.data.frame(DF)

    N = NROW(DF)

    h5group = h5file$create_group(gname)

    h5addAttr.str(h5group, "_index", "index")
    h5addAttr.str(h5group, "encoding-version", "0.1.0")
    h5addAttr.str(h5group, "encoding-type", "dataframe")

    if (0 < NCOL(DF)) {
        cat.vars = which(sapply(1:NCOL(DF), function(i) length(unique(DF[, i])) <
            256))

        noncat.vars = setdiff(1:NCOL(DF), cat.vars)
        if (length(noncat.vars) > 0) {
            noncat.num.vars = noncat.vars[sapply(noncat.vars, function(i) {
                x = as.numeric(DF[, i])
                return(sum(!is.na(x)) > 0)
            })]
        } else {
            noncat.num.vars = noncat.vars
        }

        noncat.nonnum.vars = setdiff(noncat.vars, noncat.num.vars)

        if (length(cat.vars) > 0) {
            cat.vars = setdiff(cat.vars, noncat.num.vars)
        }

        # cn = colnames(DF)[c(cat.vars, noncat.num.vars)]
        cn = colnames(DF)
        catDF = DF[, cat.vars, drop = F]
        catDF = apply(catDF, 2, as.character)
        catDF[is.na(catDF)] = "NA"

        numDF = DF[, noncat.num.vars, drop = F]
        numDF = apply(numDF, 2, as.numeric)
        numDF[is.na(numDF)] = NA

        nonNumDF = DF[, noncat.nonnum.vars, drop = F]
        nonNumDF = apply(nonNumDF, 2, as.character)
        nonNumDF[is.na(nonNumDF)] = NA

        if (length(cn) == 0) {
            dtype = H5T_STRING$new(type = "c", size = Inf)
            dtype = dtype$set_cset(cset = "UTF-8")
            space = H5S$new(type = "simple", dims = 0, maxdims = 10)

            h5group$create_attr(attr_name = "column-order", dtype = dtype, space = space)

            # attr = h5group$attr_open_by_name(attr_name = 'column-order', '.') attr$write()
        } else {
            h5addAttr.str_array(h5group, "column-order", cn)
        }

        if (length(cat.vars) > 0) {
            cat = h5group$create_group("__categories")

            for (i in 1:length(cat.vars)) {
                x = catDF[, i]  #DF[, cat.vars[i]]
                if (class(x) == "factor") {
                  l = as.character(levels(x))
                  v = as.numeric(x) - 1
                } else {
                  x = as.character(x)
                  l = sort(unique(x))
                  v = match(x, l) - 1
                }

                dtype = H5T_STRING$new(type = "c", size = Inf)
                dtype = dtype$set_cset(cset = "UTF-8")
                l.enum = cat$create_dataset(colnames(DF)[cat.vars[i]], l, gzip_level = compression_level,
                  dtype = dtype)


                dtype = H5T_ENUM$new(labels = c("FALSE", "TRUE"), values = 0:1)
                space = H5S$new(type = "scalar")
                res = l.enum$create_attr(attr_name = "ordered", dtype = dtype, space = space)

                attr = l.enum$attr_open_by_name(attr_name = "ordered", ".")
                attr$write(0)

                l.vec = h5group$create_dataset(colnames(DF)[cat.vars[i]], as.integer(v),
                  gzip_level = compression_level, dtype = h5types$H5T_NATIVE_INT8)

                ref = cat$create_reference(name = colnames(DF)[cat.vars[i]])

                dtype = guess_dtype(ref)
                space = H5S$new(type = "scalar")
                res = l.vec$create_attr(attr_name = "categories", dtype = dtype,
                  space = space)
                attr = l.vec$attr_open_by_name(attr_name = "categories", ".")
                attr$write(ref)

            }
        }
        if (length(noncat.num.vars) > 0) {
            for (i in 1:NCOL(numDF)) {
                x = numDF[, i]
                nn = colnames(numDF)[i]
                h5group$create_dataset(nn, as.single(x), gzip_level = compression_level,
                  dtype = h5types$H5T_IEEE_F32LE)
            }
        }

        if (length(noncat.nonnum.vars) > 0) {
            for (i in 1:NCOL(nonNumDF)) {
                x = nonNumDF[, i]
                nn = colnames(nonNumDF)[i]
                dtype = H5T_STRING$new(type = "c", size = Inf)
                dtype = dtype$set_cset(cset = "UTF-8")
                h5group$create_dataset(nn, x, gzip_level = compression_level,
                  dtype = string.dtype)
            }
        }
    } else {
        dtype = H5T_STRING$new(type = "c", size = Inf)
        dtype = dtype$set_cset(cset = "UTF-8")
        space = H5S$new(type = "simple", dims = 0, maxdims = 10)

        h5group$create_attr(attr_name = "column-order", dtype = dtype, space = space)
    }

    index = rownames(DF)
    if (length(unique(index)) < length(index)) {
        index = make.names(index, unique = TRUE)
    }
    h5group$create_dataset("index", index, gzip_level = compression_level, dtype = string.dtype)
}

#' @import hdf5r
write.HD5SpMat <- function(
  h5file,
  gname,
  X,
  compression_level = 0
) {

    X = Matrix::t(as(X, "dgCMatrix"))
    Xgroup = h5file$create_group(gname)


    Xgroup$create_dataset("indices", X@i, gzip_level = compression_level, dtype = h5types$H5T_NATIVE_INT32)
    Xgroup$create_dataset("indptr", X@p, gzip_level = compression_level, dtype = h5types$H5T_NATIVE_INT32)
    Xgroup$create_dataset("data", as.single(X@x), gzip_level = compression_level,
        dtype = h5types$H5T_IEEE_F32LE)

    h5addAttr.str(Xgroup, "encoding-type", "csc_matrix")
    h5addAttr.str(Xgroup, "encoding-version", "0.1.0")
    h5attr(Xgroup, "shape") = dim(X)
}

#' @import hdf5r
write.HD5List <- function(
  h5file,
  gname,
  obj_list,
  depth = 1,
  max_depth = 5,
  compression_level = 0
) {

    h5group = h5file$create_group(gname)

    obj_list = as.list(obj_list)

    for (nn in names(obj_list)) {
        obj = obj_list[[nn]]
        if ((sum(sapply(c("list", "SimpleList"), function(x) return(length(which(is(obj) ==
            x)) != 0))) != 0) & (depth < max_depth)) {
            write.HD5List(h5group, nn, obj, depth = depth + 1, max_depth = max_depth,
                compression_level = compression_level)
        } else if (sum(sapply(c("data.frame", "DataFrame", "DFrame"), function(x) return(length(which(is(obj) ==
            x)) != 0))) != 0) {
            write.HD5DF(h5group, nn, obj, compression_level = compression_level)
        } else if (is.sparseMatrix(obj)) {
            write.HD5SpMat(h5group, nn, obj, compression_level = compression_level)
        } else if (is.matrix(obj) | is.numeric(obj)) {
            h5group$create_dataset(nn, obj, gzip_level = compression_level, dtype = h5types$H5T_IEEE_F32LE)
        } else {
            h5group[[nn]] = obj
        }
    }
}

#' @import hdf5r
read.HD5DF <- function(
  h5file,
  gname,
  compression_level = 0
) {

    h5group = h5file[[gname]]

    if (!(h5group$attr_open_by_name("encoding-type", ".")$read() == "dataframe")) {
        err = sprintf("%s is not a dataframe. Abort.\n", gname)
        stop(err)
    }

    rn = h5group[[h5group$attr_open_by_name("_index", ".")$read()]]$read()

    if (h5file$attr_exists_by_name(obj_name = gname, attr_name = "column-order")) {
        cn = h5group$attr_open_by_name("column-order", ".")
        if (cn$get_storage_size() == 0) {
            DF = DataFrame(row.names = rn)
            return(DF)
        }

        column.names = cn$read()
        vars = vector("list", length(column.names))
        names(vars) = column.names
        for (vn in names(vars)) {
            vars[[vn]] = h5group[[vn]]$read()
        }

        if ("__categories" %in% names(h5group)) {
            cat = h5group[["__categories"]]
            for (nn in names(cat)) {
                l = cat[[nn]]$read()
                l = setdiff(l, "NA")
                if (length(l) < 2)
                  next

                vars[[nn]] = factor(l[vars[[nn]] + 1], l)
            }
        }
        DF = DataFrame(vars)
        rownames(DF) = rn
    } else {
        DF = DataFrame(row.names = rn)
    }

    return(DF)
}

#' @import hdf5r
read.HD5SpMat <- function(
  h5file,
  gname,
  compression_level = 0
) {

    h5group = h5file[[gname]]
    attr = h5attributes(h5group)
    if (!(("encoding-type" %in% names(attr)) & (attr[["encoding-type"]] %in% c("csc_matrix",
        "csr_matrix")))) {
        err = sprintf("%s is not a sparse matrix. Abort.\n", gname)
        stop(err)
    }

    data = h5group[["data"]]$read()
    indices = h5group[["indices"]]$read()
    indptr = h5group[["indptr"]]$read()

    Dims = attr$shape
    if (attr[["encoding-type"]] == "csc_matrix") {
        csc_sort_indices_inplace(indptr, indices, data)
        Xt = Matrix::sparseMatrix(i = indices + 1, p = indptr, x = data, dims = Dims)
        X = Matrix::t(Xt)
    } else if (attr[["encoding-type"]] == "csr_matrix") {
        # csr_sort_indices_inplace(indptr, indices, data) Xt = new('dgRMatrix', j =
        # indices, p = indptr, x = data, Dim = Dims) X = Matrix::t(Xt)
        csc_sort_indices_inplace(indptr, indices, data)
        Xt = Matrix::sparseMatrix(j = indices + 1, p = indptr, x = data, dims = Dims)
        X = Matrix::t(Xt)
    }

    return(X)
}

#' @import hdf5r
read.HD5List <- function(
  h5file,
  gname,
  depth = 1,
  max_depth = 5,
  compression_level = 0
) {

    h5group = h5file[[gname]]


    obj_names = names(h5group)
    L.out = vector("list", length(obj_names))
    names(L.out) = obj_names

    if (length(obj_names) > 0) {
        for (nn in obj_names) {
            attr = h5attributes(h5group[[nn]])
            if (length(attr) > 0 & ("encoding-type" %in% names(attr))) {
                if ((attr[["encoding-type"]] == "csc_matrix") | (attr[["encoding-type"]] ==
                  "csr_matrix")) {
                  obj = read.HD5SpMat(h5group, nn, compression_level = compression_level)
                } else if ((attr[["encoding-type"]] == "dataframe")) {
                  obj = read.HD5DF(h5group, nn, compression_level = compression_level)
                } else {
                  warning(sprintf("Unknown encoding %s", attr[["encoding-type"]]))
                  next
                }
            } else if (h5group[[nn]]$get_obj_type() == 2 & (depth < max_depth)) {
                obj = read.HD5List(h5group, nn, compression_level = compression_level,
                  depth = depth + 1)
            } else {
                obj = h5group[[nn]]$read()
            }
            L.out[[nn]] = obj
        }
        filter.mask = sapply(L.out, function(x) is.null(x))
        if (sum(filter.mask) > 0) {
            L.out = L.out[!filter.mask]
        }
    }

    return(L.out)
}

#' @import hdf5r
#' @export
ACE2AnnData <- function(
  ace,
  file,
  main_assay = "logcounts",
  full.export = TRUE,
  compression_level = 0
) {

    .check_and_load_package("hdf5r")

    # Ensure it can be case as an ACE object
    ace = as(ace, "ACTIONetExperiment")

    if (file.exists(file)) {
        file.remove(file)
    }

    if (is.null(colnames(ace))) {
        colnames(ace) = .default_colnames(NCOL(ace))
    }
    if (is.null(rownames(ace))) {
        rownames(ace) = .default_rownames(NROW(ace))
    }

    # Error is no assay specified.
    # assay_opts = c(main_assay, raw_assay)
    # if (!(main_assay %in% names(SummarizedExperiment::assays(ace)))) {
    #     err = sprintf("Invalid valid assay selection.\n")
    #     stop(err)
    # }

    # Make row/column-names unique
    colnames(ace) = ucn = .make_chars_unique(colnames(ace))
    rownames(ace) = urn = .make_chars_unique(rownames(ace))

    for (nn in names(SummarizedExperiment::assays(ace))) {
        dimnames(SummarizedExperiment::assays(ace)[[nn]]) = list(urn, ucn)
    }

    h5file = H5File$new(file, mode = "w")

    ## Write X (assays in ace, in either sparse or dense format)
    if(is.null(main_assay)){
      main_mat = Matrix::sparseMatrix(i = c(), j = c(), dims = dim(ace))
      write.HD5SpMat(h5file, gname = "X", main_mat, compression_level = compression_level)
    } else{
      if (!(main_assay %in% names(SummarizedExperiment::assays(ace)))) {
        err = sprintf("'main_assay' is not in assays of ace'.\n")
        stop(err)
      }
      main_mat = SummarizedExperiment::assays(ace)[[main_assay]]
      if (is.sparseMatrix(main_mat)) {
          write.HD5SpMat(h5file, gname = "X", main_mat, compression_level = compression_level)
      } else {
          h5file$create_dataset("X", main_mat, gzip_level = compression_level, dtype = h5types$H5T_IEEE_F32LE)
      }
    }

    remaining.assays = setdiff(names(SummarizedExperiment::assays(ace)), main_assay)
    if ((full.export == T) & (0 < length(remaining.assays))) {
        layers = h5file$create_group("layers")

        for (an in remaining.assays) {
          Xr = SummarizedExperiment::assays(ace)[[an]]
          if (is.sparseMatrix(Xr)) {
            write.HD5SpMat(layers, gname = an, Xr, compression_level = compression_level)
          } else {
            layers$create_dataset(an, Xr, gzip_level = compression_level, dtype = h5types$H5T_IEEE_F32LE)
          }
        }
    }

    uns = h5file$create_group("uns")
    obsm_annot = uns$create_group("obsm_annot")
    varm_annot = uns$create_group("varm_annot")

    obj_list = metadata(ace)
    write.HD5List(uns, "metadata", obj_list, depth = 1, max_depth = 10, compression_level = compression_level)

    ## Write obs (colData() in ace)
    obs.DF = as.data.frame(SummarizedExperiment::colData(ace))
    if (0 < NCOL(obs.DF)) {
        obs.DF = as.data.frame(lapply(SummarizedExperiment::colData(ace), function(x) {
            if (is.numeric(x) & (!is.null(names(x)))) {
                return(factor(names(x), names(x)[match(unique(x), x)]))
            } else {
                return(x)
            }
        }))
    }
    rownames(obs.DF) = colnames(ace)

    write.HD5DF(h5file, gname = "obs", obs.DF, compression_level = compression_level)

    ## Write var (matching rowData() in ace)
    var.DF = as.data.frame(SummarizedExperiment::rowData(ace))
    rownames(var.DF) = rownames(ace)
    if (class(SummarizedExperiment::rowRanges(ace)) == "GRanges") {
        GR = SummarizedExperiment::rowRanges(ace)
        BED = data.frame(chr = as.character(GenomeInfoDb::seqnames(GR)), start = start(GR), end = end(GR))
        var.DF = cbind(BED, var.DF)
    }
    write.HD5DF(h5file, "var", var.DF, compression_level = compression_level)

    ## Write subset of obsm related to the cell embeddings (Dim=2 or 3)
    obsm = h5file$create_group("obsm")
    obsm.mats = colMaps(ace)
    obsm.public.idx = which(colMapTypes(ace) != "internal")
    if (length(obsm.public.idx) > 0) {
        obsm.subset = obsm.mats[obsm.public.idx]
        for (i in 1:length(obsm.subset)) {
            nn = names(obsm.subset)[[i]]
            Y = Matrix::t(obsm.subset[[i]])
            if (NROW(Y) <= 3) {
                # AD_nn = paste("X", nn, sep = "_")
                AD_nn = nn
            } else {
                AD_nn = nn
            }

            if (is.matrix(Y)) {
                obsm$create_dataset(AD_nn, Y, gzip_level = compression_level, dtype = h5types$H5T_IEEE_F32LE)
            } else {
                write.HD5SpMat(obsm, AD_nn, Y, compression_level)
            }

            factor_info = obsm_annot$create_group(AD_nn)
            factor_info[["type"]] = colMapTypes(ace)[[nn]]
            factor.meta.DF = colMapMeta(ace)[[nn]]
            if (NCOL(factor.meta.DF) > 0) {
                write.HD5DF(factor_info, "annotatation", factor.meta.DF, compression_level = 0)
            }
        }
    }

    varm = h5file$create_group("varm")
    varm.mats = rowMaps(ace)
    varm.public.idx = which(rowMapTypes(ace) != "internal")
    if (length(varm.public.idx) > 0) {
        varm.subset = varm.mats[varm.public.idx]
        for (i in 1:length(varm.subset)) {
            nn = names(varm.subset)[[i]]
            Y = Matrix::t(varm.subset[[i]])

            if (is.matrix(Y)) {
                varm$create_dataset(nn, Y, gzip_level = compression_level, dtype = h5types$H5T_IEEE_F32LE)
            } else {
                write.HD5SpMat(varm, nn, Y, compression_level)
            }


            factor_info = varm_annot$create_group(nn)
            factor_info[["type"]] = rowMapTypes(ace)[[nn]]
            factor.meta.DF = rowMapMeta(ace)[[nn]]
            if (NCOL(factor.meta.DF) > 0) {
                write.HD5DF(factor_info, "annotatation", factor.meta.DF, compression_level = 0)
            }
        }
    }

    if (full.export) {
        print("Full export mode")

        if (length(obsm.public.idx) < length(obsm.mats)) {
            obsm.private.idx = setdiff(1:length(obsm.mats), obsm.public.idx)
            if (length(obsm.private.idx) > 0) {
                obsm.subset = obsm.mats[obsm.private.idx]
                for (i in 1:length(obsm.subset)) {
                  nn = names(obsm.subset)[[i]]
                  Y = Matrix::t(obsm.subset[[i]])
                  if (is.matrix(Y)) {
                    obsm$create_dataset(nn, Y, gzip_level = compression_level, dtype = h5types$H5T_IEEE_F32LE)
                  } else {
                    write.HD5SpMat(obsm, nn, Y, compression_level)
                  }

                  factor_info = obsm_annot$create_group(nn)
                  factor_info[["type"]] = colMapTypes(ace)[[nn]]
                  factor.meta.DF = colMapMeta(ace)[[nn]]
                  if (NCOL(factor.meta.DF) > 0) {
                    write.HD5DF(factor_info, "annotation", factor.meta.DF, compression_level = 0)
                  }
                }
            }
        }

        if (length(varm.public.idx) < length(varm.mats)) {
            varm.private.idx = setdiff(1:length(varm.mats), varm.public.idx)
            if (length(varm.private.idx) > 0) {
                varm.subset = varm.mats[varm.private.idx]
                for (i in 1:length(varm.subset)) {
                  nn = names(varm.subset)[[i]]
                  Y = Matrix::t(varm.subset[[i]])

                  if (is.matrix(Y)) {
                    varm$create_dataset(nn, Y, gzip_level = compression_level, dtype = h5types$H5T_IEEE_F32LE)
                  } else {
                    write.HD5SpMat(varm, nn, Y, compression_level)
                  }

                  factor_info = varm_annot$create_group(nn)
                  factor_info[["type"]] = rowMapTypes(ace)[[nn]]
                  factor.meta.DF = rowMapMeta(ace)[[nn]]
                  if (NCOL(factor.meta.DF) > 0) {
                    write.HD5DF(factor_info, "annotatation", factor.meta.DF, compression_level = 0)
                  }
                }
            }
        }

        # Export 'obsp'-associated matrices, i.e. colNets(): obs in AnnData ~ cols in SCE
        # ~ cells => cell-cell networks (such as ACTIONet)
        CN = colNets(ace)
        if ((length(CN) > 0)) {
            obsp = h5file$create_group("obsp")
            CN = lapply(CN, function(x) as(x, "dgCMatrix"))

            for (i in 1:length(CN)) {
                write.HD5SpMat(obsp, gname = names(CN)[[i]], CN[[i]], compression_level = compression_level)
            }
        }

        # Export 'varp'-associated matrices, i.e. rowNets(): var in AnnData ~ rows in SCE
        # ~ genes => gene-gene networks (such as SCINET)
        RN = rowNets(ace)
        if ((length(RN) > 0)) {
            varp = h5file$create_group("varp")
            RN = lapply(RN, function(x) as(x, "dgCMatrix"))

            for (i in 1:length(RN)) {
                write.HD5SpMat(varp, gname = names(RN)[[i]], RN[[i]], compression_level = compression_level)
            }
        }
    }

    h5file$close_all()
}

#' @import hdf5r
#' @export
AnnData2ACE <- function(
  file,
  main_assay = "X",
  import_X = TRUE
) {

    .check_and_load_package("hdf5r")

    h5file = H5File$new(file, mode = "r")

    objs = names(h5file)

    if(import_X == TRUE){
      X.attr = h5attributes(h5file[["X"]])
      if (length(X.attr) == 0) {
          # Full matrix
          X = h5file[["X"]]$read()
      } else {
          X = read.HD5SpMat(h5file = h5file, gname = "X")
      }

      input_assays = list(X)
      names(input_assays) = main_assay
    }

    if ("layers" %in% objs) {
        layers = h5file[["layers"]]
        additional_assays = vector("list", length(names(layers)))
        names(additional_assays) = names(layers)

        for (an in names(layers)) {
            attr = h5attributes(layers[[an]])
            if (length(attr) == 0) {
                # Dense matrix
                additional_assays[[an]] = layers[[an]]$read()
            } else {
                additional_assays[[an]] = read.HD5SpMat(h5file = layers, gname = an)
            }
        }
        input_assays = c(input_assays, additional_assays)
    }


    if ("obs" %in% objs) {
        obs.DF = read.HD5DF(h5file = h5file, gname = "obs")
    } else {
        obs.DF = DataFrame(row.names = paste("Cell", 1:NCOL(X), sep = ""))
    }

    if ("var" %in% objs) {
        var.DF = read.HD5DF(h5file = h5file, gname = "var")
    } else {
        var.DF = DataFrame(row.names = paste("Gene", 1:NROW(X), sep = ""))
    }

    input_assays = lapply(input_assays, function(X) {
        rownames(X) = rownames(var.DF)
        colnames(X) = rownames(obs.DF)
        return(X)
    })

    ace = ACTIONetExperiment(assays = input_assays, rowData = var.DF, colData = obs.DF)

    var.DF = rowData(ace)
    if (sum(colnames(var.DF) %in% c("chr", "start", "end")) == 3) {
        cols = match(c("chr", "start", "end"), colnames(obs.DF))
        BED = var.DF[, cols]
        var.DF = var.DF[, -cols]
        GR = GenomicRanges::makeGRangesFromDataFrame(BED)
        elementMetadata(GR) = var.DF

        SummarizedExperiment::rowRanges(ace) = GR
    }

    if ("obsm" %in% objs) {
        obsm = h5file[["obsm"]]
        for (mn in names(obsm)) {
            attr = h5attributes(obsm[[mn]])
            if ("encoding-type" %in% names(attr)) {
                if (attr[["encoding-type"]] %in% c("csc_matrix", "csr_matrix")) {
                  Xr = read.HD5SpMat(obsm, mn, compression_level)
                } else {
                  err = sprintf("Error reading obsm %s", mn)
                  h5file$close_all()
                  message(attr)
                  stop(msg)
                }
            } else {
                Xr = obsm[[mn]]$read()
            }

            if (sum(grepl(pattern = "^X_", mn))) {
                nn = stringr::str_sub(mn, start = 3)
            } else {
                nn = mn
            }
            colMaps(ace)[[nn]] = Matrix::t(Xr)
        }
    }

    if ("varm" %in% objs) {
        varm = h5file[["varm"]]
        for (nn in names(varm)) {
            attr = h5attributes(varm[[nn]])
            if ("encoding-type" %in% names(attr)) {
                if (attr[["encoding-type"]] %in% c("csc_matrix", "csr_matrix")) {
                  Xr = read.HD5SpMat(varm, nn, compression_level)
                } else {
                  err = sprintf("Error reading obsm %s", nn)
                  h5file$close_all()
                  message(attr)
                  stop(err)
                }
            } else {
                Xr = varm[[nn]]$read()
            }



            rowMaps(ace)[[nn]] = Matrix::t(Xr)
        }
    }

    if ("obsp" %in% objs) {
        obsp = h5file[["obsp"]]
        for (pn in names(obsp)) {
            Net = read.HD5SpMat(obsp, pn)
            colNets(ace)[[pn]] = Net
        }
    }

    if ("varp" %in% objs) {
        varp = h5file[["varp"]]
        for (pn in names(varp)) {
            Net = read.HD5SpMat(varp, pn)
            rowNets(ace)[[pn]] = Net
        }
    }


    if ("uns" %in% objs) {
        uns = h5file[["uns"]]
        if ("obsm_annot" %in% names(uns)) {
            # Import obs annotations
            obsm_annot = uns[["obsm_annot"]]
            for (nn in names(obsm_annot)) {
                factor_annot = obsm_annot[[nn]]
                colMapTypes(ace)[[nn]] = factor_annot[["type"]]$read()
                if ("annotation" %in% names(factor_annot)) {
                  DF = read.HD5DF(factor_annot, "annotation")
                  colMapMeta(ace)[[nn]] = DF
                }
            }
        }
        if ("varm_annot" %in% names(uns)) {
            # Import obs annotations
            var_annot = uns[["varm_annot"]]
            for (nn in names(var_annot)) {
                factor_annot = var_annot[[nn]]
                rowMapTypes(ace)[[nn]] = factor_annot[["type"]]$read()
                if ("annotation" %in% names(factor_annot)) {
                  DF = read.HD5DF(factor_annot, "annotation")
                  rowMapMeta(ace)[[nn]] = DF
                }
            }
        }
        if ("metadata" %in% names(uns)) {
            meta_data_objs = read.HD5List(uns, "metadata", depth = 1, max_depth = 10,
                compression_level = compression_level)
        } else {
            meta_data_objs = NULL
        }



        meta = read.HD5List(h5file = h5file, gname = "uns")
        meta = meta[setdiff(names(meta), c("obsm_annot", "varm_annot", "metadata"))]
        if (!is.null(meta_data_objs)) {
            if (length(meta) == 0) {
                meta = meta_data_objs
            } else {
                meta = c(meta_data_objs, meta)
            }
        }
        metadata(ace) = meta
    }

    h5file$close_all()

    return(ace)
}
