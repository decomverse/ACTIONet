# From: SnapATAC
CPal88 = c("#E31A1C", "#FFD700", "#771122", "#777711", "#1F78B4", "#68228B", "#AAAA44", "#60CC52", "#771155", "#DDDD77", 
    "#774411", "#AA7744", "#AA4455", "#117744", "#000080", "#44AA77", "#AA4488", "#DDAA77", "#D9D9D9", "#BC80BD", "#FFED6F", "#7FC97F", 
    "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
    "#E6AB02", "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", 
    "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", 
    "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFFF33", "#A65628", 
    "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", 
    "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5")

# From: Scanpy
CPal20 = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31')
    

preprocess.labels <- function(labels, ace = NULL) {
	if(is.null(labels)) {
		return(NULL)
	}
	if( (length(labels) == 1) & is.character(labels)) {
		if(is.null(ace)) {
			R.utils::printf('Error preprocess.labels: annotation.name %s not found (no ace object provided)\n', labels)
			return(NULL)
		}
		idx = which(names(colData(ace)) == labels)
		if(length(idx) == 0) {
			R.utils::printf('Error preprocess.labels: annotation.name %s not found\n', labels)
			return(NULL)
		}		
		labels = colData(ace)[, idx]
	}

	if((length(labels) > 1) & is.logical(labels)) {
		labels = factor(as.numeric(labels), levels = c(0, 1), labels = c("No", "Yes"))
	}
	
	if((length(labels) > 1) & is.character(labels)) {
		labels = factor(labels)
	}
	
	if(is.factor(labels)) {
		v = as.numeric(labels)
		names(v) = levels(labels)[v]		
		labels = v
	}
	if(is.matrix(labels)) {
		L = as.numeric(labels)
		names(L) = names(labels)
		labels = L
	} 
	
	if( is.null(names(labels)) | length(unique(names(labels))) > 100 )  {
		names(labels) = as.character(labels)
	}
	
	return(labels)
}

layout.labels <- function(x, y, labels, col = "white", bg = "black", r = 0.1, cex = 1.0, ...) {
	require(wordcloud)
    lay <- wordlayout(x, y, words = labels, cex = 1.25*cex, ...)

    x = lay[, 1] + 0.5 * lay[, 3]
    y = lay[, 2] + 0.5 * lay[, 4]
    
    theta = seq(0, 2 * pi, length.out = 50)
    xy <- xy.coords(x, y)
    xo <- r * strwidth("A")
    yo <- r * strheight("A")
    
    for (i in theta) {
        text(xy$x + cos(i) * xo, xy$y + sin(i) * yo, labels, col = bg, cex = cex, ...)
    }
    text(xy$x, xy$y, labels, col = col, cex = cex, ...)
}

#' Main ACTIONet plotting functions
#'
#' @param ace ACTIONet output object
#' @param labels Annotation of interest (clusters, celltypes, etc.) to be projected on the ACTIONet plot
#' @param transparency.attr Additional continuous attribute to project onto the transparency of nodes
#' @param trans.z.threshold, trans.fact Control the effect of transparency mapping
#' @param node.size Size of nodes in the ACTIONet plot
#' @param CPal Color palette (named vector or a name for a given known palette)
#' @param add.text Whether or not to add labels on the top of ACTIONet plot
#' @param suppress.legend Whether to suppress legend or include one
#' @param legend.pos If suppress.legend == F, where should we put the legend
#' @param title Main title of the plot
#' @param border.contrast.factor How much the node and its border should contrast
#' @param add.states Whether or not to include interpolated cell state positions in the plot
#' @param coordinate.slot Entry in reducedDims(ace) containing the plot coordinates (default:"ACTIONet2D")
#' 
#' @return Visualized ACTIONet
#' 
#' @examples
#' ace = run.ACTIONet(sce)
#' plot.ACTIONet(ace, ace$assigned_archetype, transparency.attr = ace$node_centrality)
plot.ACTIONet <- function(ace, labels = NULL, transparency.attr = NULL, trans.z.threshold = -0.5, trans.fact = 1, 
	node.size = 1, CPal = CPal20, add.text = TRUE, suppress.legend = TRUE, legend.pos = "bottomright", title = "", border.contrast.factor = 0.1, coordinate.slot = "ACTIONet2D") {
    
    text.halo.width = 0.1
    label.text.size = 0.8
    
    node.size = node.size * 0.25
    
    if(class(ace) == "ACTIONetExperiment") {
		labels = preprocess.labels(labels, ace)
		if(is.character(coordinate.slot)) {
			coors = reducedDims(ace)[[coordinate.slot]]
		} else {
			coors = as.matrix(coordinate.slot)
		}
	} else {
		if(is.matrix(ace) | is.sparseMatrix(ace)) {
			coors = as.matrix(ace)
			labels = preprocess.labels(labels)
		}
		else {
			print("Unknown type for ace")
			return()
		}
	}
	if(is.null(labels)) {
		if(class(ace) == "ACTIONetExperiment") {			
			vCol = rgb(reducedDims(ace)$denovo_color)
		} else {
			vCol = rep("tomato", nrow(coors))
		}
		Annot = NULL
	} else {
		Annot = names(labels)[match(sort(unique(labels)), labels)]
		if(length(CPal) > 1) {
			if(length(CPal) < length(Annot)) {
				if(length(Annot) <= 20) {
					CPal = CPal20
					message("Not enough colors. Switching to CPal20")
				} else {
					CPal = CPal88
					message("Not enough colors. Switching to CPal88")
				}
			}            
            if(is.null(names(CPal))) {					
				Pal = CPal[1:length(Annot)]
			} else {            
				Pal = CPal[Annot]		
			}
		} else {
            Pal = ggpubr::get_palette(CPal, length(Annot))
		}
			
        names(Pal) = Annot
        vCol = Pal[names(labels)]
	}



    if (!is.null(transparency.attr)) {
        z = scale(transparency.attr) # (transparency.attr - median(transparency.attr))/mad(transparency.attr)
        beta = 1/(1 + exp(-trans.fact * (z - trans.z.threshold)))
        beta[z > trans.z.threshold] = 1
        beta = beta^trans.fact
        
        vCol.border = scales::alpha(colorspace::darken(vCol, border.contrast.factor), beta)
        vCol = scales::alpha(vCol, beta)
    } else {
        vCol.border = colorspace::darken(vCol, border.contrast.factor)
    }

	x = coors[, 1]
	y = coors[, 2]
	x.min = min(x)
	x.max = max(x)
	y.min = min(y)
	y.max = max(y)
	x.min = x.min - (x.max-x.min)/20
	x.max = x.max + (x.max-x.min)/20
	y.min = y.min - (y.max-y.min)/20
	y.max = y.max + (y.max-y.min)/20
	XL = c(x.min, x.max)
	YL = c(y.min, y.max)


	rand.perm = sample(nrow(coors))
	graphics::plot(coors[rand.perm, c(1, 2)], pch = 21, cex = node.size, bg = vCol[rand.perm], col = vCol.border[rand.perm], axes = F, xlab = "", ylab = "", main = title, xlim = XL, ylim = YL) 

    if ( add.text == T & (!is.null(Annot)) ) {
		par(xpd = T, mar = par()$mar * c(1.1,1.1,1.1,1.1))
		
    	require(wordcloud)
        centroids = t(sapply(Annot, function(l) {
            idx = which(names(labels) == l)
            if(length(idx) == 1) {
				return(as.numeric(coors[idx, ]))
			} 

            sub.coors = coors[idx, ]
            anchor.coor = as.numeric(apply(sub.coors, 2, function(x) mean(x, trim = 0.80)))
#             sub.coors.sq = sub.coors^2
# 			norm.sq = Matrix::rowSums(sub.coors.sq)
# 			anchor.idx = which.min(sapply(1:nrow(sub.coors.sq), function(i) { 
# 				dd = norm.sq[i] + norm.sq - 2* sub.coors %*% sub.coors[i, ]
# 				mean.dist.sq = median(dd)
# 				return(mean.dist.sq)
# 			}))
            
                        
            # D = as.matrix(dist(sub.coors))
            # stats = Matrix::rowMeans(D)
            # anchor.idx = which.min(stats)
            
			# anchor.coor = as.numeric(sub.coors[anchor.idx, ])            
            
			return(anchor.coor)
        }))
        layout.labels(x = centroids[, 1], y = centroids[, 2], labels = Annot, col = colorspace::darken(Pal, 0.5), bg = "#eeeeee", r = text.halo.width, cex = label.text.size)
        #wordcloud::textplot(x = centroids[, 1], y = centroids[, 2], new = F, words = Annot, col = colorspace::darken(Pal, 0.5), bg = "#eeeeee") 
    }
    
    if ( (suppress.legend == FALSE) & !is.null(Annot) ) {
		xmin <- par("usr")[1]
		xmax <- par("usr")[2]
		ymin <- par("usr")[3]
		ymax <- par("usr")[4]

		lgd <- legend(x = mean(c(xmin,xmax)), y =  mean(c(ymin,ymax)), legend = Annot, fill = Pal, cex = 0.5, bty = "n", plot = F)

    	par(xpd = T, mai = c(0, 0, 0, lgd$rect$w))

		legend(x = xmax, y = ymin+lgd$rect$h, legend = Annot, fill = Pal, cex = 0.5, bty = "n", plot = T)

		#legend("bottom", legend = Annot, fill = Pal, cex = 0.5, plot = T)
    }    
}

#' Main ACTIONet 3D plotting functions
#'
#' @param ace ACTIONet output object
#' @param labels Annotation of interest (clusters, celltypes, etc.) to be projected on the ACTIONet plot
#' @param transparency.attr Additional continuous attribute to project onto the transparency of nodes
#' @param trans.z.threshold, trans.fact Control the effect of transparency mapping
#' @param node.size Size of nodes in the ACTIONet plot
#' @param CPal Color palette (named vector or a name for a given known palette)
#' @param title Main title of the plot
#' @param coordinate.slot Entry in reducedDims(ace) containing the plot coordinates (default:"ACTIONet2D")
#' 
#' @return Visualized ACTIONet
#' 
#' @examples
#' ace = run.ACTIONet(sce)
#' plot.ACTIONet.3D(ace, ace$assigned_archetype, transparency.attr = ace$node_centrality)
plot.ACTIONet.3D <- function(ace, labels = NULL, transparency.attr = NULL, trans.z.threshold = -1, trans.fact = 1, node.size = 1, CPal = CPal20, coordinate.slot = "ACTIONet3D") {
    require(ggplot2)
    require(ggpubr)
    require(threejs)
    

    nV = length(ncol(ace))
    
    node.size = node.size * 0.2
    

    if(class(ace) == "ACTIONetExperiment") {
		labels = preprocess.labels(labels, ace)
		if(is.character(coordinate.slot)) {
			coors = reducedDims(ace)[[coordinate.slot]]
		} else {
			coors = as.matrix(coordinate.slot)
		}
	} else {
		if(is.matrix(ace) | is.sparseMatrix(ace)) {
			coors = as.matrix(ace)
			labels = preprocess.labels(labels)
		}
		else {
			print("Unknown type for ace")
			return()
		}
	}
	    
	if(is.null(labels)) {
		if(class(ace) == "ACTIONetExperiment") {			
			vCol = rgb(reducedDims(ace)$denovo_color)
		} else {
			vCol = rep("tomato", nrow(coors))
		}
		Annot = NULL
	} else {
		Annot = names(labels)[match(sort(unique(labels)), labels)]
		if(length(CPal) > 1) {
			if(length(CPal) < length(Annot)) {
				if(length(Annot) <= 20) {
					CPal = CPal20
					message("Not enough colors. Switching to CPal20")
				} else {
					CPal = CPal88
					message("Not enough colors. Switching to CPal88")
				}
			}
			if(is.null(names(CPal))) {					
				Pal = CPal[1:length(Annot)]
			} else {            
				Pal = CPal[Annot]		
			}
		} else {
            Pal = ggpubr::get_palette(CPal, length(Annot))
		}

        names(Pal) = Annot
        vCol = Pal[names(labels)]
	}



    if (!is.null(transparency.attr)) {
        z = scale(transparency.attr) # (transparency.attr - median(transparency.attr))/mad(transparency.attr)
        beta = 1/(1 + exp(-trans.fact * (z - trans.z.threshold)))
        beta[z > trans.z.threshold] = 1
        beta = beta^trans.fact
        
        vCol.border = scales::alpha(colorspace::darken(vCol, 0.5), beta)
        vCol = scales::alpha(vCol, beta)
    } else {
        vCol.border = colorspace::darken(vCol, 0.5)
    }
        
    scatterplot3js(x = coors[, 1], y = coors[, 2], z = coors[, 3], axis.scales = FALSE, size = node.size, axis = F, grid = F, color = as.character(vCol), 
        stroke = as.character(vCol.border), bg = "black")
}

#' Plots heatmap of the top-ranked features of an enrichment table
#'
#' @param feature.enrichment.table An arbitrary enrichment table with columns being archetypes and rows being specificity of features
#' @param top.features Number of features to return
#' @param reorder.columns Whether to optimally re-order columns of the enrichment table
#' 
#' @return Enrichment heatmap
#' 
#' @examples
#' feature.enrichment.table = as.matrix(rowFactors(ace)[["H_unified_upper_significance"]])
#' plot.top.k.features(feature.enrichment.table, 3)
plot.top.k.features <- function(feature.enrichment.table, top.features = 3, normalize = T, reorder.columns = T, row.title = "Archetypes", column.title = "Genes", rowPal = "black") {

	W = select.top.k.features(feature.enrichment.table, top.features = 3, normalize = T)
	
	gradPal =  grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)

	M = apply(W, 1, max)
	
	Z = Matrix::t(apply(W, 1, function(w) (w-min(w)) / (max(w) - min(w))))
	# Z = t(scale(t(W)))
	# Z[Z > 3] = 3
	# Z[Z < -3] = -3
	ht =  Heatmap(Z, name = "Expression (scaled)", cluster_rows = F, cluster_columns = F, col = gradPal,
					row_title = row.title, column_title = column.title,
					column_names_gp = gpar(fontsize = 8, fontface = "bold"),
					row_names_gp = gpar(fontsize = 8, fontface = "bold", col = rowPal), column_title_gp = gpar(fontsize = 14,
					  fontface = "bold"), row_title_gp = gpar(fontsize = 14,
					  fontface = "bold"), row_names_side = "left", rect_gp = gpar(col = "black"))
	
	return(ht)
}


#' Plots projection of features onto the ACTIONet plot
#'
#' @param ace ACTIONet output
#' @param feature.enrichment.table An arbitrary enrichment table with columns being archetypes and rows being specificity of features
#' @param top.features Number of features to return
#' @param reorder.columns Whether to optimally re-order columns of the enrichment table
#' @param CPal Color palette to use
#' @param title Main title of the plot
#' 
#' @return Featur view
#' 
#' @examples
#' feature.enrichment.table = as.matrix(rowFactors(ace)[["H_unified_upper_significance"]])
#' plot.ACTIONet.feature.view(ace, feature.enrichment.table, 5)
plot.ACTIONet.feature.view <- function(ace, feature.enrichment.table, top.features = 5, CPal = NULL, title = "Feature view", label.text.size = 1, renormalize = F) {
	M = Matrix::t(as(colFactors(ace)[["H_unified"]], 'sparseMatrix'))
	cs = Matrix::colSums(M)
	M = scale(M, center = FALSE, scale = cs)
	
	if(ncol(feature.enrichment.table) != nrow(colFactors(ace)[["H_unified"]])) {
		feature.enrichment.table = Matrix::t(feature.enrichment.table)
	}
	
	if(max(feature.enrichment.table) > 50)
		feature.enrichment.table = log1p(feature.enrichment.table)

	X = t(select.top.k.features(feature.enrichment.table, 3, normalize = renormalize, reorder.columns = F))
	selected.features = colnames(X)
	
	core.coors = t(t(reducedDims(ace)[["ACTIONet2D"]]) %*% M)
	cs = colSums(X)
	cs[cs == 0] = 1
	X = scale(X, center = F, scale = cs)
	feature.coors = scale(t(X) %*% core.coors)

    if (is.null(CPal)) {
        #Pal = ace$unification.out$Pal
			cells.Lab = grDevices::convertColor(color = reducedDims(ace)$denovo_color, from = "sRGB", to = "Lab")
			arch.Lab = Matrix::t(M) %*% cells.Lab			
			arch.RGB = grDevices::convertColor(color = arch.Lab, from = "Lab", to = "sRGB")
			core.Pal = rgb(arch.RGB)
    } else {
    	if(length(CPal) == 1) {
            core.Pal = ggpubr::get_palette(CPal, length(unique(ace$archetype.assignment)))
    	} else {
            core.Pal = CPal[1:length(unique(ace$archetype.assignment))]
    	}
    }
    core.Lab = grDevices::convertColor(color = t(col2rgb(core.Pal)/256), from = "sRGB", to = "Lab")
    
    feature.color.Lab = t(X) %*% core.Lab
    feature.colors = rgb(grDevices::convertColor(color = feature.color.Lab, from = "Lab", to = "sRGB"))
    names(feature.colors) = selected.features


	x = feature.coors[, 1]
	y = feature.coors[, 2]
    plot (x, y, type = "n", col = feature.colors, axes = FALSE, xlab = "", ylab = "", main = title)

	words = selected.features
    lay <- wordlayout(x, y, words, label.text.size)
    for (i in 1:length(x)) {
        xl <- lay[i, 1]
        yl <- lay[i, 2]
        w <- lay[i, 3]
        h <- lay[i, 4]
        if (x[i] < xl || x[i] > xl + w || y[i] < yl || y[i] > 
            yl + h) {
            points(x[i], y[i], pch = 16, col = colorspace::darken(feature.colors[[i]], 0.6), cex = 0.75*label.text.size)
            nx <- xl + 0.5 * w
            ny <- yl + 0.5 * h
            lines(c(x[i], nx), c(y[i], ny), col = colorspace::darken(feature.colors[[i]], 0.5))
        }
    }
    #plot.new()
    loc.x = lay[, 1] + 0.5 * lay[, 3]
    loc.y = lay[, 2] + 0.5 * lay[, 4]
    text(loc.x, loc.y, words, col = feature.colors, cex = label.text.size)

}

#' Plots projection of genes onto the ACTIONet plot
#'
#' @param ace ACTIONet output
#' @param top.genes Number of genes to return
#' @param blacklist.pattern List of genes to filter-out
#' @param reorder.columns Whether to optimally re-order columns of the enrichment table
#' @param CPal Color palette to use
#' @param title Main title of the plot
#' 
#' @return Featur view
#' 
#' @examples
#' plot.ACTIONet.gene.view(ace, 5)
plot.ACTIONet.gene.view <- function(ace, top.genes = 5, CPal = NULL, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH", title = "", label.text.size = 0.8, renormalize = F) {
	require(wordcloud)
	
	feature.enrichment.table = as.matrix(rowFactors(ace)[["H_unified_upper_significance"]])
	
	filtered.rows = grep(blacklist.pattern, rownames(feature.enrichment.table))
	if(length(filtered.rows) > 0)
		feature.enrichment.table = feature.enrichment.table[-filtered.rows, ]
	
	plot.ACTIONet.feature.view(ace, feature.enrichment.table, title = "Gene view", renormalize = renormalize)
}


#' Main ACTIONet interactive plotting function
#'
#' @param ace ACTIONet output object
#' @param labels Annotation of interest (clusters, celltypes, etc.) to be projected on the ACTIONet plot
#' @param transparency.attr Additional continuous attribute to project onto the transparency of nodes
#' @param trans.z.threshold, trans.fact Control the effect of transparency mapping
#' @param node.size Size of nodes in the ACTIONet plot
#' @param CPal Color palette (named vector or a name for a given known palette)
#' @param enrichment.table To project the top-ranked features interactively.
#' @param top.features Number of features to show per cell
#' @param blacklist.pattern List of genes to filter-out
#' @param title Main title of the plot
#' @param threeD Whether to show the plot in 3D
#' 
#' @return Visualized ACTIONet
#' 
#' @examples
#' ace = run.ACTIONet(sce)
#' plot.ACTIONet.interactive(ace, ace$assigned_archetype)
plot.ACTIONet.interactive <- function(ace, labels = NULL, transparency.attr = NULL, trans.z.threshold = -1, trans.fact = 1, 
	require(igraph)
	require(plotly)
	
	node.size = 1, CPal = CPal20, enrichment.table = NULL, top.features = 7, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH", threeD = FALSE, title = "ACTIONet", coordinate.slot = "ACTIONet2D") {
    require(plotly)
    require(ACTIONet)

    nV = ncol(ace)
    node.size = node.size * 3
    if(coordinate.slot == "ACTIONet2D" & threeD == T)
    	coordinate.slot = "ACTIONet3D"
    
    if(class(ace) == "ACTIONetExperiment") {
		labels = preprocess.labels(labels, ace)
		if(is.character(coordinate.slot)) {
			coors = reducedDims(ace)[[coordinate.slot]]
		} else {
			coors = as.matrix(coordinate.slot)
		}
	} else {
		if(is.matrix(ace) | is.sparseMatrix(ace)) {
			coors = as.matrix(ace)
			labels = preprocess.labels(labels)
		}
		else {
			print("Unknown type for ace")
			return()
		}
	}
		
	if(is.null(labels)) {
		if(class(ace) == "ACTIONetExperiment") {			
			vCol = rgb(reducedDims(ace)$denovo_color)
		} else {
			vCol = rep("tomato", nrow(coors))
		}
		Annot = NULL
	} else {
		Annot = names(labels)[match(sort(unique(labels)), labels)]
		if(length(CPal) > 1) {
			if(length(CPal) < length(Annot)) {
				if(length(Annot) <= 20) {
					CPal = CPal20
					message("Not enough colors. Switching to CPal20")
				} else {
					CPal = CPal88
					message("Not enough colors. Switching to CPal88")
				}
			}
			if(is.null(names(CPal))) {					
				Pal = CPal[1:length(Annot)]
			} else {            
				Pal = CPal[Annot]		
			}
		} else {
            Pal = ggpubr::get_palette(CPal, length(Annot))
		}
		
        names(Pal) = Annot
        vCol = Pal[names(labels)]
	}

    if (!is.null(transparency.attr)) {
        z = scale(transparency.attr) # (transparency.attr - median(transparency.attr))/mad(transparency.attr)
        beta = 1/(1 + exp(-trans.fact * (z - trans.z.threshold)))
        beta[z > trans.z.threshold] = 1
        beta = beta^trans.fact
        
        vCol.border = scales::alpha(colorspace::darken(vCol, 0.5), beta)
        vCol = scales::alpha(vCol, beta)
    } else {
        vCol.border = colorspace::darken(vCol, 0.5)
    }

    
	
	if( !is.null(enrichment.table) ) {
		if(ncol(enrichment.table) == nV) {
			cell.scores = Matrix::t(enrichment.table)
		} else if( (nrow(enrichment.table) != nV) ) {
			H = colFactors(ace)[["H_unified"]]
			if( (nrow(enrichment.table) == nrow(H)) | (ncol(enrichment.table) == nrow(H)) ) {
				cell.scores = map.cell.scores.from.archetype.enrichment(ace, enrichment.table)				
			} else {
				cell.scores = NULL
			}
		} else {
			cell.scores = enrichment.table
		}
	} else {
		temp.enrichment.table = as.matrix(rowFactors(ace)[["H_unified_upper_significance"]])			
		if( !is.null(row.names(temp.enrichment.table)) ) {
			filtered.rows = grep(blacklist.pattern, rownames(temp.enrichment.table))
			if(length(filtered.rows) > 0)
				enrichment.table = temp.enrichment.table[-filtered.rows, ]
			else
				enrichment.table = temp.enrichment.table

			GT = apply(enrichment.table, 2, function(x) rownames(enrichment.table)[order(x, decreasing = T)[1:min(100, nrow(enrichment.table))]])
			selected.features = sort(unique(as.character(GT)))
			
			cell.scores = Matrix::t(enrichment.table[selected.features, ] %*% colFactors(ace)[["H_unified"]])
		} else {
			cell.scores = NULL
		}	
	}

    if ( !is.null(cell.scores) ) {
		selected.features = colnames(cell.scores)
		node.annotations = apply(cell.scores, 1, function(x) paste(selected.features[order(x, decreasing = T)[1:top.features]], collapse = '\n'))
		# node.annotations = sapply(1:length(ace$log$cells), function(i) sprintf('(%s) %s', ace$log$cells[[i]], node.annotations[[i]]) )
	} else {
		# node.annotations = sapply(1:length(ace$log$cells), function(i) sprintf('(%s)', ace$log$cells[[i]]) )
		node.annotations = rep('', nV)
	}
    
    # Setup visualization parameters
    sketch.graph = graph_from_adjacency_matrix(colNets(ace)$ACTIONet, mode = "undirected", weighted = TRUE)
    
    if(threeD == F) {
		V(sketch.graph)$x = coors[, 1]
		V(sketch.graph)$y = coors[, 2]
	} else {
		V(sketch.graph)$x3D = coors[, 1]
		V(sketch.graph)$y3D = coors[, 2]
		V(sketch.graph)$z3D = coors[, 3]
	}
	
    sketch.graph = delete.edges(sketch.graph, E(sketch.graph))
    
    node.data <- get.data.frame(sketch.graph, what = "vertices")
    edge.data <- get.data.frame(sketch.graph, what = "edges")
    
    Nv <- dim(node.data)[1]
    Ne <- dim(edge.data)[1]
    
    edge_shapes <- list()
    
    # Adjust parameters
    node.data$size = node.size
    
    axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)

    if(threeD == TRUE) {
		if(is.null(Annot)) {
		    node.data$vCol = vCol
			node.data$vCol.border = vCol.border
			network <- plot_ly(node.data, x = ~x3D, y = ~y3D, z = ~z3D, opacity = 1, marker = list(color = ~vCol, size = ~size, opacity = 1, alpha = 1, line = list(width = 0.1 * node.size, alpha = 0.5, color = ~vCol.border)), text = node.annotations, mode = "markers", hoverinfo = "text", type = "scatter3d", showlegend = FALSE)		
			p <- plotly::layout(network, title = title, shapes = edge_shapes, acene = list(xaxis=axis,yaxis=axis,zaxis=axis))
			
		} else {
			# border.Pal = colorspace::darken(Pal, 0.5)
			# names(border.Pal) = names(Pal)
			node.data$vCol.border = vCol.border

		    node.data$type = factor(names(labels), levels = Annot)
		    
			network <- plot_ly(node.data, x = ~x3D, y = ~y3D, z = ~z3D, opacity = 1, color = ~type, colors = Pal, marker = list(size = ~size, opacity = 1, alpha = 1, line = list(width = 0.1 * node.size, alpha = 0.5, color = ~vCol.border)), text = node.annotations, mode = "markers", hoverinfo = "text", type = "scatter3d")		
			p <- plotly::layout(network, title = title, shapes = edge_shapes, acene = list(xaxis=axis,yaxis=axis,zaxis=axis), showlegend = TRUE, legend = list(marker = list(marker.size = 10)))
		}		
	} else {				
		if(is.null(Annot)) {
		    node.data$vCol = vCol
			node.data$vCol.border = vCol.border
			network <- plot_ly(node.data, x = ~x, y = ~y, marker = list(color = ~vCol, size = ~size, opacity = 1, alpha = 1, line = list(width = 0.1 * node.size, alpha = 0.5, color = ~vCol.border)), text = node.annotations, mode = "markers", type = "scatter", hoverinfo = "text", showlegend = FALSE)
			p <- plotly::layout(network, title = title, shapes = edge_shapes, xaxis = axis, yaxis = axis)
		} else {
			# border.Pal = colorspace::darken(Pal, 0.5)
			# names(border.Pal) = names(Pal)
			node.data$vCol.border = vCol.border

		    node.data$type = factor(names(labels), levels = Annot)
		    
			network <- plot_ly(node.data, x = ~x, y = ~y, color = ~type, colors = Pal, marker = list(size = ~size, line = list(width = 0.1 * node.size, color = ~vCol.border)), text = node.annotations, mode = "markers", type = "scatter", hoverinfo = "text")
			p <- plotly::layout(network, title = title, shapes = edge_shapes, xaxis = axis, yaxis = axis, showlegend = TRUE, legend = list(marker = list(marker.size = 10)))
		}
	}
		
    p
    
}


#' Plot gene expression violin plot
#'
#' @param ace ACTIONet output object
#' @param labels Annotation of interest (clusters, celltypes, etc.) to be projected on the ACTIONet plot
#' @param gene.name Name of the gene to plot
#' @param CPal Color palette (named vector or a name for a given known palette)
#' 
#' @return Visualized ACTIONet
#' 
#' @examples
#' plot.individual.gene(ace, ace$assigned_archetype, "CD14")
plot.individual.gene <- function(ace, labels, gene.name, CPal = CPal20) {
    require(igraph)
    require(ACTIONet)
    require(ggpubr)

	clusters = preprocess.labels(ace, labels)

	Labels = names(clusters)
	Annot = sort(unique(Labels))
	Annot = Annot[order(clusters[match(Annot, Labels)], decreasing = F)]
	Labels = factor(Labels, levels = Annot)
	
	if(length(CPal) > 1) {
		if(length(CPal) < length(Annot)) {
			if(length(Annot) <= 20) {
				CPal = CPal20
				message("Not enough colors. Switching to CPal20")
			} else {
				CPal = CPal88
				message("Not enough colors. Switching to CPal88")
			}
		}		
        if(is.null(names(CPal))) {
				
			Pal = CPal[1:length(Annot)]
		} else {            
			Pal = CPal[Annot]		
		}
	} else {
        Pal = ggpubr::get_palette(CPal, length(Annot))
	}
	
    names(Pal) = Annot

	if( !(gene.name %in% rownames(ace)) ) {
		R.utils::printf("Gene %s not found\n", gene.name)
		return()
	}
	
	x = SummarizedExperiment::assays(ace)$logcounts[gene.name, ]
	if(sum(x) == 0) {
		print("Gene is not expressed")
		return()
	}
	
	
	require(ggpubr)
	require(ggplot2)
	df = data.frame(Annotation = Labels, Expression = x)	
    gp = ggviolin(df, "Annotation", "Expression", fill = "Annotation", palette = Pal, add = "boxplot", add.params = list(fill = "white"))    
    print(gp)
}


#' Projects a given continuous score on the ACTIONet plot
#'
#' @param ace ACTIONet output object
#' @param x Score vector
#' @param transparency.attr Additional continuous attribute to project onto the transparency of nodes
#' @param trans.z.threshold, trans.fact Control the effect of transparency mapping
#' @param node.size Size of nodes in the ACTIONet plot
#' @param CPal Color palette (named vector or a name for a given known palette)
#' @param coordinate.slot Entry in reducedDims(ace) containing the plot coordinates (default:"ACTIONet2D")
#' @param alpha_val Between [0, 1]. If it is greater than 0, smoothing of scores would be performed
#' 
#' @return Visualized ACTIONet with projected scores
#' 
#' @examples
#' ace = run.ACTIONet(sce)
#' x = logcounts(ace)["CD14", ]
#' plot.ACTIONet.gradient(ace, x, transparency.attr = ace$node_centrality)
plot.ACTIONet.gradient <- function(ace, x, transparency.attr = NULL, trans.z.threshold = -0.5, trans.fact = 3, node.size = 1, CPal = "magma", title = "", alpha_val = 0.85, nonparameteric = FALSE, coordinate.slot = "ACTIONet2D") {

    node.size = node.size * 0.3

    if(class(ace) == "ACTIONetExperiment") {
		if(is.character(coordinate.slot)) {
			coors = reducedDims(ace)[[coordinate.slot]]
		} else {
			coors = as.matrix(coordinate.slot)
		}
	} else {
		alpha_val = 0
		if(is.matrix(ace) | is.sparseMatrix(ace)) {
			coors = as.matrix(ace)
		}
		else {
			print("Unknown type for ace")
			return()
		}
	}    
    

    NA.col = "#eeeeee"
        
    ## Create color gradient generator
    if (CPal %in% c("greys", "inferno", "magma", "viridis", "BlGrRd", "RdYlBu", "Spectral")) {
		require(viridis)
        Pal_grad = switch(CPal, greys = gray.colors(100), inferno = viridis::inferno(500, alpha = 0.8), magma = viridis::magma(500, alpha = 0.8), viridis = viridis::viridis(500, alpha = 0.8), 
            BlGrRd = colorRampPalette(c("blue", "grey", "red"))(500), Spectral = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, 
                name = "Spectral"))))(100), RdYlBu = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"))))(100))
    } else {
        Pal_grad = colorRampPalette(c(NA.col, CPal))(500)
    }	

    ## Scale/prune scorees, if needed
    x[x < 0] = 0
	if(max(x) > 50)
		x = log1p(x)

	if(alpha_val > 0) {
		x = as.numeric(compute_network_diffusion(colNets(ace)$ACTIONet, as(x, 'sparseMatrix')))
	}
	
    if (nonparameteric == TRUE) {
        vCol = scales::col_bin(Pal_grad, domain = NULL, na.color = NA.col, bins = 7)(rank(x))
    } else {
        vCol = scales::col_bin(Pal_grad, domain = NULL, na.color = NA.col, bins = 7)(x)
    }	

    if (!is.null(transparency.attr)) {
        z = scale(transparency.attr) # (transparency.attr - median(transparency.attr))/mad(transparency.attr)
        beta = 1/(1 + exp(-trans.fact * (z - trans.z.threshold)))
        beta[z > trans.z.threshold] = 1
        beta = beta^trans.fact
        
        vCol = scales::alpha(vCol, beta)
        vCol.border = scales::alpha(colorspace::darken(vCol, 0.1), beta)
    } else {
        vCol.border = colorspace::darken(vCol, 0.1)
    }
	
    idx = order(x, decreasing = F)
    plot(coors[idx, 1], coors[idx, 2], bg = vCol[idx], col = vCol.border[idx], cex = node.size, pch = 21, axes = FALSE, xlab = "", ylab = "", main = title)
	
}


#' Projects a set of markers on the ACTIONet
#' It also optionally imputes the markers.
#'
#' @param ace ACTIONet output object
#' @param x Score vector
#' @param transparency.attr Additional continuous attribute to project onto the transparency of nodes
#' @param trans.z.threshold, trans.fact Control the effect of transparency mapping
#' @param node.size Size of nodes in the ACTIONet plot
#' @param CPal Color palette (named vector or a name for a given known palette)
#' @param coordinate.slot Entry in reducedDims(ace) containing the plot coordinates (default:"ACTIONet2D")
#' @param alpha_val Between [0, 1]. If it is greater than 0, smoothing of scores would be performed
#' 
#' @return Visualized ACTIONet with projected scores
#' 
#' @examples
#' ace = run.ACTIONet(sce)
#' x = logcounts(ace)["CD14", ]
#' visualize.markers(ace, c("CD14", "CD19", "CD3G"), transparency.attr = ace$node_centrality)
visualize.markers <- function(ace, marker.genes, transparency.attr = NULL, trans.z.threshold = -0.5, trans.fact = 3, node.size = 1, CPal = "magma",  alpha_val = 0.95) {
    require(igraph)
    
    
    if (!sum(sapply(marker.genes, length) != 1) & is.null(names(marker.genes))) {
        names(marker.genes) = marker.genes
    }
    
    gg = intersect(unique(unlist(marker.genes)), rownames(ace))
    all.marker.genes = sort(intersect(gg, rownames(ace)))
    if(length(all.marker.genes) == 0) {
		return()
	}
    
    if(alpha_val > 0)
		imputed.marker.expression = impute.genes.using.ACTIONet(ace, all.marker.genes, alpha_val = alpha_val)
	else {
		imputed.marker.expression = Matrix::t(assays(ace)[["logcounts"]])
	}
    
    lapply(all.marker.genes, function(gene) {
        print(gene)
        if (!(gene %in% colnames(imputed.marker.expression))) 
            return()
        
        idx = which(sapply(marker.genes, function(gs) gene %in% gs))[1]
        celltype.name = names(marker.genes)[idx]
        
        x = imputed.marker.expression[, gene]

		plot.ACTIONet.gradient(ace, x, transparency.attr, trans.z.threshold, trans.fact, node.size, CPal = CPal, title = gene, alpha_val = 0)
    })
}


select.top.k.genes <- function(ace, top.genes = 5, CPal = NULL, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH", top.features = 3, normalize = F, reorder.columns = F, specificity.slot.name = "H_unified_upper_significance") {	
	feature.enrichment.table = as.matrix(rowFactors(ace)[[specificity.slot.name]])
	
	filtered.rows = grep(blacklist.pattern, rownames(feature.enrichment.table))
	if(length(filtered.rows) > 0)
		feature.enrichment.table = feature.enrichment.table[-filtered.rows, ]
	
	tbl = select.top.k.features(feature.enrichment.table, top.features = top.features, normalize = normalize, reorder.columns = reorder.columns)
	
	return(tbl)
}


plot.top.k.genes <- function(ace, top.genes = 5, CPal = NULL, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH", top.features = 3, normalize = F, reorder.columns = T, row.title = "Archetypes", column.title = "Genes", rowPal = "black", specificity.slot.name = "H_unified_upper_significance") {
	require(ComplexHeatmap)
	
	feature.enrichment.table = as.matrix(rowFactors(ace)[[specificity.slot.name]])
	
	filtered.rows = grep(blacklist.pattern, rownames(feature.enrichment.table))
	if(length(filtered.rows) > 0)
		feature.enrichment.table = feature.enrichment.table[-filtered.rows, ]
	
	ht = plot.top.k.features(feature.enrichment.table, top.features = top.features, normalize = normalize, reorder.columns = reorder.columns, row.title = row.title, column.title = column.title, rowPal = rowPal)
	
	return(ht)
	
}


plot.archetype.selected.genes <- function(ace, genes, CPal = NULL, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH", top.features = 3, normalize = F, reorder.columns = T, row.title = "Archetypes", column.title = "Genes", rowPal = "black", specificity.slot.name = "H_unified_upper_significance") {
	feature.enrichment.table = as.matrix(rowFactors(ace)[["H_unified_upper_significance"]])
	
	filtered.rows = match(intersect(rownames(ace), genes), rownames(ace))
	if(length(filtered.rows) > 0)
		feature.enrichment.table = feature.enrichment.table[-filtered.rows, ]
	
	ht = plot.top.k.features(feature.enrichment.table, top.features = top.features, normalize = normalize, reorder.columns = reorder.columns, row.title = row.title, column.title = column.title, rowPal = rowPal)
	
	return(ht)
	
}
plot.ACTIONet.archetype.footprint <- function(ace, node.size = 1, CPal = "magma", title = "", arch.labels = NULL, coordinate.slot = "ACTIONet2D", alpha_val = 0.9) {
	Ht = Matrix::t(colFactors(ace)[["H_unified"]])
	cs = Matrix::colSums(Ht)
	cs[cs == 0] = 1
	U = as(scale(Ht, center = F, scale = cs), 'dgTMatrix')
	U.pr = compute_network_diffusion(colNets(ace)$ACTIONet, U, alpha = alpha_val)
	
	node.size = node.size * 0.3
    coors = reducedDims(ace)[[coordinate.slot]]

    if (CPal %in% c("inferno", "magma", "viridis", "BlGrRd", "RdYlBu", "Spectral")) {
		require(viridis)
        Pal_grad = switch(CPal, inferno = viridis::inferno(500, alpha = 0.8), magma = viridis::magma(500, alpha = 0.8), viridis = viridis::viridis(500, alpha = 0.8), 
            BlGrRd = colorRampPalette(c("blue", "grey", "red"))(500), Spectral = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, 
                name = "Spectral"))))(100), RdYlBu = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"))))(100))
    } else {
		NA.col = "#cccccc"
        Pal_grad = colorRampPalette(c(NA.col, CPal))(500)
    }	


	k1 = k2 = round(sqrt(ncol(Ht)))
	if(k1*k2 < ncol(Ht)) {
		k2 = k2+1
	}

	if(is.null(arch.labels))
		arch.labels = sapply(1:ncol(Ht), function(i) sprintf("Archetype %d", i))
	 
	#par(mfrow = c(k1, k2), mar = c(0, 0, 1, 0))
	sapply(1:ncol(Ht), function(i) {
		print(i)
		x = U.pr[, i]

		xs = sort(x, decreasing = T)
		nnz = round( (sum(xs)^2) / (sum(xs^2)) )
		threshold = xs[nnz]

		x[x < threshold] = threshold
		x = log(x)
		


		
	    vCol = scales::col_bin(Pal_grad, domain = NULL, bins = 10)(x)
        vCol = scales::alpha(vCol, 0.05 + 0.95*x/max(x))
        vCol = colorspace::lighten(vCol, 0.2)
	    idx = order(x, decreasing = F)
		plot(coors[idx, 1], coors[idx, 2], bg = vCol[idx], col = vCol[idx], cex = node.size, pch = 21, axes = FALSE, xlab = "", ylab = "", main = arch.labels[[i]])
	    
	})
}



#' Report the top-rank features from a given enrichment table
#'
#' @param top.features Number of features to return
#' @param reorder.columns Whether to optimally re-order columns of the enrichment table
#' 
#' @return Sorted table with the selected top-ranked
#' 
#' @examples
#' feature.enrichment.table = as.matrix(rowFactors(ace)[["H_unified_upper_significance"]])
#' enrichment.table.top = select.top.k.features(feature.enrichment.table, 3)
select.top.k.features <- function(feature.enrichment.table, top.features = 3, normalize = F, reorder.columns = T) {
	
	W0 = (feature.enrichment.table)
	if(normalize == T)
		W0 = doubleNorm(W0)
	
	IDX = matrix(0, nrow = top.features, ncol = ncol(W0))
	VV = matrix(0, nrow = top.features, ncol = ncol(W0))
	W = (W0)
	for(i in 1:nrow(IDX)) {
	  W.m = as(MWM_hungarian(W), 'dgTMatrix')
	  IDX[i, W.m@j+1] = W.m@i + 1
	  VV[i, W.m@j+1] = W.m@x
	  
	  W[IDX[i, W.m@j+1], ] = 0
	}

	if(reorder.columns == T) {
		feature.enrichment.table.aggregated = apply(IDX, 2, function(perm) as.numeric(Matrix::colMeans(W0[perm, ])))
		CC = cor(feature.enrichment.table.aggregated)
		D = as.dist(1 - CC)
		cols = seriation::get_order(seriation::seriate(D, "OLO"))	
		rows = as.numeric(IDX[, cols])
	} else {
		cols = 1:ncol(W0)
		rows = unique(as.numeric(IDX))
	}
	W = feature.enrichment.table[rows, cols]
	
	return(W)
}


