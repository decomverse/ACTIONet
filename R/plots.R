# From: https://github.com/r3fang/SnapATAC/blob/master/R/plottings-utilities.R
Pals.Snap88 = c("#E31A1C", "#FFD700", "#771122", "#777711", "#1F78B4", "#68228B", "#AAAA44", "#60CC52", "#771155", "#DDDD77", 
    "#774411", "#AA7744", "#AA4455", "#117744", "#000080", "#44AA77", "#AA4488", "#DDAA77", "#D9D9D9", "#BC80BD", "#FFED6F", "#7FC97F", 
    "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
    "#E6AB02", "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", 
    "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", 
    "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFFF33", "#A65628", 
    "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", 
    "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5")



plot.ACTIONet <- function(ACTIONet.out, labels = NULL, transparency.attr = NULL, trans.z.threshold = -0.5, trans.fact = 2, 
	node.size = 1, CPal = Pals.Snap88, add.text = TRUE, text.halo.width = 0.1, label.text.size = 0.8, 
    suppress.legend = TRUE, legend.pos = "bottomright", add.states = F, title = "", highlight = NULL, border.contrast.factor = 0.5) {
    
    node.size = node.size * 0.5
    
    if (is.igraph(ACTIONet.out)) 
        ACTIONet = ACTIONet.out else ACTIONet = ACTIONet.out$ACTIONet
    
    coors = cbind(V(ACTIONet)$x, V(ACTIONet)$y)

	if( (length(labels) == 1) & (is.character(labels)) ) {
		annotation_name = labels
	} else {
		annotation_name = "plot.ACTIONet-tmp"
	}	
    
	labels = preprocess.labels(ACTIONet.out, labels)

	if(!is.null(highlight)) {
		if( !is.null(labels) & (tolower(highlight) == "labels") ) {
			if( annotation_name == "plot.ACTIONet-tmp") {
				ACTIONet.out = add.cell.annotations(ACTIONet.out, cell.annotations = labels, annotation.name = annotation_name)
			}
			
			if( is.null(ACTIONet.out$annotations[[annotation_name]]$highlight) ) {
				print("Labels are not highlighted ... generating highlight on the fly")
				ACTIONet.out = highlight.annotations(ACTIONet.out, annotation.name = annotation_name)			
			}
			label.hightlight = ACTIONet.out$annotations[[annotation_name]]$highlight

			transparency.attr = label.hightlight$connectivity.scores		
		} else if(tolower(highlight) == "connectivity") {
			transparency.attr = V(ACTIONet.out$ACTIONet)$connectivity
		}
	}
	
	if(is.null(labels)) {
        vCol = V(ACTIONet)$color
        Annot = NULL
	} else {
		Annot = names(labels)[match(sort(unique(labels)), labels)]
		if(length(CPal) > 1) {
			#idx = labels[match(sort(unique(labels)), labels)]
            #Pal = CPal[idx]	
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

	if(add.states == T) {
		par(new=TRUE)
		
		M = as(ACTIONet.out$unification.out$C.core, 'sparseMatrix')
		cs = Matrix::colSums(M)
		M = scale(M, center = FALSE, scale = cs)
		
	    cell.Lab = grDevices::convertColor(color = t(col2rgb(vCol)/256), from = "sRGB", to = "Lab")	    
	    core.Lab = t(t(cell.Lab) %*% M)
	    core.colors = rgb(grDevices::convertColor(color = core.Lab, from = "Lab", to = "sRGB"))
		core.colors = colorspace::lighten(core.colors, 0.1)
		
		core.coors = t(t(coors) %*% M)		

	    graphics::plot(core.coors, pch = 25, add = T, cex = 3*node.size, bg = core.colors, col = "#eeeeee", axes = FALSE, xlab = "", ylab = "", xlim = XL, ylim = YL)
	    
        layout.labels(x = core.coors[, 1], y = core.coors[, 2]-strheight("A"), labels = 1:ncol(M), col = colorspace::darken(core.colors, 0.5), bg = "#eeeeee", r = text.halo.width, cex = label.text.size)

	}    
    
    
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


plot.ACTIONet.3D <- function(ACTIONet.out, labels = NULL, transparency.attr = NULL, trans.z.threshold = -1, trans.fact = 2, node.size = 1, CPal = Pals.Snap88, highlight = NULL) {
    require(ggplot2)
    require(ggpubr)
    require(threejs)
    
    if (is.igraph(ACTIONet.out)) 
        ACTIONet = ACTIONet.out else ACTIONet = ACTIONet.out$ACTIONet
    
    nV = length(V(ACTIONet))
    coor = cbind(V(ACTIONet)$x3D, V(ACTIONet)$y3D, V(ACTIONet)$z3D)
    
    node.size = node.size * 0.2
    

	if( (length(labels) == 1) & (is.character(labels)) ) {
		annotation_name = labels
	} else {
		annotation_name = "plot.ACTIONet-tmp"
	}	
	    
	labels = preprocess.labels(ACTIONet.out, labels)
	
	if(!is.null(highlight)) {
		if( !is.null(labels) & (tolower(highlight) == "labels") ) {
			if( annotation_name == "plot.ACTIONet-tmp") {
				ACTIONet.out = add.cell.annotations(ACTIONet.out, cell.annotations = labels, annotation.name = annotation_name)
			}
			
			if( is.null(ACTIONet.out$annotations[[annotation_name]]$highlight) ) {
				print("Labels are not highlighted ... generating highlight on the fly")
				ACTIONet.out = highlight.annotations(ACTIONet.out, annotation.name = annotation_name)			
			}
			label.hightlight = ACTIONet.out$annotations[[annotation_name]]$highlight

			transparency.attr = label.hightlight$connectivity.scores		
		} else if(tolower(highlight) == "connectivity") {
			transparency.attr = V(ACTIONet.out$ACTIONet)$connectivity
		}
	}
	
	
	if(is.null(labels)) {
        vCol = V(ACTIONet)$color
        add.text = F
        suppress.legend = T
	} else {
		Annot = names(labels)[match(sort(unique(labels)), labels)]
		if(length(CPal) > 1) {
			idx = labels[match(sort(unique(labels)), labels)]
            Pal = CPal[idx]			
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
        
    scatterplot3js(x = coor[, 1], y = coor[, 2], z = coor[, 3], axis.scales = FALSE, size = node.size, axis = F, grid = F, color = as.character(vCol), 
        stroke = as.character(vCol.border), bg = "black")
}




plot.ACTIONet.feature.view <- function(ACTIONet.out, feature.enrichment.table, top.features = 5, CPal = NULL, title = "Feature view", label.text.size = 1) {
	if(ncol(feature.enrichment.table) != nrow(ACTIONet.out$unification.out$H.core)) {
		feature.enrichment.table = Matrix::t(feature.enrichment.table)
	}
	
	if(max(feature.enrichment.table) > 50)
		feature.enrichment.table = log1p(feature.enrichment.table)

	feature.enrichment.table = doubleNorm(feature.enrichment.table)
	
	selected.features = sort(unique(as.character(apply(feature.enrichment.table, 2, function(x) rownames(feature.enrichment.table)[order(x, decreasing = T)[1:top.features]]))))
	
	M = Matrix::t(as(ACTIONet.out$unification.out$H.core, 'sparseMatrix'))
	cs = Matrix::colSums(M)
	M = scale(M, center = FALSE, scale = cs)
	
	core.coors = t(t(ACTIONet.out$vis.out$coordinates) %*% M)
	X = t(feature.enrichment.table[selected.features, ])
	cs = colSums(X)
	cs[cs == 0] = 1
	X = scale(X, center = F, scale = cs)
	feature.coors = t(X) %*% core.coors

    if (is.null(CPal)) {
        Pal = ACTIONet.out$unification.out$Pal
    } else {
    	if(length(CPal) == 1) {
            Pal = ggpubr::get_palette(CPal, length(ACTIONet.out$unification.out$Pal))
    	} else {
            Pal = CPal[1:length(ACTIONet.out$unification.out$Pal)]
    	}
    }
    
    core.Lab = grDevices::convertColor(color = t(col2rgb(Pal)/256), from = "sRGB", to = "Lab")
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

plot.ACTIONet.gene.view <- function(ACTIONet.out, top.genes = 5, CPal = NULL, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH", title = "", label.text.size = 0.8) {
    if( !('unification.out' %in% names(ACTIONet.out)) ) {
		print('Error in plot.ACTIONet.gene.view: "unification.out" is not in ACTIONet.out. Please run unify.cell.states() first.')
		return()
	}
	gene.enrichment.table = as.matrix(SummarizedExperiment::assays(ACTIONet.out$unification.out$DE.core)[["significance"]])
	
	filtered.rows = grep(blacklist.pattern, rownames(gene.enrichment.table))
	if(length(filtered.rows) > 0)
		gene.enrichment.table = gene.enrichment.table[-filtered.rows, ]
	
	GT = apply(gene.enrichment.table, 2, function(x) rownames(gene.enrichment.table)[order(x, decreasing = T)[1:top.genes]])
	selected.genes = sort(unique(as.character(GT)))
	
	M = Matrix::t(as(ACTIONet.out$unification.out$H.core, 'sparseMatrix'))
	cs = Matrix::colSums(M)
	M = scale(M, center = FALSE, scale = cs)
	
	core.coors = t(t(ACTIONet.out$vis.out$coordinates) %*% M)
	X = t(gene.enrichment.table[selected.genes, ])
	cs = colSums(X)
	cs[cs == 0] = 1
	X = scale(X, center = F, scale = cs)
	gene.coors = t(X) %*% core.coors

    if (is.null(CPal)) {
        Pal = ACTIONet.out$unification.out$Pal
    } else {
    	if(length(CPal) == 1) {
            Pal = ggpubr::get_palette(CPal, length(ACTIONet.out$unification.out$Pal))
    	} else {
            Pal = CPal[1:length(ACTIONet.out$unification.out$Pal)]
    	}
    }
    
    core.Lab = grDevices::convertColor(color = t(col2rgb(Pal)/256), from = "sRGB", to = "Lab")
    gene.color.Lab = t(X) %*% core.Lab
    gene.colors = rgb(grDevices::convertColor(color = gene.color.Lab, from = "Lab", to = "sRGB"))
    names(gene.colors) = selected.genes

	x = gene.coors[, 1]
	y = gene.coors[, 2]
	x.min = min(x)
	x.max = max(x)
	y.min = min(y)
	y.max = max(y)
	x.min = x.min - (x.max-x.min)/10
	x.max = x.max + (x.max-x.min)/10
	y.min = y.min - (y.max-y.min)/10
	y.max = y.max + (y.max-y.min)/10
	
	words = selected.genes
	plot(x = NULL, y = NULL, xlim = c(x.min, x.max), ylim = c(y.min, y.max), axes = FALSE, xlab = "", ylab = "", main = title)
    lay <- wordlayout(x, y, words, 1.5*label.text.size)
	# x0 <- lay[, 1]
	# y0 <- lay[, 2]
	# w0 <- lay[, 3]
	# h0 <- lay[, 4]
	# x.min = min(x0-w0/2)
	# x.max = max(x0+w0/2)
	# y.min = min(y0-h0/2)
	# y.max = max(y0+h0/2)
    
    for (i in 1:length(x)) {
        xl <- lay[i, 1]
        yl <- lay[i, 2]
        w <- lay[i, 3]
        h <- lay[i, 4]
        if (x[i] < xl || x[i] > xl + w || y[i] < yl || y[i] > 
            yl + h) {
            points(x[i], y[i], pch = 16, col = colorspace::darken(gene.colors[[i]], 0.35), cex = 0.75*label.text.size)
            nx <- xl + 0.5 * w
            ny <- yl + 0.5 * h
            lines(c(x[i], nx), c(y[i], ny), col = colorspace::darken(gene.colors[[i]], 0.3))
        }
    }
    text(lay[, 1] + 0.5 * lay[, 3], lay[, 2] + 0.5 * lay[, 4], words, col = gene.colors, cex = label.text.size, xlab = "", ylab = "", main = title)
    
}

plot.ACTIONet.interactive <- function(ACTIONet.out, labels = NULL, transparency.attr = NULL, trans.z.threshold = -1, trans.fact = 2, 
	node.size = 1, CPal = Pals.Snap88, enrichment.table = NULL, top.features = 7, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH", threeD = FALSE, title = "ACTIONet") {
    require(plotly)
    require(ACTIONet)
    
    if (is.igraph(ACTIONet.out)) 
        ACTIONet = ACTIONet.out else ACTIONet = ACTIONet.out$ACTIONet
    
    nV = length(V(ACTIONet))
    node.size = node.size * 3
    
	labels = preprocess.labels(ACTIONet.out, labels)
	if(is.null(labels)) {
        vCol = V(ACTIONet)$color
        Annot = NULL
	} else {
		Annot = names(labels)[match(sort(unique(labels)), labels)]
		if(length(CPal) > 1) {
			idx = labels[match(sort(unique(labels)), labels)]
            Pal = CPal[idx]			
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
			H = ACTIONet.out$unification.out$H.core
			if( (nrow(enrichment.table) == nrow(H)) | (ncol(enrichment.table) == nrow(H)) ) {
				cell.scores = map.cell.scores.from.archetype.enrichment(ACTIONet.out, enrichment.table)				
			} else {
				cell.scores = NULL
			}
		} else {
			cell.scores = enrichment.table
		}
	} else {
		if( ('unification.out' %in% names(ACTIONet.out)) ) {
			temp.enrichment.table = as.matrix(SummarizedExperiment::assays(ACTIONet.out$unification.out$DE.core)[["significance"]])			
			if( !is.null(row.names(temp.enrichment.table)) ) {
				filtered.rows = grep(blacklist.pattern, rownames(temp.enrichment.table))
				if(length(filtered.rows) > 0)
					enrichment.table = temp.enrichment.table[-filtered.rows, ]
				else
					enrichment.table = temp.enrichment.table

				GT = apply(enrichment.table, 2, function(x) rownames(enrichment.table)[order(x, decreasing = T)[1:min(100, nrow(enrichment.table))]])
				selected.features = sort(unique(as.character(GT)))
				
				cell.scores = t(enrichment.table[selected.features, ] %*% ACTIONet.out$unification.out$H.core)
			} else {
				cell.scores = NULL
			}
		} else {
			cell.scores = NULL
		}	
	}

    if ( !is.null(cell.scores) ) {
		selected.features = colnames(cell.scores)
		node.annotations = apply(cell.scores, 1, function(x) paste(selected.features[order(x, decreasing = T)[1:top.features]], collapse = '\n'))
		# node.annotations = sapply(1:length(ACTIONet.out$log$cells), function(i) sprintf('(%s) %s', ACTIONet.out$log$cells[[i]], node.annotations[[i]]) )
	} else {
		# node.annotations = sapply(1:length(ACTIONet.out$log$cells), function(i) sprintf('(%s)', ACTIONet.out$log$cells[[i]]) )
		node.annotations = rep('', nV)
	}
    
    # Setup visualization parameters
    sketch.graph = ACTIONet.out$ACTIONet
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
			p <- plotly::layout(network, title = title, shapes = edge_shapes, scene = list(xaxis=axis,yaxis=axis,zaxis=axis))
			
		} else {
			# border.Pal = colorspace::darken(Pal, 0.5)
			# names(border.Pal) = names(Pal)
			node.data$vCol.border = vCol.border

		    node.data$type = factor(names(labels), levels = Annot)
		    
			network <- plot_ly(node.data, x = ~x3D, y = ~y3D, z = ~z3D, opacity = 1, color = ~type, colors = Pal, marker = list(size = ~size, opacity = 1, alpha = 1, line = list(width = 0.1 * node.size, alpha = 0.5, color = ~vCol.border)), text = node.annotations, mode = "markers", hoverinfo = "text", type = "scatter3d")		
			p <- plotly::layout(network, title = title, shapes = edge_shapes, scene = list(xaxis=axis,yaxis=axis,zaxis=axis), showlegend = TRUE, legend = list(marker = list(marker.size = 10)))
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

plot.individual.gene <- function(ACTIONet.out, annotation.name, sce, gene.name, imputation = F, alpha_val = 0.85, node.size = 3, CPal = Pals.Snap88) {
    require(igraph)
    require(ACTIONet)
    require(ggpubr)
	cl.idx = which(names(ACTIONet.out$annotations) == annotation.name)
	if(length(cl.idx) == 0) {
		R.utils::printf('Error in plot.annotations.individual.gene: annotation.name "%s" not found\n', annotation.name)
		return()
	}	

	clusters = ACTIONet.out$annotations[[cl.idx]]$Labels
	Labels = names(clusters)
	Annot = sort(unique(Labels))
	Annot = Annot[order(clusters[match(Annot, Labels)], decreasing = F)]
	Labels = factor(Labels, levels = Annot)
	
	if(length(CPal) > 1) {
        Pal = CPal[1:length(Annot)]			
	} else {
        Pal = ggpubr::get_palette(CPal, length(Annot))
	}
    names(Pal) = Annot

	if( !(gene.name %in% rownames(sce)) ) {
		R.utils::printf("Gene %s not found\n", gene.name)
		return()
	}
	
	x = SummarizedExperiment::assays(sce)$logcounts[gene.name, ]
	if(sum(x) == 0) {
		print("Gene is not expressed")
		return()
	}
	
	if(imputation == T) {
		x = x / sum(x)
		x = igraph::page.rank(ACTIONet.out$ACTIONet, personalized = x, damping = alpha_val)$vector
	}
	x = x / mean(x)
	
	require(ggpubr)
	require(ggplot2)
	df = data.frame(Annotation = Labels, Expression = x)	
    gp = ggboxplot(df, x = "Annotation", y = "Expression", fill = "Annotation", palette = Pal) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    print(gp)
}


plot.ACTIONet.gradient <- function(ACTIONet.out, x, transparency.attr = NULL, trans.z.threshold = -0.5, trans.fact = 3, node.size = 1, CPal = "magma", title = "", alpha_val = 0.85, nonparameteric = FALSE, highlight = NULL) {

    node.size = node.size * 0.3

    if (is.igraph(ACTIONet.out)) 
        ACTIONet = ACTIONet.out else ACTIONet = ACTIONet.out$ACTIONet
    coors = cbind(V(ACTIONet)$x, V(ACTIONet)$y)

	if(!is.null(highlight)) {
		if(tolower(highlight) == "connectivity") {
			transparency.attr = V(ACTIONet.out$ACTIONet)$connectivity
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
		x = as.numeric(PageRank_iter(ACTIONet.out$build.out$ACTIONet, as(x, 'sparseMatrix')))
	}
	
    if (nonparameteric == TRUE) {
        vCol = scales::col_numeric(Pal_grad, domain = NULL, na.color = NA.col)(rank(x))
    } else {
        vCol = scales::col_numeric(Pal_grad, domain = NULL, na.color = NA.col)(x)
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

visualize.markers <- function(ACTIONet.out, sce, marker.genes, transparency.attr = NULL, trans.z.threshold = -0.5, trans.fact = 3, node.size = 1, CPal = "magma",  export_path = NA, alpha_val = 0.85, highlight = NULL) {
    require(igraph)
    
    
    if (!sum(sapply(marker.genes, length) != 1) & is.null(names(marker.genes))) {
        names(marker.genes) = marker.genes
    }
    
    gg = intersect(unique(unlist(marker.genes)), rownames(sce))
    all.marker.genes = sort(intersect(gg, rownames(sce)))
    if(length(all.marker.genes) == 0) {
		return()
	}
    
    if(alpha_val > 0)
		imputed.marker.expression = impute.genes.using.ACTIONet(ACTIONet.out, sce, all.marker.genes, alpha_val = alpha_val)
	else {
		imputed.marker.expression = Matrix::t(assays(sce)[["logcounts"]])
	}
    
    lapply(all.marker.genes, function(gene) {
        print(gene)
        if (!(gene %in% colnames(imputed.marker.expression))) 
            return()
        
        idx = which(sapply(marker.genes, function(gs) gene %in% gs))[1]
        celltype.name = names(marker.genes)[idx]
        
        x = imputed.marker.expression[, gene]

		plot.ACTIONet.gradient(ACTIONet.out, x, transparency.attr, trans.z.threshold, trans.fact, node.size, CPal = CPal, title = gene, alpha_val = alpha_val, highlight = highlight)
    })
}

#  [1] "ARSA"          "BBURCG"        "BBWRCG"        "GW"            "GW_average"    "GW_complete"   "GW_single"    
#  [8] "GW_ward"       "HC"            "HC_average"    "HC_complete"   "HC_single"     "HC_ward"       "Identity"     
# [15] "MDS"           "MDS_angle"     "MDS_metric"    "MDS_nonmetric" "OLO"           "OLO_average"   "OLO_complete" 
# [22] "OLO_single"    "OLO_ward"      "QAP_2SUM"      "QAP_BAR"       "QAP_Inertia"   "QAP_LS"        "R2E"          
# [29] "Random"        "SA"            "Spectral"      "Spectral_norm" "SPIN_NH"       "SPIN_STS"      "TSP"          
# [36] "VAT"     
plot.enrichment.list <- function(Enrichment.list, row.title, seriation.method = "OLO_complete", scale.rows = T, shared.columns = F, name = "z-score", sort.rows = T, sort.cols = T, rowPal = NULL, border.col = "black") {
	require(ComplexHeatmap)
	require(seriation)
	require(RColorBrewer)
	
	if(is.null(names(Enrichment.list))) {
		names(Enrichment.list) = 1:length(Enrichment.list)
	}
	
	if(sort.rows == T) {
		set.seed(0)
		CC = Reduce("+", lapply(Enrichment.list, function(E) cor(as.matrix(Matrix::t(E))))) / length(Enrichment.list)
		CC[is.na(CC)] = 0
		D = as.dist(1-CC)
		row.perm = seriation::get_order(seriation::seriate(D, seriation.method))
	} else {
		row.perm = 1:nrow(Enrichment.list[[1]])
	}
	
	if(sort.cols == T) {
		if(shared.columns == T) {
			set.seed(0)
			CC = Reduce("+", lapply(Enrichment.list, function(E) cor(as.matrix(E)))) / length(Enrichment.list)
			CC[is.na(CC)] = 0
			D = as.dist(1-CC)
			col.perm = seriation::get_order(seriation::seriate(D, seriation.method))
			col.perms = lapply(1:length(Enrichment.list), function(i) return(col.perm))
		} else {
			col.perms = lapply(Enrichment.list, function(E) {
				set.seed(0)
				CC = cor(as.matrix(E))
				CC[is.na(CC)] = 0		
				D = as.dist(1-CC)
				col.perm = seriation::get_order(seriation::seriate(D, seriation.method))
				return(col.perm)
			})
		}
	} else {
		col.perms = lapply(Enrichment.list, function(E) {
			1:ncol(E)
		})		
	}
	
	gradPal = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
	
	if(sort.rows) {
		rowPal = rowPal[row.perm]
	}	
		
	ht_list = NULL  ## Heatmap(...) + NULL gives you a HeatmapList object
	for(i in 1:length(Enrichment.list)) {
		Enrichment = Enrichment.list[[i]]
		if(scale.rows == T) {
			Enrichment = Matrix::t(scale(Matrix::t(Enrichment)))
		}
		if(sort.cols) {
			Enrichment = Enrichment[, col.perms[[i]]]
		}
		if(sort.rows) {
			Enrichment = Enrichment[row.perm, ]
		}
		if(is.null(rowPal)) {
			ht_list = ht_list + Heatmap(Enrichment, name = name, cluster_rows = F, cluster_columns = F, col = gradPal, row_title = row.title, column_title = names(Enrichment.list)[[i]], column_names_gp = gpar(fontsize = 10, fontface="bold"), row_names_gp = gpar(fontsize = 10, fontface="bold"), column_title_gp = gpar(fontsize = 12, fontface="bold"), row_title_gp = gpar(fontsize = 12, fontface="bold"), row_names_side = "left", rect_gp = gpar(col = border.col), row_names_max_width = unit(10, "cm"), column_names_max_height = unit(10, "cm"))
		} else {
			ht_list = ht_list + Heatmap(Enrichment, name = name, cluster_rows = F, cluster_columns = F, col = gradPal, row_title = row.title, column_title = names(Enrichment.list)[[i]], column_names_gp = gpar(fontsize = 10, fontface="bold"), row_names_gp = gpar(fontsize = 10, fontface="bold", col = rowPal), column_title_gp = gpar(fontsize = 12, fontface="bold"), row_title_gp = gpar(fontsize = 12, fontface="bold"), row_names_side = "left", rect_gp = gpar(col = border.col), row_names_max_width = unit(10, "cm"), column_names_max_height = unit(10, "cm"))
		}
	}
	
	draw(ht_list)
}

plot.annotated.heatmap <-function(W, row.annotation, column.annotation, plot_title = NA, row_title = NA, colnames_title = NA, CPal = "d3") {
  modules = V(W)$module

  Annot = sort(unique(modules))
	modCelltype.Pal = ggpubr::get_palette("d3", length(Annot))
  names(modCelltype.Pal) = Annot
    
  Annot = sort(unique(rownames(W)))
	if(1 < length(CPal)) {
		rowCelltype.Pal = CPal[1:length(Annot)]
	} else {
		rowCelltype.Pal = ggpubr::get_palette(CPal, length(Annot))
	}
	names(rowCelltype.Pal) = Annot
  
  ha_rows = HeatmapAnnotation(df = list(Celltype = rownames(W), Module = modules), col = list(Celltype = rowCelltype.Pal, Module = modCelltype.Pal), annotation_legend_param = list(Celltype=list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 5)), Module=list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 5))), which = "row")


  
  Annot = sort(unique(colnames(W)))
	if(1 < length(CPal)) {
		colCelltype.Pal = CPal[1:length(Annot)]
	} else {
		colCelltype.Pal = ggpubr::get_palette(CPal, length(Annot))
	}
	names(colCelltype.Pal) = Annot
  
  ha_cols = HeatmapAnnotation(df = list(Celltype = colnames(W), Module = modules), col = list(Celltype = colCelltype.Pal, Module = modCelltype.Pal), annotation_legend_param = list(Celltype=list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 5)), Module=list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 5))), which = "column")  
  

  gradPal = grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = "YlOrRd"))(100)

  Heatmap(W, col = gradPal, row_names_gp = gpar(fontsize = 0), column_names_gp = gpar(fontsize = 0), left_annotation = ha_rows, top_annotation = ha_cols, name = "Correlation", row_title = ds1.name, column_title = ds1.name, cluster_rows = FALSE, cluster_columns = FALSE)
}


plot.annotations.selected.genes <- function(ACTIONet.out, sce, annotation.name, markers, type = "heatmap", seriation.method = "OLO", sort.rows = T) {
	cl.idx = which(names(ACTIONet.out$annotations) == annotation.name)
	if(length(cl.idx) == 0) {
		R.utils::printf('Error in plot.annotations.differential.heatmap: annotation.name "%s" not found\n', annotation.name)
		return()
	}		
	
	Labels = ACTIONet.out$annotations[[cl.idx]]$Labels
	UL = sort(unique(Labels))
	Annot = names(Labels)[match(UL, Labels)]
	Labels = factor(Labels, levels = UL, labels = Annot)

	IDX = split(1:length(Labels), Labels)


	A = assays(sce)$logcounts
	Avg.profile = sapply(IDX, function(idx) Matrix::rowMeans(A[markers, idx]))

	if(sort.rows == T) {
		set.seed(0)
		CC = cor(as.matrix(Matrix::t(A[markers, ])))
		CC[is.na(CC)] = 0
		D = as.dist(1-CC)
		row.perm = seriation::get_order(seriation::seriate(D, seriation.method))
	} else {
		row.perm = length(markers)
	}
	
	set.seed(0)
	CC = cor(as.matrix(Avg.profile))
	CC[is.na(CC)] = 0
	D = as.dist(1-CC)
	col.perm = seriation::get_order(seriation::seriate(D, seriation.method))
	
	gradPal = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
	
	Enrichment = Matrix::t(scale(Matrix::t(Avg.profile[row.perm, col.perm])))

	if(type == "dotplot") {
		library(corrplot)
		corrplot::corrplot(Enrichment, is.corr = F, tl.col = "black", tl.cex = 0.7, cl.pos = "n")
	} else if(type == "heatmap") {
		library(ComplexHeatmap)
		ComplexHeatmap::Heatmap(Enrichment, name = "z-score", cluster_rows = F, cluster_columns = F, col = gradPal, row_title = "", column_title = annotation.name, column_names_gp = gpar(fontsize = 8, fontface="bold"), row_names_gp = gpar(fontsize = 8, fontface="bold"), column_title_gp = gpar(fontsize = 10, fontface="bold"), row_names_side = "left", rect_gp = gpar(col = "black"))		
	}	
		
}

plot.annotations.top.genes <- function(ACTIONet.out, sce, annotation.name, top.genes, type = "heatmap", seriation.method = "OLO") {
	cl.idx = which(names(ACTIONet.out$annotations) == annotation.name)
	if(length(cl.idx) == 0) {
		R.utils::printf('Error in plot.annotations.differential.heatmap: annotation.name "%s" not found\n', annotation.name)
		return()
	}		
    if (is.null(ACTIONet.out$annotations[[cl.idx]]$DE.profile)) {
		R.utils::printf('Error in plot.annotations.differential.heatmap: annotation.name "%s" does not DE.profile. Please run compute.annotations.feature.specificity() first.\n', annotation.name)
		return()
    }
    X = log1p(as.matrix(SummarizedExperiment::assays(ACTIONet.out$annotations[[cl.idx]]$DE.profile)$significance))
	rows = unique(as.numeric(apply(X, 2, function(x) order(x, decreasing=T)[1:top.genes])))
	markers = rownames(X)[rows]

	plot.annotations.selected.genes(ACTIONet.out, sce = sce, annotation.name = annotation.name, markers = markers, type = type, seriation.method = seriation.method, sort.rows = F)
}

plot.archetype.top.genes <- function(ACTIONet.out, top.genes, type = "heatmap", seriation.method = "OLO_complete", blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH") {
    X = log1p(as.matrix(SummarizedExperiment::assays(ACTIONet.out$unification.out$DE.core)[["significance"]]))
	filtered.rows = grep(blacklist.pattern, rownames(X))
	if(length(filtered.rows) > 0)
		X = X[-filtered.rows, ]

	selected.rows = unique(as.numeric(apply(X, 2, function(x) order(x, decreasing = T)[1:top.genes])))
	markers = rownames(X)[selected.rows]
	
	plot.archetype.selected.genes(ACTIONet.out= ACTIONet.out, selected.genes = markers, type = type, seriation.method = seriation.method, sort.rows = F)
}

plot.archetype.selected.genes <- function(ACTIONet.out, selected.genes, type = "heatmap", seriation.method = "OLO", sort.rows = T) {
	selected.rows = match(selected.genes, rownames(ACTIONet.out$unification.out$cellstates.core))
	
	Avg.profile = ACTIONet.out$unification.out$cellstates.core[selected.rows, ]

	if(sort.rows) {
		set.seed(0)
		CC = cor(as.matrix(Matrix::t(Avg.profile)))
		CC[is.na(CC)] = 0
		D = as.dist(1-CC)
		row.perm = seriation::get_order(seriation::seriate(D, seriation.method))
	}  else {
		row.perm = 1:length(selected.genes)
	}
	set.seed(0)
	CC = cor(as.matrix(Avg.profile))
	CC[is.na(CC)] = 0
	D = as.dist(1-CC)
	col.perm = seriation::get_order(seriation::seriate(D, seriation.method))
	
	gradPal = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
	
	Enrichment = Matrix::t(scale(Matrix::t(Avg.profile[row.perm, col.perm])))

	if(type == "dotplot") {
		library(corrplot)
		corrplot::corrplot(Enrichment, is.corr = F, tl.col = "black", tl.cex = 0.7, cl.pos = "n")
	} else if(type == "heatmap") {
		library(ComplexHeatmap)
		ComplexHeatmap::Heatmap(Enrichment, name = "z-score", cluster_rows = F, cluster_columns = F, col = gradPal, row_title = "", column_title = "Archeypes", column_names_gp = gpar(fontsize = 8, fontface="bold"), row_names_gp = gpar(fontsize = 8, fontface="bold"), column_title_gp = gpar(fontsize = 10, fontface="bold"), row_names_side = "left", rect_gp = gpar(col = "black"))		
	}	
}


plot.ACTIONet.archetype.footprint <- function(ACTIONet.out, node.size = 1, CPal = "magma", title = "", arch.labels = NULL) {
	Ht = Matrix::t(ACTIONet.out$unification.out$H.core)
	
	node.size = node.size * 0.3
	ACTIONet = ACTIONet.out$ACTIONet
    coors = cbind(V(ACTIONet)$x, V(ACTIONet)$y)

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
		h = Ht[, i]

		hs = sort(h, decreasing = T)
		nnz = round( (sum(hs^2)^2) / (sum(hs^4)) )
		
		h = (h / hs[nnz])^3		
		
	    vCol = scales::col_bin(Pal_grad, domain = NULL, bins = 10)(h)
        vCol = scales::alpha(vCol, 0.05 + 0.95*h/max(h))
	    idx = order(h, decreasing = F)
		plot(coors[idx, 1], coors[idx, 2], bg = vCol[idx], col = vCol[idx], cex = node.size, pch = 21, axes = FALSE, xlab = "", ylab = "", main = arch.labels[[i]])
	    
	})
}

