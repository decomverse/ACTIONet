
#' Main ACTIONet plotting functions
#'
#' @param ace ACTIONet output object
#' @param labels Annotation of interest (clusters, celltypes, etc.) to be projected on the ACTIONet plot
#' @param transparency.attr Additional continuous attribute to project onto the transparency of nodes
#' @param trans.z.threshold, trans.fact Control the effect of transparency mapping
#' @param node_size Size of nodes in the ACTIONet plot
#' @param CPal Color palette (named vector or a name for a given known palette)
#' @param add.text Whether or not to add labels on the top of ACTIONet plot
#' @param suppress.legend Whether to suppress legend or include one
#' @param legend.pos If suppress.legend == F, where should we put the legend
#' @param title Main title of the plot
#' @param border.contrast.factor How much the node and its border should contrast
#' @param add.states Whether or not to include interpolated cell state positions in the plot
#' @param coordinate_slot Entry in colMaps(ace) containing the plot coordinates (default:'ACTIONet2D')
#'
#' @return Visualized ACTIONet
#'
#' @examples
#' ace = run.ACTIONet(sce)
#' plot.ACTIONet(ace, ace$assigned_archetype, transparency.attr = ace$node_centrality)
#' @export
.plot.ACTIONet2D <- function(ace, labels = NULL, transparency.attr = NULL, trans.z.threshold = -0.5, 
    trans.fact = 1.5, node_size = 0.1, CPal = CPal20, add.text = TRUE, suppress.legend = TRUE, 
    legend.pos = "bottomright", title = "", border.contrast.factor = 0.1, add.backbone = FALSE, 
    arch.size.factor = 3, coordinate_slot = "ACTIONet2D") {
    
    text.halo.width = 0.1
    label_size = 0.8
    
    node_size = node_size * 0.25
    
    if (class(ace) == "ACTIONetExperiment") {
        labels = .preprocess_annotation_labels(labels, ace)
        if (is.character(coordinate_slot)) {
            coors = as.matrix(colMaps(ace)[[coordinate_slot]])
            coor.mu = apply(coors, 2, mean)
            coor.sigma = apply(coors, 2, sd)
            coors = scale(coors)
        } else {
            coors = as.matrix(coordinate_slot)
            coor.mu = apply(coors, 2, mean)
            coor.sigma = apply(coors, 2, sd)
            coors = scale(coors)
        }
    } else {
        if (is.matrix(ace) | is.sparseMatrix(ace)) {
            coors = as.matrix(ace)
            coor.mu = apply(coors, 2, mean)
            coor.sigma = apply(coors, 2, sd)
            coors = scale(coors)
            labels = .preprocess_annotation_labels(labels)
        } else {
            err = sprintf("Unknown type for object 'ace'.\n")
            stop(err)
        }
    }
    
    if (is.null(labels) || length(unique(labels)) == 1) {
        if (class(ace) == "ACTIONetExperiment") {
            vCol = grDevices::rgb(colMaps(ace)$denovo_color)
        } else {
            vCol = rep("tomato", nrow(coors))
        }
        Annot = NULL
    } else {
        Annot = names(labels)[match(sort(unique(labels)), labels)]
        if (length(CPal) > 1) {
            if (length(CPal) < length(Annot)) {
                if (length(Annot) <= 20) {
                  CPal = CPal20
                } else {
                  CPal = CPal88
                }
            }
            if (is.null(names(CPal))) {
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
        z = scale(transparency.attr)  # (transparency.attr - median(transparency.attr))/mad(transparency.attr)
        beta = 1/(1 + exp(-trans.fact * (z - trans.z.threshold)))
        beta[z > trans.z.threshold] = 1
        beta = beta^trans.fact
        
        vCol.border = scales::alpha(colorspace::darken(vCol, border.contrast.factor), 
            beta)
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
    x.min = x.min - (x.max - x.min)/20
    x.max = x.max + (x.max - x.min)/20
    y.min = y.min - (y.max - y.min)/20
    y.max = y.max + (y.max - y.min)/20
    XL = c(x.min, x.max)
    YL = c(y.min, y.max)
    
    
    rand.perm = sample(nrow(coors))
    graphics::plot(coors[rand.perm, c(1, 2)], pch = 21, cex = node_size, bg = vCol[rand.perm], 
        col = vCol.border[rand.perm], axes = FALSE, xlab = "", ylab = "", main = title, 
        xlim = XL, ylim = YL)
    
    if (add.text == TRUE & (!is.null(Annot))) {
        graphics::par(xpd = TRUE, mar = graphics::par()$mar * c(1.1, 1.1, 1.1, 1.1))
        
        centroids = Matrix::t(sapply(Annot, function(l) {
            idx = which(names(labels) == l)
            if (length(idx) == 1) {
                return(as.numeric(coors[idx, ]))
            }
            
            sub.coors = coors[idx, ]
            anchor.coor = as.numeric(apply(sub.coors, 2, function(x) mean(x, trim = 0.8)))
            
            return(anchor.coor)
        }))
        
        layout.labels(x = centroids[, 1], y = centroids[, 2], labels = Annot, col = colorspace::darken(Pal, 
            0.5), bg = "#eeeeee", r = text.halo.width, cex = label_size)
    }
    
    if ((suppress.legend == FALSE) & !is.null(Annot)) {
        xmin <- graphics::par("usr")[1]
        xmax <- graphics::par("usr")[2]
        ymin <- graphics::par("usr")[3]
        ymax <- graphics::par("usr")[4]
        
        lgd <- graphics::legend(x = mean(c(xmin, xmax)), y = mean(c(ymin, ymax)), 
            legend = Annot, fill = Pal, cex = 0.5, bty = "n", plot = FALSE)
        
        graphics::par(xpd = TRUE, mai = c(0, 0, 0, lgd$rect$w))
        
        graphics::legend(x = xmax, y = ymin + lgd$rect$h, legend = Annot, fill = Pal, 
            cex = 0.5, bty = "n", plot = TRUE)
        
    }
}
