# From: SnapATAC
CPal88 = c("#E31A1C", "#FFD700", "#771122", "#777711", "#1F78B4", "#68228B", "#AAAA44",
    "#60CC52", "#771155", "#DDDD77", "#774411", "#AA7744", "#AA4455", "#117744",
    "#000080", "#44AA77", "#AA4488", "#DDAA77", "#D9D9D9", "#BC80BD", "#FFED6F",
    "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17",
    "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
    "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
    "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928", "#FBB4AE", "#B3CDE3",
    "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2",
    "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC",
    "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFFF33", "#A65628",
    "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
    "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
    "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5")

# From: Scanpy
CPal20 = c("#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc", "#8c564b", "#e377c2",
    "#b5bd61", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5", "#ad494a", "#8c6d31")

layout.labels <- function(
  x,
  y,
  labels,
  col = "white",
  bg = "black",
  r = 0.1,
  cex = 1,
  ...
) {

    lay <- wordcloud::wordlayout(x, y, words = labels, cex = 1.25 * cex, ...)

    x = lay[, 1] + 0.5 * lay[, 3]
    y = lay[, 2] + 0.5 * lay[, 4]

    theta = seq(0, 2 * pi, length.out = 50)
    xy <- grDevices::xy.coords(x, y)
    xo <- r * graphics::strwidth("A")
    yo <- r * graphics::strheight("A")

    for (i in theta) {
        graphics::text(
          x = xy$x + cos(i) * xo,
          y = xy$y + sin(i) * yo,
          labels = labels,
          col = bg,
          cex = cex,
          ...
        )
    }
    graphics::text(
      x = xy$x,
      y = xy$y,
      labels = labels,
      col = col,
      cex = cex,
      ...
    )
}

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
plot.ACTIONet <- function(
  ace,
  labels = NULL,
  transparency.attr = NULL,
  trans.z.threshold = -0.5,
  trans.fact = 1.5,
  node_size = 0.1,
  CPal = CPal20,
  add.text = TRUE,
  suppress.legend = TRUE,
  legend.pos = "bottomright",
  title = "",
  border.contrast.factor = 0.1,
  add.backbone = FALSE,
  arch.size.factor = 3,
  coordinate_slot = "ACTIONet2D"
) {

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
    x.min = x.min - (x.max - x.min)/20
    x.max = x.max + (x.max - x.min)/20
    y.min = y.min - (y.max - y.min)/20
    y.max = y.max + (y.max - y.min)/20
    XL = c(x.min, x.max)
    YL = c(y.min, y.max)


    rand.perm = sample(nrow(coors))
    graphics::plot(
      coors[rand.perm, c(1, 2)],
      pch = 21,
      cex = node_size,
      bg = vCol[rand.perm],
      col = vCol.border[rand.perm],
      axes = FALSE,
      xlab = "",
      ylab = "",
      main = title,
      xlim = XL,
      ylim = YL
    )

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

        layout.labels(
          x = centroids[, 1],
          y = centroids[, 2],
          labels = Annot,
          col = colorspace::darken(Pal, 0.5),
          bg = "#eeeeee",
          r = text.halo.width,
          cex = label_size
        )
    }

    if ((suppress.legend == FALSE) & !is.null(Annot)) {
        xmin <- graphics::par("usr")[1]
        xmax <- graphics::par("usr")[2]
        ymin <- graphics::par("usr")[3]
        ymax <- graphics::par("usr")[4]

        lgd <- graphics::legend(
          x = mean(c(xmin, xmax)),
          y = mean(c(ymin, ymax)),
          legend = Annot,
          fill = Pal,
          cex = 0.5,
          bty = "n",
          plot = FALSE
        )

        graphics::par(xpd = TRUE, mai = c(0, 0, 0, lgd$rect$w))

        graphics::legend(
          x = xmax,
          y = ymin + lgd$rect$h,
          legend = Annot,
          fill = Pal, cex = 0.5,
          bty = "n",
          plot = TRUE
        )

    }
}

#' Main ACTIONet 3D plotting functions
#'
#' @param ace ACTIONet output object
#' @param labels Annotation of interest (clusters, celltypes, etc.) to be projected on the ACTIONet plot
#' @param transparency.attr Additional continuous attribute to project onto the transparency of nodes
#' @param trans.z.threshold, trans.fact Control the effect of transparency mapping
#' @param node_size Size of nodes in the ACTIONet plot
#' @param CPal Color palette (named vector or a name for a given known palette)
#' @param title Main title of the plot
#' @param coordinate_slot Entry in colMaps(ace) containing the plot coordinates (default:'ACTIONet2D')
#'
#' @return Visualized ACTIONet
#'
#' @examples
#' ace = run.ACTIONet(sce)
#' plot.ACTIONet.3D(ace, ace$assigned_archetype, transparency.attr = ace$node_centrality)
#' @export
plot.ACTIONet.3D <- function(
  ace,
  labels = NULL,
  transparency.attr = NULL,
  trans.z.threshold = -1,
  trans.fact = 1,
  node_size = 1,
  CPal = CPal20,
  coordinate_slot = "ACTIONet3D"
) {

    nV = length(ncol(ace))

    node_size = node_size * 0.2

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

    if (is.null(labels)) {
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

        vCol.border = scales::alpha(colorspace::darken(vCol, 0.5), beta)
        vCol = scales::alpha(vCol, beta)
    } else {
        vCol.border = colorspace::darken(vCol, 0.5)
    }

    threejs::scatterplot3js(
      x = coors[, 1],
      y = coors[, 2],
      z = coors[, 3],
      axis.scales = FALSE,
      size = node_size,
      axis = FALSE,
      grid = FALSE,
      color = as.character(vCol),
      stroke = as.character(vCol.border),
      bg = "black"
    )
}

#' Plots heatmap of the top-ranked features of an enrichment table
#'
#' @param feat_scores An arbitrary enrichment table with columns being archetypes and rows being specificity of features
#' @param top_features Number of features to return
#' @param reorder_columns Whether to optimally re-order columns of the enrichment table
#'
#' @return Enrichment heatmap
#'
#' @examples
#' feat_scores = as.matrix(rowMaps(ace)[['unified_feature_specificity']])
#' plot.top.k.features(feat_scores, 3)
#' @export
plot.top.k.features <- function(
  feat_scores,
  top_features = 3,
  normalize = TRUE,
  reorder_columns = TRUE,
  row.title = "Archetypes",
  column.title = "Genes",
  rowPal = "black"
) {

    W = select.top.k.features(
      feat_scores = feat_scores,
      top_features = top_features,
      normalize = normalize
    )

    gradPal = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"))))(100)

    M = apply(W, 1, max)

    Z = W
    ht = ComplexHeatmap::Heatmap(
      matrix = Z,
      name = "Expression (scaled)",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      col = gradPal,
      row_title = row.title,
      column_title = column.title,
      column_names_gp = grid::gpar(fontsize = 8, fontface = "bold"),
      row_names_gp = grid::gpar(fontsize = 8, fontface = "bold", col = rowPal),
      column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
      row_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
      row_names_side = "left",
      rect_gp = grid::gpar(col = "black"),
      row_names_max_width = grid::unit(100, "cm"),
      column_names_max_height = grid::unit(100, "cm")
    )

    return(ht)
}


#' Plots projection of features onto the ACTIONet plot
#'
#' @param ace ACTIONet output
#' @param feat_scores An arbitrary enrichment table with columns being archetypes and rows being specificity of features
#' @param top_features Number of features to return
#' @param reorder_columns Whether to optimally re-order columns of the enrichment table
#' @param CPal Color palette to use
#' @param title Main title of the plot
#'
#' @return Featur view
#'
#' @examples
#' feat_scores = as.matrix(rowMaps(ace)[['unified_feature_specificity']])
#' plot.ACTIONet.feature.view(ace, feat_scores, 5)
#' @export
plot.ACTIONet.feature.view <- function(
  ace,
  feat_scores,
  top_features = 5,
  CPal = NULL,
  title = "Feature view",
  label_size = 1,
  renormalize = FALSE,
  footprint_slot = "H_unified"
) {

    M = as(colMaps(ace)[[footprint_slot]], "sparseMatrix")

    if (ncol(feat_scores) != ncol(colMaps(ace)[["H_unified"]])) {
        feat_scores = Matrix::t(feat_scores)
    }

    if (max(feat_scores) > 50)
        feat_scores = log1p(feat_scores)

    X = t(select.top.k.features(
      feat_scores = feat_scores,
      top_features = top_features,
      normalize = renormalize,
      reorder_columns = FALSE
    ))
    selected.features = colnames(X)

    X = exp(scale(X))

    core.coors = Matrix::t(metadata(ace)$backbone$coordinates)
    cs = fastColSums(X)
    cs[cs == 0] = 1
    X = scale(X, center = FALSE, scale = cs)

    feature.coors = Matrix::t(core.coors %*% X)

    if (is.null(CPal)) {
        core.Pal = grDevices::rgb(S4Vectors::metadata(ace)$backbone$colors)
    } else {
        if (length(CPal) == 1) {
            core.Pal = ggpubr::get_palette(CPal, length(unique(ace$archetype.assignment)))
        } else {
            core.Pal = CPal[1:length(unique(ace$archetype.assignment))]
        }
    }
    core.Lab = grDevices::convertColor(
      color = Matrix::t(grDevices::col2rgb(core.Pal)/256),
      from = "sRGB",
      to = "Lab"
    )

    feature.color.Lab = Matrix::t(X) %*% core.Lab
    feature.colors = grDevices::rgb(grDevices::convertColor(
      color = feature.color.Lab,
      from = "Lab",
      to = "sRGB"
    ))
    names(feature.colors) = selected.features

    x = feature.coors[, 1]
    y = feature.coors[, 2]
    x.min = min(x)
    x.max = max(x)
    y.min = min(y)
    y.max = max(y)
    x.min = x.min - (x.max - x.min)/4
    x.max = x.max + (x.max - x.min)/4
    y.min = y.min - (y.max - y.min)/4
    y.max = y.max + (y.max - y.min)/4
    XL = c(x.min, x.max)
    YL = c(y.min, y.max)

    graphics::plot(
      x = x,
      y = y,
      type = "n",
      col = feature.colors,
      axes = FALSE,
      xlab = "",
      ylab = "",
      main = title,
      xlim = XL,
      ylim = YL
    )

    words = selected.features
    lay <- wordcloud::wordlayout(x, y, words, label_size)

    for (i in 1:length(x)) {
        xl <- lay[i, 1]
        yl <- lay[i, 2]
        w <- lay[i, 3]
        h <- lay[i, 4]
        if (x[i] < xl || x[i] > xl + w || y[i] < yl || y[i] > yl + h) {

            graphics::points(
              x = x[i],
              y = y[i],
              pch = 16,
              col = colorspace::darken(feature.colors[[i]], 0.6),
              cex = 0.75 * label_size
            )

            nx <- xl + 0.5 * w
            ny <- yl + 0.5 * h

            graphics::lines(
              x = c(x[i], nx),
              y = c(y[i], ny),
              col = colorspace::darken(feature.colors[[i]], 0.5)
            )
        }
    }

    loc.x = lay[, 1] + 0.5 * lay[, 3]
    loc.y = lay[, 2] + 0.5 * lay[, 4]

    graphics::text(
      x = loc.x,
      y = loc.y,
      labels = words,
      col = feature.colors,
      cex = label_size
    )
}


#' Plots projection of genes onto the ACTIONet plot
#'
#' @param ace ACTIONet output
#' @param top_genes Number of genes to return
#' @param blacklist_pattern List of genes to filter-out
#' @param reorder_columns Whether to optimally re-order columns of the enrichment table
#' @param CPal Color palette to use
#' @param title Main title of the plot
#'
#' @return Featur view
#'
#' @examples
#' plot.ACTIONet.gene.view(ace, 5)
#' @export
plot.ACTIONet.gene.view <- function(
  ace,
  top_genes = 5,
  CPal = NULL,
  blacklist_pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|MALAT1|B2M|GAPDH",
  title = "",
  label_size = 0.8,
  renormalize = FALSE
) {

    feat_scores = as.matrix(rowMaps(ace)[["unified_feature_specificity"]])
    filtered.rows = grep(blacklist_pattern, rownames(feat_scores))
    if (length(filtered.rows) > 0)
        feat_scores = feat_scores[-filtered.rows, ]

    plot.ACTIONet.feature.view(
      ace = ace,
      feat_scores = feat_scores,
      title = "Gene view",
      renormalize = renormalize,
      top_features = top_genes
    )
}


#' Main ACTIONet interactive plotting function
#'
#' @param ace ACTIONet output object
#' @param labels Annotation of interest (clusters, celltypes, etc.) to be projected on the ACTIONet plot
#' @param transparency.attr Additional continuous attribute to project onto the transparency of nodes
#' @param trans.z.threshold, trans.fact Control the effect of transparency mapping
#' @param node_size Size of nodes in the ACTIONet plot
#' @param CPal Color palette (named vector or a name for a given known palette)
#' @param enrichment.table To project the top-ranked features interactively.
#' @param top_features Number of features to show per cell
#' @param blacklist_pattern List of genes to filter-out
#' @param title Main title of the plot
#' @param threeD Whether to show the plot in 3D
#'
#' @return Visualized ACTIONet
#'
#' @examples
#' ace = run.ACTIONet(sce)
#' plot.ACTIONet.interactive(ace, ace$assigned_archetype)
#' @import plotly
#' @export
plot.ACTIONet.interactive <- function(
  ace,
  labels = NULL,
  transparency.attr = NULL,
  trans.z.threshold = -1,
  trans.fact = 1,
  node_size = 1,
  CPal = CPal20,
  enrichment.table = NULL,
  top_features = 7,
  Alt_Text = NULL,
  blacklist_pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|MALAT1|B2M|GAPDH",
  threeD = FALSE,
  title = "ACTIONet",
  coordinate_slot = "ACTIONet2D"
) {

    nV = ncol(ace)
    node_size = node_size * 3
    if (coordinate_slot == "ACTIONet2D" & threeD == TRUE)
        coordinate_slot = "ACTIONet3D"

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

    if (is.null(labels)) {
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

        vCol.border = scales::alpha(colorspace::darken(vCol, 0.5), beta)
        vCol = scales::alpha(vCol, beta)
    } else {
        vCol.border = colorspace::darken(vCol, 0.5)
    }

    if (!is.null(enrichment.table)) {
        if (ncol(enrichment.table) == nV) {
            cell.scores = Matrix::t(enrichment.table)
        } else if ((nrow(enrichment.table) != nV)) {
            H = colMaps(ace)[["H_unified"]]
            if ((nrow(enrichment.table) == nrow(H)) | (ncol(enrichment.table) ==
                nrow(H))) {
                cell.scores = map.cell.scores.from.archetype.enrichment(ace, enrichment.table)
            } else {
                cell.scores = NULL
            }
        } else {
            cell.scores = enrichment.table
        }
    } else {
        temp.enrichment.table = as.matrix(rowMaps(ace)[["unified_feature_specificity"]])
        if (!is.null(row.names(temp.enrichment.table))) {
            filtered.rows = grep(blacklist_pattern, rownames(temp.enrichment.table))
            if (length(filtered.rows) > 0){
              enrichment.table = temp.enrichment.table[-filtered.rows, ]
            } else {
              enrichment.table = temp.enrichment.table
            }

            GT = apply(enrichment.table, 2, function(x) rownames(enrichment.table)[order(x,
                decreasing = TRUE)[1:min(100, nrow(enrichment.table))]])
            selected.features = sort(unique(as.character(GT)))

            W = exp(scale(Matrix::t(colMaps(ace)[["H_unified"]])))
            cs = fastColSums(W)
            W = Matrix::t(scale(W, center = FALSE, scale = cs))

            cell.scores = W %*% Matrix::t(enrichment.table[selected.features, ])
        } else {
            cell.scores = NULL
        }
    }

    if (!is.null(Alt_Text)) {
        node.annotations = Alt_Text
    } else {
        if (!is.null(cell.scores)) {
            selected.features = colnames(cell.scores)
            node.annotations = apply(cell.scores, 1, function(x) paste(selected.features[order(x,
                decreasing = TRUE)[1:top_features]], collapse = "\n"))
        } else {
            node.annotations = rep("", nV)
        }
    }
    # Setup visualization parameters
    sketch.graph = igraph::graph_from_adjacency_matrix(
      adjmatrix = colNets(ace)$ACTIONet,
      mode = "undirected",
      weighted = TRUE
    )

    if (threeD == FALSE) {
        igraph::V(sketch.graph)$x = coors[, 1]
        igraph::V(sketch.graph)$y = coors[, 2]
    } else {
        igraph::V(sketch.graph)$x3D = coors[, 1]
        igraph::V(sketch.graph)$y3D = coors[, 2]
        igraph::V(sketch.graph)$z3D = coors[, 3]
    }

    sketch.graph = igraph::delete_edges(sketch.graph, igraph::E(sketch.graph))

    node.data <- igraph::get.data.frame(sketch.graph, what = "vertices")
    edge.data <- igraph::get.data.frame(sketch.graph, what = "edges")

    Nv <- dim(node.data)[1]
    Ne <- dim(edge.data)[1]

    edge_shapes <- list()

    # Adjust parameters
    node.data$size = node_size

    axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)

    if (threeD == TRUE) {
        if (is.null(Annot)) {
            node.data$vCol = vCol
            node.data$vCol.border = vCol.border

            network <- plotly::plot_ly(
              data = node.data,
              x = ~x3D,
              y = ~y3D,
              z = ~z3D,
              opacity = 1,
              marker = list(
                color = ~vCol,
                size = ~size,
                opacity = 1,
                alpha = 1,
                line = list(
                  width = 0.1 * node_size,
                  alpha = 0.5,
                  color = ~vCol.border
                )
              ),
              text = node.annotations,
              mode = "markers",
              hoverinfo = "text",
              type = "scatter3d",
              showlegend = FALSE
            )

            p <- plotly::layout(
              p = network,
              title = title,
              shapes = edge_shapes,
              scene = list(
                xaxis = axis,
                yaxis = axis,
                zaxis = axis)
              )

        } else {
            node.data$vCol.border = vCol.border
            node.data$type = factor(names(labels), levels = Annot)

            network <- plotly::plot_ly(
              data = node.data,
              x = ~x3D,
              y = ~y3D,
              z = ~z3D,
              opacity = 1,
              color = ~type,
              colors = Pal,
              marker = list(
                size = ~size,
                opacity = 1,
                alpha = 1,
                line = list(
                  width = 0.1 * node_size,
                  alpha = 0.5,
                  color = ~vCol.border
                )
              ),
              text = node.annotations,
              mode = "markers",
              hoverinfo = "text",
              type = "scatter3d"
            )

            p <- plotly::layout(
              p = network,
              title = title,
              shapes = edge_shapes,
              scene = list(
                xaxis = axis,
                yaxis = axis,
                zaxis = axis
              ),
              showlegend = TRUE,
              legend = list(
                marker = list(
                  marker.size = 10
                )
              )
            )
        }
    } else {
        if (is.null(Annot)) {
            node.data$vCol = vCol
            node.data$vCol.border = vCol.border
            network <- plotly::plot_ly(
              data = node.data,
              x = ~x,
              y = ~y,
              marker = list(
                color = ~vCol,
                size = ~size,
                opacity = 1,
                alpha = 1,
                line = list(
                  width = 0.1 * node_size,
                  alpha = 0.5,
                  color = ~vCol.border
                )
              ),
              text = node.annotations,
              mode = "markers",
              type = "scattergl",
              hoverinfo = "text",
              showlegend = FALSE
            )

            p <- plotly::layout(
              p = network,
              title = title,
              shapes = edge_shapes,
              xaxis = axis,
              yaxis = axis
            )
        } else {
            node.data$vCol.border = vCol.border
            node.data$type = factor(names(labels), levels = Annot)

            network <- plotly::plot_ly(
              data = node.data,
              x = ~x,
              y = ~y,
              color = ~type,
              colors = Pal,
              marker = list(
                size = ~size,
                line = list(
                  width = 0.1 * node_size,
                  color = ~vCol.border
                )
              ),
              text = node.annotations,
              mode = "markers",
              type = "scattergl",
              hoverinfo = "text"
            )

            p <- plotly::layout(
              p = network,
              title = title,
              shapes = edge_shapes,
              xaxis = axis,
              yaxis = axis,
              showlegend = TRUE,
              legend = list(
                marker = list(
                  marker.size = 10
                )
              )
            )
        }
    }

    p
}


#' Plot gene expression violin plot
#'
#' @param ace ACTIONet output object
#' @param labels Annotation of interest (clusters, celltypes, etc.) to be projected on the ACTIONet plot
#' @param gene_name Name of the gene to plot
#' @param CPal Color palette (named vector or a name for a given known palette)
#'
#' @return Visualized ACTIONet
#'
#' @examples
#' plot.individual.gene(ace, ace$assigned_archetype, 'CD14')
#' @export
plot.individual.gene <- function(
  ace,
  labels,
  gene_name,
  features_use = NULL,
  assay_name = "logcounts",
  CPal = CPal20
) {

    clusters = .preprocess_annotation_labels(ace, labels)
    features_use = .preprocess_annotation_features(ace, features_use)

    Labels = names(clusters)
    Annot = sort(unique(Labels))
    Annot = Annot[order(clusters[match(Annot, Labels)], decreasing = FALSE)]
    Labels = factor(Labels, levels = Annot)

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

    if (!(gene_name %in% features_use)) {
        err = sprintf("Gene %s not found\n", gene_name)
        stop(err)
    }

    x = SummarizedExperiment::assays(ace)[[assay_name]][gene_name, ]
    if (sum(x) == 0) {
        err = sprintf("Gene must have non-zero expression.\n")
        stop(err)
    }


    df = data.frame(Annotation = Labels, Expression = x)
    gp = ggpubr::ggviolin(
      data = df,
      x = "Annotation",
      y = "Expression",
      fill = "Annotation",
      palette = Pal,
      add = "boxplot",
      add.params = list(fill = "white")
    )
    print(gp)
}


#' Projects a given continuous score on the ACTIONet plot
#'
#' @param ace ACTIONet output object
#' @param x Score vector
#' @param transparency.attr Additional continuous attribute to project onto the transparency of nodes
#' @param trans.z.threshold, trans.fact Control the effect of transparency mapping
#' @param node_size Size of nodes in the ACTIONet plot
#' @param CPal Color palette (named vector or a name for a given known palette)
#' @param coordinate_slot Entry in colMaps(ace) containing the plot coordinates (default:'ACTIONet2D')
#' @param alpha_val Between [0, 1]. If it is greater than 0, smoothing of scores would be performed
#'
#' @return Visualized ACTIONet with projected scores
#'
#' @examples
#' ace = run.ACTIONet(sce)
#' x = logcounts(ace)['CD14', ]
#' plot.ACTIONet.gradient(ace, x, transparency.attr = ace$node_centrality)
#' @export
plot.ACTIONet.gradient <- function(
  ace,
  x,
  transparency.attr = NULL,
  trans.z.threshold = -0.5,
  trans.fact = 3,
  node_size = 0.1,
  CPal = "magma",
  title = "",
  alpha_val = 0.85,
  nonparameteric = FALSE,
  coordinate_slot = "ACTIONet2D"
) {

    node_size = node_size * 0.3

    if (class(ace) == "ACTIONetExperiment") {
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
        } else {
            err = sprintf("Unknown type for object 'ace'.\n")
            stop(err)
        }
    }


    NA_col = "#eeeeee"

    ## Create color gradient generator
    if (CPal %in% c("greys", "inferno", "magma", "viridis", "BlGrRd", "RdYlBu", "Spectral")) {

        Pal_grad = switch(CPal,
          greys = grDevices::gray.colors(100),
          inferno = viridis::inferno(500, alpha = 0.8),
          magma = viridis::magma(500, alpha = 0.8),
          viridis = viridis::viridis(500, alpha = 0.8),
          BlGrRd = grDevices::colorRampPalette(c("blue", "grey", "red"))(500),
          Spectral = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "Spectral"))))(100),
          RdYlBu = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"))))(100)
        )
    } else {
        Pal_grad = grDevices::colorRampPalette(c(NA_col, CPal))(500)
    }

    ## Scale/prune scorees, if needed
    x[x < 0] = 0
    if (max(x) > 50)
        x = log1p(x)

    if (alpha_val > 0) {
        x = as.numeric(compute_network_diffusion(
          G = colNets(ace)$ACTIONet,
          X0 = as(as.matrix(x), "sparseMatrix")
        ))
    }

    if (nonparameteric == TRUE) {

        vCol = (scales::col_bin(
          palette = Pal_grad,
          domain = NULL,
          na.color = NA_col,
          bins = 7
        ))(rank(x))

    } else {

        vCol = (scales::col_bin(
          palette = Pal_grad,
          domain = NULL,
          na.color = NA_col,
          bins = 7
        ))(x)

    }

    if (!is.null(transparency.attr)) {
        z = scale(transparency.attr)  # (transparency.attr - median(transparency.attr))/mad(transparency.attr)
        beta = 1/(1 + exp(-trans.fact * (z - trans.z.threshold)))
        beta[z > trans.z.threshold] = 1
        beta = beta^trans.fact

        vCol = scales::alpha(vCol, beta)
        vCol.border = scales::alpha(colorspace::darken(vCol, 0.1), beta)
    } else {
        vCol.border = colorspace::darken(vCol, 0.1)
    }

    idx = order(x, decreasing = FALSE)

    graphics::plot(
      x = coors[idx, 1],
      y = coors[idx, 2],
      bg = vCol[idx],
      col = vCol.border[idx],
      cex = node_size,
      pch = 21,
      axes = FALSE,
      xlab = "",
      ylab = "",
      main = title
    )

}


#' Projects a set of markers on the ACTIONet
#' It also optionally imputes the markers.
#'
#' @param ace ACTIONet output object
#' @param markers Set of row features (e.g. genes) to visualize .
#' @param transparency.attr Additional continuous attribute to project onto the transparency of nodes
#' @param trans.z.threshold, trans.fact Control the effect of transparency mapping
#' @param node_size Size of nodes in the ACTIONet plot
#' @param CPal Color palette (named vector or a name for a given known palette)
#' @param coordinate_slot Entry in colMaps(ace) containing the plot coordinates (default:'ACTIONet2D')
#' @param alpha_val Between [0, 1]. If it is greater than 0, smoothing of scores would be performed
#'
#' @return Visualized ACTIONet with projected scores
#'
#' @examples
#' ace = run.ACTIONet(sce)
#' visualize.markers(ace, markers = c('CD14', 'CD19', 'CD3G'), transparency.attr = ace$node_centrality)
visualize.markers <- function(
  ace,
  markers,
  features_use = NULL,
  assay_name = "logcounts",
  transparency.attr = NULL,
  trans.z.threshold = -0.5,
  trans.fact = 3,
  node_size = 1,
  CPal = "magma",
  alpha_val = 0.9,
  export_path = NA
) {

    features_use = .preprocess_annotation_features(ace, features_use = features_use)
    markers_all = sort(unique(unlist(markers)))
    marker_set = intersect(markers_all, features_use)

    if (length(marker_set) == 0) {
        err = sprintf("No given markers found in feature set.\n")
        stop(err, call. = FALSE)
    }

    if (length(marker_set) == 1)
        alpha_val = 0

    if (alpha_val > 0) {

        expression_profile = impute.genes.using.ACTIONet(
          ace = ace,
          genes = marker_set,
          features_use = features_use,
          alpha_val = alpha_val
        )

    } else {
        expression_profile = assays(ace)[[assay_name]][match(marker_set, features_use),
            , drop = FALSE]
        expression_profile = Matrix::t(expression_profile)
        colnames(expression_profile) = marker_set
    }

    print(sprintf("Markers Visualized: %s", paste0(marker_set, collapse = ", ")))
    markers_missing = setdiff(markers_all, marker_set)
    if (length(markers_missing) > 0)
        print(sprintf("Markers Missing: %s", paste0(markers_missing, collapse = ", ")))

    for (i in 1:ncol(expression_profile)) {
        feat_name = colnames(expression_profile)[i]
        x = expression_profile[, i]

        nnz = round(sum(x^2)^2/sum(x^4))
        x.threshold = sort(x, decreasing = TRUE)[nnz]
        x[x < x.threshold] = 0
        x = x/max(x)

        plot.ACTIONet.gradient(
          ace = ace,
          x = x,
          transparency.attr = transparency.attr,
          trans.z.threshold = trans.z.threshold,
          trans.fact = trans.fact,
          node_size = node_size,
          CPal = CPal,
          title = feat_name,
          alpha_val = 0
        )

        if (!is.na(export_path)) {
            fname = sprintf("%s/%s.pdf", export_path, feat_name)
            dir.create(dirname(fname), showWarnings = FALSE, recursive = TRUE)
            pdf(fname)
            par(mar = c(0, 0, 1, 0))
            plot.ACTIONet.gradient(
              ace = ace,
              x = x,
              transparency.attr = transparency.attr,
              trans.z.threshold = trans.z.threshold,
              trans.fact = trans.fact,
              node_size = node_size,
              CPal = CPal,
              title = feat_name,
              alpha_val = 0
            )
            dev.off()
        }
    }
}


select.top.k.genes <- function(
  ace,
  top_genes = 5,
  CPal = NULL,
  blacklist_pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH",
  top_features = 3,
  normalize = FALSE,
  reorder_columns = FALSE,
  slot_name = "unified_feature_specificity"
) {

    feat_scores = as.matrix(rowMaps(ace)[[slot_name]])
    filtered.rows = grep(blacklist_pattern, rownames(feat_scores))
    if (length(filtered.rows) > 0)
        feat_scores = feat_scores[-filtered.rows, ]

    tbl = select.top.k.features(
      feat_scores = feat_scores,
      top_features = top_features,
      normalize = normalize,
      reorder_columns = reorder_columns
    )

    return(tbl)
}


plot.top.k.genes <- function(
  ace,
  top_genes = 5,
  CPal = NULL,
  blacklist_pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH",
  top_features = 3,
  normalize = FALSE,
  reorder_columns = TRUE,
  row.title = "Archetypes",
  column.title = "Genes",
  rowPal = "black",
  slot_name = "unified_feature_specificity"
) {

    feat_scores = as.matrix(rowMaps(ace)[[slot_name]])
    filtered.rows = grep(blacklist_pattern, rownames(feat_scores))
    if (length(filtered.rows) > 0)
        feat_scores = feat_scores[-filtered.rows, ]

    ht = plot.top.k.features(
      feat_scores = feat_scores,
      top_features = top_features,
      normalize = normalize,
      reorder_columns = reorder_columns,
      row.title = row.title,
      column.title = column.title,
      rowPal = rowPal
    )

    return(ht)
}


plot.archetype.selected.genes <- function(
  ace,
  genes,
  CPal = NULL,
  blacklist_pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH",
  top_features = 3,
  normalize = FALSE,
  reorder_columns = TRUE,
  row.title = "Archetypes",
  column.title = "Genes",
  rowPal = "black",
  slot_name = "unified_feature_specificity"
) {

    feat_scores = as.matrix(rowMaps(ace)[["unified_feature_specificity"]])
    filtered.rows = match(intersect(rownames(ace), genes), rownames(ace))

    if (length(filtered.rows) > 0)
        feat_scores = feat_scores[-filtered.rows, ]

    ht = plot.top.k.features(
      feat_scores = feat_scores,
      top_features = top_features,
      normalize = normalize,
      reorder_columns = reorder_columns,
      row.title = row.title,
      column.title = column.title,
      rowPal = rowPal
    )

    return(ht)

}


plot.ACTIONet.archetype.footprint <- function(
  ace,
  node_size = 0.1,
  CPal = "magma",
  title = "",
  arch.labels = NULL,
  coordinate_slot = "ACTIONet2D",
  alpha_val = 0.9
) {

    Ht = colMaps(ace)[["H_unified"]]
    cs = fastColSums(Ht)
    cs[cs == 0] = 1

    U = as(scale(Ht, center = FALSE, scale = cs), "dgTMatrix")
    U.pr = compute_network_diffusion(
      G = colNets(ace)$ACTIONet,
      X0 = U,
      alpha = alpha_val
    )

    node_size = node_size * 0.3
    coors = scale(colMaps(ace)[[coordinate_slot]])

    if (CPal %in% c("inferno", "magma", "viridis", "BlGrRd", "RdYlBu", "Spectral")) {

        Pal_grad = switch(CPal,
          inferno = viridis::inferno(500, alpha = 0.8),
          magma = viridis::magma(500, alpha = 0.8),
          viridis = viridis::viridis(500, alpha = 0.8),
          BlGrRd = grDevices::colorRampPalette(c("blue", "grey", "red"))(500),
          Spectral = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "Spectral"))))(100),
          RdYlBu = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"))))(100))

    } else {
        NA_col = "#cccccc"
        Pal_grad = grDevices::colorRampPalette(c(NA_col, CPal))(500)
    }

    k1 = k2 = round(sqrt(NCOL(Ht)))
    if (k1 * k2 < NCOL(Ht)) {
        k2 = k2 + 1
    }

    if (is.null(arch.labels))
        arch.labels = sapply(1:NCOL(Ht), function(i) sprintf("Archetype %d", i))

    # par(mfrow = c(k1, k2), mar = c(0, 0, 1, 0))
    sapply(1:NCOL(Ht), function(i) {
        print(i)
        x = U.pr[, i]

        xs = sort(x, decreasing = TRUE)
        nnz = round((sum(xs)^2)/(sum(xs^2)))
        threshold = xs[nnz]

        x[x < threshold] = threshold
        x = log(x)

        vCol = (scales::col_bin(Pal_grad, domain = NULL, bins = 10))(x)
        vCol = scales::alpha(vCol, 0.05 + 0.95 * x/max(x))
        vCol = colorspace::lighten(vCol, 0.2)
        idx = order(x, decreasing = FALSE)

        graphics::plot(
          x = coors[idx, 1],
          y = coors[idx, 2],
          bg = vCol[idx],
          col = vCol[idx],
          cex = node_size,
          pch = 21,
          axes = FALSE,
          xlab = "",
          ylab = "",
          main = arch.labels[[i]]
        )
    })
}


#' Report the top-rank features from a given enrichment table
#'
#' @param top_features Number of features to return
#' @param reorder_columns Whether to optimally re-order columns of the enrichment table
#'
#' @return Sorted table with the selected top-ranked
#'
#' @examples
#' feat_scores = as.matrix(rowMaps(ace)[['unified_feature_specificity']])
#' enrichment.table.top = select.top.k.features(feat_scores, 3)
#' @export
select.top.k.features <- function(
  feat_scores,
  top_features = 3,
  normalize = FALSE,
  reorder_columns = TRUE
) {

    W0 = (feat_scores)
    if (normalize == TRUE)
        W0 = doubleNorm(W0)

    IDX = matrix(0, nrow = top_features, ncol = NCOL(W0))
    VV = matrix(0, nrow = top_features, ncol = NCOL(W0))
    W = (W0)

    for (i in 1:NROW(IDX)) {
        W.m = as(MWM_hungarian(W), "dgTMatrix")
        IDX[i, W.m@j + 1] = W.m@i + 1
        VV[i, W.m@j + 1] = W.m@x
        W[IDX[i, W.m@j + 1], ] = 0
    }

    if (reorder_columns == TRUE) {
        feat_scores_agg = apply(IDX, 2, function(perm) as.numeric(fastColMeans(W0[perm, ])))
        CC = cor(feat_scores_agg)
        D = stats::as.dist(1 - CC)
        cols = seriation::get_order(seriation::seriate(D, "OLO"))
        rows = as.numeric(IDX[, cols])
    } else {
        cols = 1:NCOL(W0)
        rows = unique(as.numeric(IDX))
    }

    W = feat_scores[rows, cols]

    return(W)
}


#' @export
plot.ACTIONet.backbone <- function(
  ace,
  labels = NULL,
  arch.labels = NULL,
  transparency.attr = NULL,
  trans.z.threshold = -0.5,
  trans.fact = 1.5,
  node_size = 0.1,
  CPal = CPal20,
  title = "",
  border.contrast.factor = 0.1,
  arch.size.factor = 1,
  label_size = 1,
  coordinate_slot = "ACTIONet2D"
) {

    if (!("backbone" %in% names(metadata(ace)))) {
        err = sprintf("Cannot find backbone in metadata(ace). Please run construct.backbone() first.")
        stop(err)
    }
    backbone = metadata(ace)$backbone

    node_size = node_size * 0.3

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

    if (is.null(labels)) {
        if (class(ace) == "ACTIONetExperiment") {
            vCol = grDevices::rgb(colMaps(ace)$denovo_color)
        } else {
            vCol = rep("tomato", NROW(coors))
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

        if (is.null(arch.labels)) {
            arch.annot = annotate.archetypes.using.labels(ace, labels)
            arch.labels = arch.annot$Labels
        }
    }
    if (is.null(arch.labels)) {
        arch.labels = paste("A", 1:nrow(backbone$G), sep = "")
    } else {
        arch.labels = paste("A", 1:nrow(backbone$G), "-", arch.labels, sep = "")
    }

    if (!is.null(transparency.attr)) {
        z = scale(transparency.attr)  # (transparency.attr - median(transparency.attr))/mad(transparency.attr)
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
    x.min = x.min - (x.max - x.min)/20
    x.max = x.max + (x.max - x.min)/20
    y.min = y.min - (y.max - y.min)/20
    y.max = y.max + (y.max - y.min)/20
    XL = c(x.min, x.max)
    YL = c(y.min, y.max)


    rand.perm = sample(nrow(coors))
    graphics::plot(
      coors[rand.perm, c(1, 2)],
      pch = 21,
      cex = node_size,
      bg = vCol[rand.perm],
      col = vCol.border[rand.perm],
      axes = FALSE,
      xlab = "",
      ylab = "",
      main = title,
      xlim = XL,
      ylim = YL
    )

    ## Add backbone anchors
    cell.RGB = Matrix::t(grDevices::col2rgb(vCol))/255
    cells.Lab = grDevices::convertColor(color = cell.RGB, from = "sRGB", to = "Lab")
    arch.Lab = Matrix::t(ace$archetype_footprint) %*% cells.Lab
    arch.RGB = grDevices::convertColor(color = arch.Lab, from = "Lab", to = "sRGB")
    aCol = grDevices::rgb(arch.RGB)

    w = fastColSums(colMaps(ace)$H_unified)
    w = 0.3 + 0.7 * (w - min(w))/(max(w) - min(w))
    arch.sizes = (0.25 + arch.size.factor * w)

    cell.coors = colMaps(ace)[[coordinate_slot]]
    arch.coors = backbone$coordinates[, c(1, 2)]
    arch.coors = sapply(1:2, function(i) (arch.coors[, i] - mean(cell.coors[, i]))/sd(cell.coors[,
        i]))

    text.halo.width = 0.1
    label_size = label_size * 0.6


    x = arch.coors[, 1]
    y = arch.coors[, 2] - graphics::strheight("A")
    theta = seq(0, 2 * pi, length.out = 50)
    xy <- grDevices::xy.coords(x, y)
    xo <- text.halo.width * graphics::strwidth("A")
    yo <- text.halo.width * graphics::strheight("A")
    for (i in theta) {

        graphics::text(
          x = xy$x + cos(i) * xo,
          y = xy$y + sin(i) * yo,
          labels = arch.labels,
          col = "#dddddd",
          cex = label_size
        )

    }

    graphics::text(
      x = xy$x,
      y = xy$y,
      labels = arch.labels,
      col = colorspace::darken(aCol, 0.5),
      cex = label_size
    )


    graphics::par(new = TRUE)
    graphics::plot(
      arch.coors,
      pch = 25,
      cex = arch.sizes,
      bg = aCol,
      col = colorspace::darken(aCol, 0.5),
      axes = FALSE,
      xlab = "",
      ylab = "",
      main = title,
      xlim = XL,
      ylim = YL
    )
}

#' @export
plot.ACTIONet.backbone.graph <- function(
  ace,
  labels = NULL,
  arch.labels = NULL,
  transparency.attr = NULL,
  trans.z.threshold = -0.5,
  trans.fact = 1.5,
  node_size = 0.1,
  CPal = CPal20,
  title = "",
  border.contrast.factor = 0.1,
  arch.size.factor = 1,
  label_size = 1,
  cell_lightening = 0.5,
  coordinate_slot = "ACTIONet2D",
  stretch.factor = 10
) {

    if (!("backbone" %in% names(metadata(ace)))) {
      err = sprintf("Cannot find backbone in metadata(ace). Please run construct.backbone() first.")
      stop(err)
    }
    backbone = metadata(ace)$backbone

    node_size = node_size * 0.3

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

    if (is.null(labels)) {
        if (class(ace) == "ACTIONetExperiment") {
            vCol = grDevices::rgb(colMaps(ace)$denovo_color)
        } else {
            vCol = rep("tomato", NROW(coors))
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

        if (is.null(arch.labels)) {
            arch.annot = annotate.archetypes.using.labels(ace, labels)
            arch.labels = arch.annot$Labels
        }
    }
    if (is.null(arch.labels)) {
        arch.labels = paste("A", 1:nrow(backbone$G), sep = "")
    } else {
        arch.labels = paste("A", 1:nrow(backbone$G), "-", arch.labels, sep = "")
    }

    if (!is.null(transparency.attr)) {
        z = scale(transparency.attr)  # (transparency.attr - median(transparency.attr))/mad(transparency.attr)
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
    x.min = x.min - (x.max - x.min)/20
    x.max = x.max + (x.max - x.min)/20
    y.min = y.min - (y.max - y.min)/20
    y.max = y.max + (y.max - y.min)/20
    XL = c(x.min, x.max)
    YL = c(y.min, y.max)

    rand.perm = sample(NROW(coors))
    vCol_lightend = colorspace::lighten(vCol, cell_lightening)
    vCol.border_lightend = colorspace::lighten(vCol.border, cell_lightening)
    graphics::plot(
      coors[rand.perm, c(1, 2)],
      pch = 21,
      cex = node_size,
      bg = vCol_lightend[rand.perm],
      col = vCol.border_lightend[rand.perm],
      axes = FALSE,
      xlab = "",
      ylab = "",
      main = title,
      xlim = XL,
      ylim = YL
    )

    ## Add backbone anchors
    cell.RGB = Matrix::t(grDevices::col2rgb(vCol))/255
    cells.Lab = grDevices::convertColor(color = cell.RGB, from = "sRGB", to = "Lab")
    arch.Lab = Matrix::t(ace$archetype_footprint) %*% cells.Lab
    arch.RGB = grDevices::convertColor(color = arch.Lab, from = "Lab", to = "sRGB")
    aCol = grDevices::rgb(arch.RGB)

    w = fastColSums(colMaps(ace)$H_unified)
    w = 0.3 + 0.7 * (w - min(w))/(max(w) - min(w))
    arch.sizes = (0.25 + arch.size.factor * w)

    # cell.coors =colMaps(ace)[[coordinate_slot]]
    arch.coors = backbone$coordinates
    arch.coors[, 1] = (arch.coors[, 1] - coor.mu[1])/coor.sigma[1]
    arch.coors[, 2] = (arch.coors[, 2] - coor.mu[2])/coor.sigma[2]

    graphics::par(new = TRUE)
    graphics::plot(
      arch.coors,
      pch = 25,
      cex = 0,
      bg = aCol,
      col = colorspace::darken(aCol, 0.5),
      axes = FALSE,
      xlab = "",
      ylab = "",
      main = title,
      xlim = XL,
      ylim = YL
    )

    G = construct.tspanner(backbone$G, stretch.factor = stretch.factor)
    Adj = as(G, "dgTMatrix")

    kappa = 0.25 + 0.75/(1 + exp(-2 * scale(Adj@x)))
    graphics::segments(
      x0 = arch.coors[Adj@i + 1, 1],
      y0 = arch.coors[Adj@i + 1, 2],
      x1 = arch.coors[Adj@j + 1, 1],
      y1 = arch.coors[Adj@j + 1, 2],
      col = grDevices::rgb(0, 0, 0, 0.8 * kappa),
      lwd = kappa * 3
    )

    text.halo.width = 0.1
    label_size = label_size * 0.6


    layout.labels(
      x = arch.coors[, 1],
      y = arch.coors[, 2] - graphics::strheight("A"),
      labels = arch.labels,
      col = colorspace::darken(aCol, 0.5),
      bg = "#eeeeee",
      r = text.halo.width,
      cex = label_size
    )

    graphics::par(new = TRUE)
    graphics::plot(
      arch.coors,
      pch = 25,
      cex = arch.sizes,
      bg = aCol,
      col = colorspace::darken(aCol, 0.5),
      axes = FALSE,
      xlab = "",
      ylab = "",
      main = title,
      xlim = XL,
      ylim = YL
    )

}

#' @export
plot.backbone.graph <- function(
  ace,
  arch.labels = NULL,
  arch.colors = NULL,
  node_size = 2,
  label_size = 1,
  title = "",
  stretch.factor = 2
) {

    if (!("backbone" %in% names(metadata(ace)))) {
      err = sprintf("Cannot find backbone in metadata(ace). Please run construct.backbone() first.")
      stop(err)
    }
    backbone = metadata(ace)$backbone

    G = construct.tspanner(backbone$G, stretch.factor = stretch.factor)
    Adj = as(G, "dgTMatrix")
    arch.coors = sgd2_layout_weighted_convergent(G, Matrix::t(backbone$coordinates_3D))

    if (is.null(arch.colors)) {
        arch.colors = rgb(backbone$colors)
    }

    w = fastColSums(colMaps(ace)$H_unified)
    w = 0.3 + 0.7 * (w - min(w))/(max(w) - min(w))
    arch.sizes = (0.25 + w) * node_size

    graphics::plot(
      arch.coors,
      pch = 25,
      cex = 0,
      bg = arch.colors,
      col = colorspace::darken(arch.colors, 0.5),
      axes = FALSE,
      xlab = "",
      ylab = "",
      main = title
    )

    kappa = 0.1 + 0.9/(1 + exp(-2 * scale(Adj@x)))
    graphics::segments(
      x0 = arch.coors[Adj@i + 1, 1],
      y0 = arch.coors[Adj@i + 1, 2],
      x1 = arch.coors[Adj@j + 1, 1],
      y1 = arch.coors[Adj@j + 1, 2],
      col = grDevices::rgb(0, 0, 0, 0.8 * kappa),
      lwd = kappa * 3
    )

    text.halo.width = 0.1
    label_size = label_size * 0.6

    if (is.null(arch.labels)) {
        arch.labels = paste("A", 1:nrow(arch.coors), sep = "")
    } else {
        arch.labels = paste("A", 1:nrow(arch.coors), "-", arch.labels, sep = "")
    }

    layout.labels(
      x = arch.coors[, 1],
      y = arch.coors[, 2] - graphics::strheight("A"),
      labels = arch.labels,
      col = colorspace::darken(arch.colors, 0.5),
      bg = "#eeeeee",
      r = text.halo.width,
      cex = label_size
    )

    graphics::par(new = TRUE)
    graphics::plot(
      arch.coors,
      pch = 25,
      cex = arch.sizes,
      bg = arch.colors,
      col = colorspace::darken(arch.colors, 0.5),
      axes = FALSE,
      xlab = "",
      ylab = ""
    )

}

gate.archetypes <- function(
  ace,
  i,
  j,
  H.slot = "H_unified"
) {

    H = colMaps(ace)[[H.slot]]
    hx = H[i, ]
    hy = H[j, ]

    hx.nnz = sum(hx)^2/sum(hx^2)
    hx.threshold = sort(hx, decreasing = T)[hx.nnz]

    hy.nnz = sum(hy)^2/sum(hy^2)
    hy.threshold = sort(hx, decreasing = T)[hy.nnz]


    mask = (hx > hx.threshold) | (hy > hy.threshold)
    fig <- plotly::plot_ly(x = hx[mask], y = hy[mask])

    fig <- add_trace(p = fig, type = "histogram2dcontour")

    fig
}
