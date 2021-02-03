CPal_default = c(
  "#1F77B4", "#FF7F0E", "#279E68", "#D62728", "#AA40FC", "#8C564B", "#E377C2", "#B5BD61", "#17BECF", "#AEC7E8",
  "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5", "#C49C94", "#F7B6D2", "#DBDB8D", "#9EDAE5", "#AD494A", "#8C6D31",
  "#E31A1C", "#FFD700", "#771122", "#777711", "#1F78B4", "#68228B", "#AAAA44", "#60CC52", "#771155", "#DDDD77",
  "#774411", "#AA7744", "#AA4455", "#117744", "#000080", "#44AA77", "#AA4488", "#DDAA77", "#D9D9D9", "#BC80BD",
  "#FFED6F", "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77",
  "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#A6CEE3", "#B2DF8A", "#33A02C", "#FB9A99",
  "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6",
  "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE",
  "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFFF33", "#A65628", "#F781BF", "#999999",
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3",
  "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5"
)

.default_ggtheme <-  ggplot2::theme(axis.title = element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        legend.title = ggplot2::element_blank(),
        legend.background = ggplot2::element_blank(),
        legend.key = ggplot2::element_blank(),
        plot.margin = grid::unit(c(1,1,1,1),"lines"))

#' Main ACTIONet plotting functions
#'
#' @param ace 'ACTIONetExperiment' object or numeric matrix.
#' @param label_attr Character vector of length NROW(ace) or colname of 'colData(ace)' containing cell labels of interest (clusters, celltypes, etc.).
#' @param trans_attr Numeric vector of length NROW(ace) or colname of 'colData(ace)' used to compute point transparency (smaller values are more transparent).
#' @param trans_fac Transparency modifier (default:1.5)
#' @param trans_th Minimum Z-score for which points with 'scale(trans_attr) < trans_th' are masked.
#' @param node_size Size of nodes in the ACTIONet plot
#' @param palette Color palette. character vector of hex colors or a palette name to pass to 'ggpubr::get_palette()').
#' @param add.text Whether or not to add labels on the top of ACTIONet plot
#' @param show_legend Show legend. Ignored if 'label_attr=NULL' (default:'TRUE').
#' @param title Main title of the plot
#' @param stroke_contrast_fac Factor by which to darken point border.
#' @param add.states Whether or not to include interpolated cell state positions in the plot
#' @param coordinate_attr Name of entry in colMaps(ace) containing the plot coordinates (default:'ACTIONet2D')
#'
#' @return Visualized ACTIONet
#'
#' @examples
#' ace = run.ACTIONet(sce)
#' plot.ACTIONet(ace, ace$assigned_archetype, transparency.attr = ace$node_centrality)
#' @export
.plot_ACTIONet2D <- function(
  ace,
  label_attr = NULL,
  color_attr = NULL,
  trans_attr = NULL,
  trans_fac = 1.5,
  trans_th = -0.5,
  point_size = 1,
  stroke_size = point_size * 0.1,
  palette = CPal_default,
  add.text = TRUE,
  show_legend = TRUE,
  title = "",
  stroke_contrast_fac = 0.1,
  add.backbone = FALSE,
  arch.size.factor = 3,
  coordinate_attr = "ACTIONet2D",
  color_slot = "denovo_color"
) {




    text.halo.width = 0.1
    label_size = 0.8






    plot_coors = .get_plot_coors(ace, coordinate_attr)
    plot_labels = .get_plot_labels(label_attr, ace)
    plot_fill_col = .get_plot_colors(color_attr, plot_labels, ace, color_slot, palette)
    plot_alpha = .get_plot_transparency(trans_attr, ace, trans_fac, trans_th, TRUE)
    plot_border_col = colorspace::darken(plot_fill_col, stroke_contrast_fac)

    if(is.null(plot_labels)){
      plot_labels = "NA"
      show_legend = FALSE
      names(plot_fill_col) = "NA"
      legend_labels = NULL
      legend_fill_breaks = NULL
    } else {
      names(plot_fill_col) = plot_labels
      legend_labels = sort(unique(plot_labels))
      legend_fill_breaks = plot_fill_col[legend_labels]
    }

    plot_data = data.frame(plot_coors,
      labels = plot_labels,
      fill = plot_fill_col,
      color = plot_border_col,
      trans = plot_alpha
    )

    # x = coors[, 1]
    # y = coors[, 2]
    # x.min = min(x)
    # x.max = max(x)
    # y.min = min(y)
    # y.max = max(y)
    # x.min = x.min - (x.max - x.min)/20
    # x.max = x.max + (x.max - x.min)/20
    # y.min = y.min - (y.max - y.min)/20
    # y.max = y.max + (y.max - y.min)/20
    # XL = c(x.min, x.max)
    # YL = c(y.min, y.max)

    p_out <- ggplot2::ggplot() +
         ggplot2::geom_point(
           data = plot_data,
           mapping = aes(
             x = x,
             y = y,
             color = color,
             fill = fill,
             alpha = trans
           ),
           shape = 21,
           size = point_size,
           stroke = stroke_size,
           show.legend = show_legend
         ) + scale_fill_identity(
           guide = "legend",
           breaks = legend_fill_breaks,
           labels = legend_labels
         ) +
         scale_color_identity() +
         scale_alpha_identity() +
         .default_ggtheme

    # rand.perm = sample(nrow(coors))
    # graphics::plot(
    #   coors[rand.perm, c(1, 2)],
    #   pch = 21,
    #   cex = node_size,
    #   bg = vCol[rand.perm],
    #   col = vCol.border[rand.perm],
    #   axes = FALSE,
    #   xlab = "",
    #   ylab = "",
    #   main = title,
    #   xlim = XL,
    #   ylim = YL
    # )

    # if (add.text == TRUE & (!is.null(Annot))) {
    #     graphics::par(xpd = TRUE, mar = graphics::par()$mar * c(1.1, 1.1, 1.1, 1.1))
    #
    #     centroids = Matrix::t(sapply(Annot, function(l) {
    #         idx = which(names(labels) == l)
    #         if (length(idx) == 1) {
    #             return(as.numeric(coors[idx, ]))
    #         }
    #
    #         sub.coors = coors[idx, ]
    #         anchor.coor = as.numeric(apply(sub.coors, 2, function(x) mean(x, trim = 0.8)))
    #
    #         return(anchor.coor)
    #     }))

    #     layout.labels(
    #       x = centroids[, 1],
    #       y = centroids[, 2],
    #       labels = Annot,
    #       col = colorspace::darken(Pal, 0.5),
    #       bg = "#eeeeee",
    #       r = text.halo.width,
    #       cex = label_size
    #     )
    # }

    # if ((suppress.legend == FALSE) & !is.null(Annot)) {
    #     xmin <- graphics::par("usr")[1]
    #     xmax <- graphics::par("usr")[2]
    #     ymin <- graphics::par("usr")[3]
    #     ymax <- graphics::par("usr")[4]
    #
    #     lgd <- graphics::legend(
    #       x = mean(c(xmin, xmax)),
    #       y = mean(c(ymin, ymax)),
    #       legend = Annot,
    #       fill = Pal,
    #       cex = 0.5,
    #       bty = "n",
    #       plot = FALSE
    #     )
    #
    #     graphics::par(xpd = TRUE, mai = c(0, 0, 0, lgd$rect$w))
    #
    #     graphics::legend(
    #       x = xmax,
    #       y = ymin + lgd$rect$h,
    #       legend = Annot,
    #       fill = Pal, cex = 0.5,
    #       bty = "n",
    #       plot = TRUE
    #     )
    #
    # }
}


.get_plot_coors <- function(X, coordinate_attr = NULL, scale_coors = TRUE){

  if (class(X) == "ACTIONetExperiment") {
      if (!is.null(coordinate_attr)) {
          coors = as.matrix(colMaps(X)[[coordinate_attr]])
      } else {
          err = sprintf("'coordinate_attr' cannot be NULL if 'ace' is 'ACTIONetExperiment'.\n")
          stop(err)
      }
  } else {
      if (is.matrix(X) | is.sparseMatrix(X)) {
          coors = as.matrix(X)
      } else {
          err = sprintf("'X' must be 'ACTIONetExperiment' or matrix.\n")
          stop(err)
      }
  }

  if(scale_coors == TRUE) {
    coors = scale(coors)
  }

  coors = data.frame(coors[, c(1,2)])
  colnames(coors) = c("x", "y")

  return(coors)
}

.get_plot_labels <- function(label_attr, ace = NULL){

  # if(is.null(label_attr)){
  #   return(NULL)
  # } else if (length(label_attr) == 1) {
  #   if(is.null(ace)) {
  #     err = sprintf("'ace' cannot be NULL.\n")
  #     stop(err)
  #   }
  #   plot_labels =  = .get_attr_or_split_idx(ace, attr = label_attr, return_vec = TRUE)
  # } else {
  #   plot_labels = .preprocess_annotation_labels(label_attr)
  # }
  #
  # plot_labels =  = .get_attr_or_split_idx(ace, attr = label_attr, return_vec = TRUE)
  # return(plot_labels)

  if(is.null(label_attr)){
    return(NULL)
  }

  plot_labels = .get_attr_or_split_idx(ace, attr = label_attr, return_vec = TRUE)

  return(plot_labels)

}


.get_plot_colors <- function(
  color_attr,
  plot_labels,
  ace,
  color_slot = "denovo_color",
  palette = CPal_default
){


  if(!is.null(color_attr)) {

    if(is.matrix(color_attr) || is.data.frame(color_attr)){

      if(NCOL(color_attr) >= 3) {
        plot_colors = grDevices::rgb(red = color_attr)
      } else if(NCOL(color_attr) == 1) {
        plot_colors = c(color_attr)
      }

    } else if (is.character(color_attr)) {

      if (length(color_attr) == 1) {
        plot_colors = .get_attr_or_split_idx(ace, attr = color_attr, return_vec = TRUE)
      } else {
        plot_colors = color_attr
      }

    } else {
      err = sprint("Invalid 'color_attr'.\n")
      stop(err)
    }

  } else if(!is.null(plot_labels)) {

    label_names = sort(unique(plot_labels))
    num_unique = length(label_names)

    if (num_unique == 1) {
      plot_colors = .default_colors(NROW(ace))
    } else {

      if (length(palette) == 1) {
        plot_palette = ggpubr::get_palette(palette, num_unique)
      } else if(length(palette) < num_unique) {
        plot_palette = CPal_default[1:num_unique]
        msg = sprintf("Not enough colors in 'palette'. Using default palette.\n")
        message(msg)
      } else {
        plot_palette = palette[1:num_unique]
      }

      names(plot_palette) = label_names
      plot_colors = plot_palette[match(plot_labels, names(plot_palette))]

    }

  } else {

    if (class(ace) == "ACTIONetExperiment"){
      plot_colors = grDevices::rgb(colMaps(ace)[[color_slot]])
    } else {
      plot_colors = .default_colors(NROW(ace))
    }

  }

  return(plot_colors)
}


.get_plot_transparency <- function(
  trans_attr,
  ace,
  trans_fac = 1.5,
  trans_th = -0.5,
  scale = TRUE
){


  if(is.null(trans_attr)) {
    return(1)
  } else {

    if (length(trans_attr) == 1) {
      alpha_fac = .get_attr_or_split_idx(ace, attr = trans_attr, return_vec = TRUE)
    } else {
      alpha_fac = trans_attr
    }

    if(scale == TRUE)
      z = scale(alpha_fac)
    else
      z = alpha_fac

    beta_val = 1/(1 + exp(-trans_fac * (z - trans_th)))
    beta_val[z > trans_th] = 1
    beta_val = beta_val^trans_fac

    # plot_alpha = beta_val/max(beta_val)

    return(beta_val)
  }
}


.default_colors <- function(l){
  plot_colors = rep("tomato", l)
  return(plot_colors)
}
