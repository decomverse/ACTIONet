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


.get_plot_coors <- function(
  X,
  coordinate_attr = NULL,
  scale_coors = TRUE,
  plot_dims = 2
){

  if (is(X, "ACTIONetExperiment")) {
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

  coors = data.frame(coors[, 1:plot_dims])
  colnames(coors) = c("x", "y", "z")[1:plot_dims]

  return(coors)
}



.get_plot_labels <- function(label_attr, data = NULL){

  if(is.null(label_attr)){
    return(NULL)
  }

  if (is(data, "ACTIONetExperiment")) {
    plot_labels = .get_attr_or_split_idx(data, attr = label_attr, return_vec = TRUE)
  } else {
    plot_labels = label_attr
  }

  plot_labels = as.character(plot_labels)
  return(plot_labels)

}


.get_plot_colors <- function(
  color_attr,
  plot_labels,
  data,
  color_slot = "denovo_color",
  palette = CPal_default
){

  if (is(data, "ACTIONetExperiment")) {
    n_dim = NCOL(data)
  } else {
    n_dim = NROW(data)
  }

  if(!is.null(color_attr)) {

    if(is.matrix(color_attr) || is.data.frame(color_attr)){

      if(NCOL(color_attr) >= 3) {
        plot_colors = grDevices::rgb(red = color_attr)
      } else if(NCOL(color_attr) == 1) {
        plot_colors = c(color_attr)
      }

    } else if (is.character(color_attr)) {

      if (length(color_attr) == 1) {
        plot_colors = .get_attr_or_split_idx(data, attr = color_attr, return_vec = TRUE)
      } else {
        plot_colors = color_attr
      }

    } else {
      err = sprint("Invalid 'color_attr'.\n")
      stop(err)
    }

  } else if(!is.null(plot_labels)) {
    
    plot_labels = as.character(plot_labels)
    label_names = sort(unique(plot_labels))
    num_unique = length(label_names)

    if (num_unique == 1) {
      plot_colors = .default_colors(n_dim)
    } else {

      if (length(palette) == 1) {
        plot_palette = ggpubr::get_palette(palette, num_unique)
      } else if(length(palette) < num_unique) {
        plot_palette = CPal_default[1:num_unique]
        msg = sprintf("Not enough colors in 'palette'. Using default palette.\n")
        message(msg)
      } else {

        if(!is.null(names(palette))){
          if(all(label_names %in% names(palette))){
            plot_palette = palette[label_names]
          } else {
            plot_palette = palette[1:num_unique]
          }
        } else {
          plot_palette = palette[1:num_unique]
        }

      }

      names(plot_palette) = label_names
      plot_colors = plot_palette[match(plot_labels, names(plot_palette))]

    }

  } else {

    if (is(data, "ACTIONetExperiment")) {
      plot_colors = grDevices::rgb(colMaps(data)[[color_slot]])
    } else {
      plot_colors = .default_colors(n_dim)
    }

  }

  return(plot_colors)
}


.get_plot_transparency <- function(
  trans_attr,
  ace = NULL,
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

    alpha_val = 1/(1 + exp(-trans_fac * (z - trans_th)))
    alpha_val[z > trans_th] = 1
    alpha_val = alpha_val^trans_fac

    return(alpha_val)
  }
}


.default_colors <- function(l){
  plot_colors = rep("tomato", l)
  return(plot_colors)
}


#' @import ggplot2
.layout_plot_labels <- function(
  # p,
  plot_data = NULL,
  label_names = NULL,
  label_colors = NULL,
  darken = TRUE,
  alpha_val = 0.5,
  text_size = 3,
  constrast_fac = 0.5,
  nudge = FALSE,
  use_repel = TRUE
) {

    if(is.null(label_names)){
      label_names = sort(unique(plot_data$labels))
    }

    label_coors = split(plot_data[, 1:2], plot_data$labels)
    label_coors = label_coors[label_names]
    centroids = lapply(label_coors, function(df){
      mat = as.matrix(df)
      cent = apply(mat, 2, function(x) mean(x, trim = 0.25))
      return(cent)
    })

    cent_sd = lapply(label_coors, function(df){
      mat = as.matrix(df)
      cent_sd = apply(mat, 2, sd)
      return(cent_sd)
    })
    cent_sd = do.call(rbind, cent_sd)
    colnames(cent_sd) = c("x_sd", "y_sd")

    if(is.null(label_colors)){
      label_colors = rep("black", length(label_names))
    } else {
        if(darken == TRUE)
          label_colors = colorspace::darken(label_colors, constrast_fac)
    }

    layout_data = data.frame(do.call(rbind, centroids),
                             cent_sd,
                             labels = names(centroids),
                             color = label_colors
                           )

    layout_data[, c("x", "y")] = gplots::space(layout_data$x, layout_data$y, s=c(1/20, 1/5), na.rm=TRUE, direction="y")

    if(nudge == TRUE){
      layout_data$x = layout_data$x + (1 - exp(-0.5 * abs(layout_data$x_sd - max(layout_data$x_sd))))
      layout_data$y = layout_data$y + (1 - exp(-0.5 * abs(layout_data$y_sd - max(layout_data$y_sd))))
    }

    # if(use_repel == TRUE){
    #   p <- p + ggrepel::geom_label_repel(
    #     data = layout_data,
    #     mapping = aes(
    #       x = x,
    #       y = y,
    #       label = labels,
    #       color = color
    #     ),
    #     fill = scales::alpha(c("white"), alpha_val),
    #     size = text_size,
    #     segment.color = 'transparent'
    #   )
    # } else {
    #   p <- p + geom_label(
    #     data = layout_data,
    #     mapping = aes(
    #       x = x,
    #       y = y,
    #       label = labels,
    #       color = color
    #     ),
    #     fill = scales::alpha(c("white"), alpha_val),
    #     size = text_size)
    # }

    if(use_repel == TRUE){
      layer_out <- ggrepel::geom_label_repel(
        data = layout_data,
        mapping = aes(
          x = x,
          y = y,
          label = labels,
          color = color
        ),
        fill = scales::alpha(c("white"), alpha_val),
        size = text_size,
        segment.color = 'transparent'
      )
    } else {
      layer_out <- geom_label(
        data = layout_data,
        mapping = aes(
          x = x,
          y = y,
          label = labels,
          color = color
        ),
        fill = scales::alpha(c("white"), alpha_val),
        size = text_size)
    }

  return(layer_out)
}
