CPal_default <- c(
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


.get_plot_coors <- function(X,
                            coordinate_attr = NULL,
                            scale_coors = TRUE) {
  if (is(X, "ACTIONetExperiment")) {
    if (!is.null(coordinate_attr)) {
      if(coordinate_attr %in% names(colMaps(X))) {
        coors <- as.matrix(colMaps(X)[[coordinate_attr]])
      } else {
        err <- sprintf("Attribute '%s' is not in 'colMaps'.\n", coordinate_attr)
        stop(err)
      }
    } else {
      err <- sprintf("'coordinate_attr' cannot be NULL if 'ace' is 'ACTIONetExperiment'.\n")
      stop(err)
    }
  } else {
    if (is.matrix(X) | ACTIONetExperiment:::is.sparseMatrix(X)) {
      coors <- as.matrix(X)
    } else {
      err <- sprintf("'X' must be 'ACTIONetExperiment' or matrix.\n")
      stop(err)
    }
  }

  if (scale_coors == TRUE) {
    coors <- scale(coors)
  }

  coors <- data.frame(coors)
  colnames(coors) <- c("x", "y", "z")[1:NCOL(coors)]

  return(coors)
}


.get_plot_labels <- function(label_attr, data = NULL) {
  if (is.null(label_attr)) {
    return(NULL)
  }

  # if (is(data, "ACTIONetExperiment")) {
  #   plot_labels <- ACTIONetExperiment::get.data.or.split(data, attr = label_attr, to_return = "data")
  # } else {
  #   plot_labels <- label_attr
  # }

  plot_labels <- .validate_attr(
    obj = data,
    attr = label_attr,
    obj_name = "data",
    attr_name = "label_attr",
    match_row = TRUE
  )

  plot_labels <- as.character(plot_labels)

  # if (!is.numeric(plot_labels)) {
  #   plot_labels <- as.character(plot_labels)
  #   # plot_labels[is.na(plot_labels)] <- "NA"
  # }

  return(plot_labels)
}


.get_plot_colors <- function(color_attr,
                             plot_labels,
                             data,
                             color_slot = "denovo_color",
                             palette = CPal_default,
                             NA_color = "#CCCCCC"
                           ) {
  if (is(data, "ACTIONetExperiment")) {
    n_dim <- NCOL(data)
  } else {
    n_dim <- NROW(data)
  }

  if (!is.null(color_attr)) {
    if (is.matrix(color_attr) || is.data.frame(color_attr)) {
      if (NCOL(color_attr) >= 3) {
        plot_colors <- grDevices::rgb(red = color_attr)
      } else if (NCOL(color_attr) == 1) {
        plot_colors <- c(color_attr)
      }
    } else if (is.character(color_attr)) {
        plot_colors <- .validate_attr(
          obj = data,
          attr = color_attr,
          obj_name = "data",
          attr_name = "color_attr",
          match_row = TRUE
        )
      # if (length(color_attr) == 1) {
      #   plot_colors <- ACTIONetExperiment::get.data.or.split(data, attr = color_attr, to_return = "data")
      # } else {
      #   plot_colors <- color_attr
      # }
    } else if (is.numeric(color_attr) && length(color_attr) == n_dim) {
      plot_colors <- color_attr
    }else {
      err <- sprintf("Invalid 'color_attr'.\n")
      stop(err)
    }
  } else if (!is.null(plot_labels)) {
    plot_labels <- as.character(plot_labels)
    label_names <- sort(unique(plot_labels[!is.na(plot_labels)]))
    num_unique <- length(label_names)

    if (num_unique == 1) {
      if(startsWith(palette[1], "#")) {
        plot_colors <- .default_colors(n_dim, col = palette[1])
      } else {
        plot_colors <- .default_colors(n_dim)
      }
      plot_colors[is.na(plot_labels)] = NA_color
    } else {
      if (length(palette) == 1) {
        plot_palette <- ggpubr::get_palette(palette, num_unique)
      } else if (length(palette) < num_unique) {
        plot_palette <- CPal_default[1:num_unique]
        msg <- sprintf("Not enough colors in 'palette'.\n")
        message(msg)
      } else {
        if (!is.null(names(palette))) {
          if (all(label_names %in% names(palette))) {
            plot_palette <- palette[label_names]
          } else {
            plot_palette <- palette[1:num_unique]
          }
        } else {
          plot_palette <- palette[1:num_unique]
        }
      }

      plot_labels[is.na(plot_labels)] <- "NA"
      names(plot_palette) <- label_names
      plot_palette = c(plot_palette, "NA" = NA_color)
      plot_colors <- plot_palette[match(plot_labels, names(plot_palette))]
    }
  } else {
    if (is.null(color_slot)) {
      plot_colors <- .default_colors(n_dim)
    } else if (is(data, "ACTIONetExperiment")) {
      if (color_slot %in% names(colMaps(data))) {
        plot_colors <- grDevices::rgb(colMaps(data)[[color_slot]])
      } else {
        err <- sprintf("%s not in colMaps(ace).\n")
        stop(err)
      }
    } else {
      plot_colors <- .default_colors(n_dim)
    }
  }

  return(plot_colors)
}


.get_plot_transparency <- function(
  trans_attr,
  data,
  trans_fac = 1.5,
  trans_th = -0.5,
  scale = TRUE
) {

  if (is.null(trans_attr)) {
    return(1)
  }

  alpha_fac <- .validate_attr(
    obj = data,
    attr = trans_attr,
    obj_name = "data",
    attr_name = "trans_attr",
    match_row = TRUE
  )

  # alpha_fac <- ACTIONetExperiment::get.data.or.split(ace, attr = trans_attr, to_return = "data")

  if (scale == TRUE) {
    z <- scale(alpha_fac)
  } else {
    z <- alpha_fac
  }

  alpha_val <- 1 / (1 + exp(-trans_fac * (z - trans_th)))
  alpha_val[z > trans_th] <- 1
  alpha_val <- alpha_val^trans_fac

  return(alpha_val)
}


.default_colors <- function(l, col = NULL) {
  if(is.null(col)) {
    col_use ="#FF6347"
  } else if (startsWith(col, "#") && nchar(col) == 7) {
    col_use = col
  } else {
    col_use ="#FF6347"
  }
  plot_colors <- rep(col_use, l)
  return(plot_colors)
}
