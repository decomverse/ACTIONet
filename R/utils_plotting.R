.get_plot_coors <- function(X, coordinate_attr = NULL, scale_coors = TRUE, plot_dims = 2){

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

  coors = data.frame(coors[, 1:plot_dims])
  colnames(coors) = c("x", "y", "z")[1:plot_dims]

  return(coors)
}



.get_plot_labels <- function(label_attr, ace = NULL){

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
