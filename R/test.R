plot.ACTIONet.interactive.test <- function(
  data,
  label_attr = NULL,
  color_attr = NULL,
  trans_attr = NULL,
  trans_fac = 1.5,
  trans_th = -0.5,
  point_size = 3,
  stroke_size = point_size * 0.1,
  stroke_contrast_fac = 0.1,
  stroke_color = NULL,
  palette = CPal_default,
  show_legend = NULL,
  coordinate_attr = "ACTIONet2D",
  color_slot = "denovo_color",
  point_order = NULL,
  hover_text = NULL,
  plot_3d = FALSE
) {

  plot_coors = .get_plot_coors(data, coordinate_attr)
  plot_labels = .get_plot_labels(label_attr, data)
  plot_fill_col = .get_plot_colors(color_attr, plot_labels, data, color_slot, palette)
  plot_alpha = .get_plot_transparency(trans_attr, data, trans_fac, trans_th, TRUE)

  if(is.null(stroke_color))
    plot_border_col = colorspace::darken(plot_fill_col, stroke_contrast_fac)
  else
    plot_border_col = stroke_color

  if(plot_3d == TRUE){
    if (NCOL(plot_coors) < 3){
      if("ACTIONet3D" %in% names(colMaps(data))){
        msg = sprintf("'plot_3d == TRUE' but given coordinates have < 3 columns.\nUsing 'ACTIONet3D'.\n")
        message(msg)
        plot_coors = .get_plot_coors(data, "ACTIONet3D")
      } else {
        err = sprintf("'plot_3d == TRUE' but given coordinates have < 3 columns.\n")
        stop(err)
      }
    }
  }

  plot_data = data.frame(plot_coors,
                         fill = plot_fill_col,
                         color = plot_border_col,
                         trans = plot_alpha,
                         idx = 1:NROW(plot_coors)
  )

  if(is.null(label_attr)){
    show_legend = FALSE
    plot_data$labels = "NA"
  } else {
    plot_data$labels = plot_labels
  }

  if (!is.null(hover_text))
    plot_data$text = hover_text
  else{
    if(is.null(label_attr))
      plot_data$text = plot_data$idx
    else
      plot_data$text = plot_data$labels
  }

  if(is.null(point_order))
    pidx = sample(NROW(plot_data))
  else
    pidx = point_order

  plot_data = plot_data[pidx, ]

  cont_attr = c(color_attr, trans_attr)
  if(is.null(label_attr) | any(!sapply(cont_attr, is.null)) ){

    if(is.null(show_legend))
      show_legend = FALSE

    plot_data$fill = grDevices::rgb(t(grDevices::col2rgb(plot_data$fill)/255), alpha = plot_data$trans)
    plot_data$color = grDevices::rgb(t(grDevices::col2rgb(plot_data$color)/255), alpha = plot_data$trans)

    p <- .make_plotly_scatter_single_trace(
      x = plot_data$x,
      y = plot_data$y,
      z = plot_data$z,
      label_attr = plot_data$labels,
      cols_fill = plot_data$fill,
      cols_stroke = plot_data$color,
      point_size = point_size,
      stroke_size = stroke_size,
      show_legend = show_legend,
      hover_text = plot_data$text,
      plot_3d = plot_3d
    )

  } else {

    if(is.null(show_legend))
      show_legend = TRUE

    col_idx = which(!duplicated(plot_data$labels))
    palette_fill = plot_data$fill[col_idx]
    palette_stroke = plot_data$color[col_idx]
    names(palette_fill) = names(palette_stroke) = plot_data$labels[col_idx]

    p <- .make_plotly_scatter_split_trace(
      x = plot_data$x,
      y = plot_data$y,
      z = plot_data$z,
      label_attr = plot_data$labels,
      cols_fill = palette_fill,
      cols_stroke = palette_stroke,
      point_size = point_size,
      stroke_size = stroke_size,
      show_legend = show_legend,
      hover_text = plot_data$text,
      plot_3d = plot_3d
    )

  }

  return(p)
}
