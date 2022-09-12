.axis_params <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)


.set_default_layout_plotly <- function(
  p,
  show_legend = FALSE,
  plot_3d = FALSE
) {

  if (plot_3d == TRUE) {
    p <- plotly::layout(
      p = p,
      scene = list(
        xaxis = .axis_params,
        yaxis = .axis_params,
        zaxis = .axis_params
      ),
      showlegend = show_legend,
      legend = list(
        marker = list(
          marker.size = 10
        )
      )
    )
  } else {
    p <- plotly::layout(
      p = p,
      xaxis = .axis_params,
      yaxis = .axis_params,
      showlegend = show_legend
    )
  }
  return(p)
}


.make_plotly_scatter_single_trace <- function(
  x,
  y,
  z = NULL,
  label_attr = NULL,
  col_vals = NULL,
  cols_stroke = NULL,
  point_size = 3,
  stroke_size = point_size * 0.1,
  show_legend = FALSE,
  hover_text = NULL,
  grad_palette = NULL,
  plot_3d = FALSE
) {

  if (is.null(z)) {
    z <- NA
  }

  plot_data <- data.frame(
    x = x,
    y = y,
    z = z,
    labels = label_attr
  )

  if (is.null(hover_text)) {
    plot_data$text <- 1:NROW(plot_data)
  } else {
    plot_data$text <- hover_text
  }

  if(is.numeric(col_vals)) {
    plot_data$text <- col_vals
    cols_stroke <- col_vals
    if(!is.character(grad_palette) || length(grad_palette) > 1){
      grad_palette = "Viridis"
    }
  } else {
    col_vals <- grDevices::rgb(t(grDevices::col2rgb(col_vals) / 255), alpha = plot_data$trans)
    cols_stroke <- grDevices::rgb(t(grDevices::col2rgb(plot_data$color) / 255), alpha = plot_data$trans)
    grad_palette = NULL
  }

  if (plot_3d == TRUE) {
    p <- plotly::plot_ly(
      data = plot_data,
      x = plot_data$x,
      y = plot_data$y,
      z = plot_data$z,
      marker = list(
        color = col_vals,
        colorscale = grad_palette,
        size = point_size,
        line = list(
          width = stroke_size,
          color = cols_stroke,
          colorscale = grad_palette
        ),
        showscale = show_legend
      ),
      text = plot_data$text,
      hoverinfo = "text",
      mode = "markers",
      type = "scatter3d"
    )
  } else {
    p <- plotly::plot_ly(
      data = plot_data,
      x = plot_data$x,
      y = plot_data$y,
      marker = list(
        color = col_vals,
        colorscale = grad_palette,
        size = point_size,
        line = list(
          width = stroke_size,
          color = cols_stroke,
          colorscale = grad_palette
        ),
        showscale = show_legend
      ),
      text = plot_data$text,
      hoverinfo = "text",
      mode = "markers",
      type = "scattergl"
    )
  }

  p <- .set_default_layout_plotly(p, FALSE, plot_3d)

  # if (show_legend == FALSE){
  #     p <- p %>% hide_colorbar()
  # }

  return(p)
}


.make_plotly_scatter_split_trace <- function(
  x,
  y,
  z = NULL,
  label_attr = NULL,
  cols_fill = NULL,
  cols_stroke = NULL,
  point_size = 3,
  stroke_size = point_size * 0.1,
  show_legend = TRUE,
  hover_text = NULL,
  plot_3d = FALSE
) {

  if (is.null(z)) {
    z <- NA
  }

  plot_data <- data.frame(
    x = x,
    y = y,
    z = z,
    labels = label_attr
  )

  if (is.null(hover_text)) {
    plot_data$text <- plot_data$labels
  } else {
    plot_data$text <- hover_text
  }

  trace_names <- gtools::mixedsort(unique(plot_data$labels))
  if ("NA" %in% trace_names) {
    trace_names = c("NA", trace_names[trace_names != "NA"])
  }

  if (plot_3d == TRUE) {
    p <- plot_ly(type = "scatter3d", mode = "markers")

    for (n in trace_names) {
      sub_data <- plot_data[plot_data$labels == n, ]

      p <- add_trace(
        p = p,
        x = sub_data$x,
        y = sub_data$y,
        z = sub_data$z,
        marker = list(
          color = cols_fill[n],
          size = point_size,
          line = list(
            width = stroke_size,
            color = cols_stroke[n]
          )
        ),
        text = sub_data$text,
        hoverinfo = "text",
        mode = "markers",
        name = n
      )
    }
  } else {
    p <- plotly::plot_ly(type = "scattergl", mode = "markers")

    for (n in trace_names) {
      sub_data <- plot_data[plot_data$labels == n, ]

      p <- plotly::add_trace(p,
        x = sub_data$x,
        y = sub_data$y,
        marker = list(
          color = cols_fill[n],
          size = point_size,
          line = list(
            width = stroke_size,
            color = cols_stroke[n]
          )
        ),
        text = sub_data$text,
        hoverinfo = "text",
        mode = "markers",
        name = n
      )
    }
  }

  p <- .set_default_layout_plotly(p, show_legend, plot_3d) %>% hide_colorbar()
  return(p)
}
