.default_ggtheme <- ggplot2::theme(
  axis.title = ggplot2::element_blank(),
  axis.text = ggplot2::element_blank(),
  axis.ticks = ggplot2::element_blank(),
  panel.grid = ggplot2::element_blank(),
  panel.background = ggplot2::element_blank(),
  plot.background = ggplot2::element_rect(fill = "white", color = NA),
  legend.title = ggplot2::element_blank(),
  legend.background = ggplot2::element_blank(),
  legend.key = ggplot2::element_blank(),
  plot.margin = grid::unit(c(1, 1, 1, 1), "lines")
)


#' @import ggplot2
.layout_plot_labels <- function(
  plot_data = NULL,
  label_names = NULL,
  label_colors = NULL,
  darken = TRUE,
  alpha_val = 0.5,
  text_size = 3,
  constrast_fac = 0.5,
  nudge = FALSE,
  use_repel = TRUE,
  repel_force = 0.05
) {

  if (is.null(label_names)) {
    label_names <- sort(unique(plot_data$labels))
  }

  label_coors <- split(plot_data[, 1:2], plot_data$labels)
  label_coors <- label_coors[label_names]
  centroids <- lapply(label_coors, function(df) {
    mat <- as.matrix(df)
    cent <- apply(mat, 2, function(x) mean(x, trim = 0.45))
    return(cent)
  })

  cent_sd <- lapply(label_coors, function(df) {
    mat <- as.matrix(df)
    cent_sd <- apply(mat, 2, sd)
    return(cent_sd)
  })
  cent_sd <- do.call(rbind, cent_sd)
  colnames(cent_sd) <- c("x_sd", "y_sd")

  if (is.null(label_colors)) {
    label_colors <- rep("black", length(label_names))
  } else {
    if (darken == TRUE) {
      label_colors <- colorspace::darken(label_colors, constrast_fac)
    }
  }

  layout_data <- data.frame(do.call(rbind, centroids),
    cent_sd,
    labels = names(centroids),
    color = label_colors
  )

  layout_data[, c("x", "y")] <- gplots::space(layout_data$x, layout_data$y, s = c(1 / 20, 1 / 5), na.rm = TRUE, direction = "y")

  if (nudge == TRUE) {
    layout_data$x <- layout_data$x + (1 - exp(-0.5 * abs(layout_data$x_sd - max(layout_data$x_sd))))
    layout_data$y <- layout_data$y + (1 - exp(-0.5 * abs(layout_data$y_sd - max(layout_data$y_sd))))
  }

  if (use_repel == TRUE) {
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
      #box.padding = 0.5,
      #max.overlaps = Inf,
      #min.segment.length = 0.1,
      force = repel_force
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
      size = text_size
    )
  }

  return(layer_out)
}


.plot_arrange_dim <- function(n) {
  s <- sqrt(n)
  sf <- round(s)
  sr <- s - sf

  if (sr <= 0) {
    d <- c(sf, sf)
  } else {
    d <- c(sf, sf + 1)
  }

  return(d)
}

.get_ggplot_scale <- function(
  p,
  col_vals = NULL,
  grad_palette = NULL,
  legend_point_size = 3,
  legend_text_size = 10,
  stroke_contrast_fac = 0.1,
  legend_labels = NULL,
  legend_fill_colors = NULL
) {

  NA_color = "#CCCCCC"

  if (is.numeric(col_vals)) {
    p_out <- p +
    scale_fill_gradientn(
      colors = grad_palette,
      na.value = NA_color,
      guide = "colourbar",
      aesthetics = "fill"
    ) +
    scale_colour_gradientn(
      colors = colorspace::darken(grad_palette, stroke_contrast_fac),
      na.value = NA_color,
      guide = NULL,
      aesthetics = "colour"
    )
  } else {
    p_out <- p + scale_fill_identity(
      guide = "legend",
      labels = legend_labels,
      breaks = legend_fill_colors
    ) +
    scale_color_identity() +
    guides(
      fill = guide_legend(override.aes = list(size = legend_point_size))
    )
  }
  p_out <- p_out + theme(legend.text = element_text(size = legend_text_size))
  return(p_out)
}
