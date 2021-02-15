#' Plots gradient of (potentially imputed) values on ACTIONet scatter plot.
#'
#' @param ace 'ACTIONetExperiment' object
#' @param x Numeric vector of length NCOL(ace).
#' @param alpha_val Smoothing parameter for PageRank imputation of 'x'. No imputation if 'alpha_val=0' (default:0.85).
#' @param log_scale Logical value for whether to log-scale values of 'x' (default:'FALSE').
#' @param nonparameteric If 'FALSE', values of 'x' are used as breaks for color gradient. If 'TRUE' the ranks of the values of 'x' are used instead  (default:'FALSE').
#' @param trans_attr Numeric vector of length NROW(ace) or colname of 'colData(ace)' used to compute point transparency. Smaller values are more transparent.
#' @param trans_fac Transparency modifier (default:1.5).
#' @param trans_th Minimum Z-score for which points with 'scale(trans_attr) < trans_th' are masked (default:-0.5).
#' @param point_size Size of points in ggplot (default:1).
#' @param stroke_size Size of points outline (stroke) in ggplot (default:point_size*0.1).
#' @param stroke_contrast_fac Factor by which to darken point outline for contrast (default:0.1).
#' @param grad_palette Gradient color palette. one of ("greys", "inferno", "magma", "viridis", "BlGrRd", "RdYlBu", "Spectral") or value to pass as 'colors' argument to 'grDevices::colorRampPalette()'.
#' @param net_attr Name of entry in colNets(ace) containing the ACTIONet adjacency matrix to use for value imputation if 'alpha_val>0' (default:'ACTIONet').
#' @param coordinate_attr Name of entry in colMaps(ace) containing the 2D plot coordinates (default:'ACTIONet2D').

#'
#' @return 'ggplot' object.
#'
#' @examples
#' ace = run.ACTIONet(ace)
#' x = logcounts(ace)['CD14', ]
#' plot.ACTIONet.gradient(ace, x, trans_attr = ace$node_centrality)
#' @export
plot.ACTIONet.gradient.test <- function(
  ace,
  x,
  alpha_val = 0.85,
  log_scale = FALSE,
  nonparameteric = FALSE,
  trans_attr = NULL,
  trans_fac = 1.5,
  trans_th = -0.5,
  point_size = 1,
  stroke_size = point_size * 0.1,
  stroke_contrast_fac = 0.1,
  grad_palette = "magma",
  net_attr = "ACTIONet",
  coordinate_attr = "ACTIONet2D"
) {

    NA_col = "#eeeeee"

    ## Create color gradient generator
    if (grad_palette %in% c("greys", "inferno", "magma", "viridis", "BlGrRd", "RdYlBu", "Spectral")) {

        grad_palette = switch(grad_palette,
          greys = grDevices::gray.colors(100),
          inferno = viridis::inferno(500, alpha = 0.8),
          magma = viridis::magma(500, alpha = 0.8),
          viridis = viridis::viridis(500, alpha = 0.8),
          BlGrRd = grDevices::colorRampPalette(c("blue", "grey", "red"))(500),
          Spectral = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "Spectral"))))(100),
          RdYlBu = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"))))(100)
        )

    } else {
        grad_palette = grDevices::colorRampPalette(c(NA_col, grad_palette))(500)
    }

    ## Scale/prune scores, if needed
    if(any(x < 0)){x = x + -1*min(x)}

    x = x - min(x)

    if (log_scale == TRUE)
        x = log1p(x)

    if (alpha_val > 0) {
        x = as.numeric(compute_network_diffusion(
          G = colNets(ace)[[net_attr]],
          X0 = as(as.matrix(x), "sparseMatrix")
        ))
    }

    col_func = (scales::col_bin(
      palette = grad_palette,
      domain = NULL,
      na.color = NA_col,
      bins = 7
    ))

    if (nonparameteric == TRUE)
        plot_fill_col = col_func(rank(x))
    else
        plot_fill_col = col_func(x)

    idx = order(x, decreasing = FALSE)

    p_out <- plot.ACTIONet(
      data = ace,
      label_attr = NULL,
      color_attr = plot_fill_col,
      trans_attr = trans_attr,
      trans_fac = trans_fac,
      trans_th = trans_th,
      point_size = point_size,
      stroke_size = stroke_size,
      stroke_contrast_fac = stroke_contrast_fac,
      palette = NULL,
      add_text_labels = FALSE,
      coordinate_attr = coordinate_attr,
      plot_order = idx
    )

    return(p_out)
}
