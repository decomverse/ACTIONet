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

    coors = .get_plot_coors(ace, coordinate_attr)


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
