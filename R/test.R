#' Main ACTIONet interactive plotting function
#'
#' @param ace ACTIONet output object
#' @param labels Annotation of interest (clusters, celltypes, etc.) to be projected on the ACTIONet plot
#' @param trans_attr Additional continuous attribute to project onto the transparency of nodes
#' @param trans_th, trans_fac Control the effect of transparency mapping
#' @param point_size Size of nodes in the ACTIONet plot
#' @param palette Color palette (named vector or a name for a given known palette)
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
#' @rawNamespace import(plotly, except = 'last_plot')
#' @export
plot.ACTIONet.interactive <- function(
  ace,
  labels = NULL,
  trans_attr = NULL,
  trans_th = -1,
  trans_fac = 1,~`
  point_size = 1,
  palette = CPal_default,
  threeD = FALSE,
  title = "ACTIONet",
  coordinate_slot = "ACTIONet2D"
) {

    nV = ncol(ace)
    point_size = point_size * 3
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
        if (length(palette) > 1) {
            if (length(palette) < length(Annot)) {
              palette = CPal_default
            }
            if (is.null(names(palette))) {
                Pal = palette[1:length(Annot)]
            } else {
                Pal = palette[Annot]
            }
        } else {
            Pal = ggpubr::get_palette(palette, length(Annot))
        }

        names(Pal) = Annot
        vCol = Pal[names(labels)]
    }

    if (!is.null(trans_attr)) {
        z = scale(trans_attr)  # (trans_attr - median(trans_attr))/mad(trans_attr)
        beta = 1/(1 + exp(-trans_fac * (z - trans_th)))
        beta[z > trans_th] = 1
        beta = beta^trans_fac

        vCol.border = scales::alpha(colorspace::darken(vCol, 0.5), beta)
        vCol = scales::alpha(vCol, beta)
    } else {
        vCol.border = colorspace::darken(vCol, 0.5)
    }

    # if (!is.null(enrichment.table)) {
    #     if (ncol(enrichment.table) == nV) {
    #         cell.scores = Matrix::t(enrichment.table)
    #     } else if ((nrow(enrichment.table) != nV)) {
    #         H = colMaps(ace)[["H_unified"]]
    #         if ((nrow(enrichment.table) == nrow(H)) | (ncol(enrichment.table) ==
    #             nrow(H))) {
    #             cell.scores = map.cell.scores.from.archetype.enrichment(ace, enrichment.table)
    #         } else {
    #             cell.scores = NULL
    #         }
    #     } else {
    #         cell.scores = enrichment.table
    #     }
    # } else {
    #     temp.enrichment.table = as.matrix(rowMaps(ace)[["unified_feature_specificity"]])
    #     if (!is.null(row.names(temp.enrichment.table))) {
    #         filtered.rows = grep(blacklist_pattern, rownames(temp.enrichment.table))
    #         if (length(filtered.rows) > 0){
    #           enrichment.table = temp.enrichment.table[-filtered.rows, ]
    #         } else {
    #           enrichment.table = temp.enrichment.table
    #         }
    #
    #         GT = apply(enrichment.table, 2, function(x) rownames(enrichment.table)[order(x,
    #             decreasing = TRUE)[1:min(100, nrow(enrichment.table))]])
    #         selected.features = sort(unique(as.character(GT)))
    #
    #         W = exp(scale(Matrix::t(colMaps(ace)[["H_unified"]])))
    #         cs = fastColSums(W)
    #         W = Matrix::t(scale(W, center = FALSE, scale = cs))
    #
    #         cell.scores = W %*% Matrix::t(enrichment.table[selected.features, ])
    #     } else {
    #         cell.scores = NULL
    #     }
    # }
    node.annotations = as.character(labels)
    # if (!is.null(Alt_Text)) {
    #     node.annotations = Alt_Text
    # } else {
    #     if (!is.null(cell.scores)) {
    #         selected.features = colnames(cell.scores)
    #         node.annotations = apply(cell.scores, 1, function(x) paste(selected.features[order(x,
    #             decreasing = TRUE)[1:top_features]], collapse = "\n"))
    #     } else {
    #         node.annotations = rep("", nV)
    #     }
    # }
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
    node.data$size = point_size

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
                  width = 0.1 * point_size,
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
                  width = 0.1 * point_size,
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
                  width = 0.1 * point_size,
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
                  width = 0.1 * point_size,
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
