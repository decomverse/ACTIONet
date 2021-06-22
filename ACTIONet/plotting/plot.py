from typing import Optional, Union
import plotly as pl
import plotly.io as pio
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly.express as px
import numpy as np
from adjustText import adjust_text
from anndata import AnnData
from scipy import sparse
import pandas as pd
from .color import *
from .palettes import palette_default
from .. import misc_utils as ut
from . import plot_utils as pu
import _ACTIONet as _an

pio.orca.config.use_xvfb = True
pio.orca.config.save()


def validate_plot_params(adata, coordinate_key, label_key, transparency_key):
    if coordinate_key not in adata.obsm.keys():
        raise ValueError(
            f'Did not find adata.obsm[\'{coordinate_key}\']. '
            'Please run nt.layout_ACTIONet() first.'
        )
    if label_key is not None and label_key not in adata.obs.columns:
        raise ValueError(f'Did not find adata.obs[\'{label_key}\'].')
    if transparency_key is not None and transparency_key not in adata.obs.columns:
        raise ValueError(f'Did not find adata.obs[\'{transparency_key}\'].')
    if transparency_key is not None and pd.api.types.is_numeric_dtype(adata.obs[transparency_key].dtype) is False:
        raise ValueError(f'transparency_key must refer to a numeric values, which is not the case for[\'{transparency_key}\'].')


def plot_ACTIONet(
        data: Union[AnnData, pd.DataFrame, np.ndarray],
        label_attr: Union[str, list, pd.Series, None] = None,
        color_attr: Union[str, list, pd.Series, pd.DataFrame, np.ndarray, None] = None,
        trans_attr: Union[str, list, pd.Series, np.ndarray, None] = None,
        trans_fac: Optional[float] = 1.5,
        trans_th: Optional[float] = -0.5,
        point_size: Optional[float] = 3,
        stroke_size: Optional[float] = 0.3,
        stroke_contrast_fac: Optional[float] = 1.2,
        palette: Union[str, list, pd.Series, dict] = palette_default,
        show_legend: Optional[bool] = True,
        hover_text: Union[list, pd.Series, np.ndarray] = None,
        plot_3d: Optional[bool] = False,
        point_order: Union[list, pd.Series, np.ndarray] = None,
        coordinate_attr: Optional[str] = None,
        color_key: Optional[str] = "denovo_color",
) -> go.Figure:
    """Creates an interactive ACTIONet plot with plotly
    :param data:AnnData object with coordinates in '.obsm[coordinate_attr]' or numeric matrix of X-Y(-Z) coordinates. \
        If data is AnnData, 'coordinate_attr' defaults to 'ACTIONet3D' if 'plot_3d=True' or 'ACTIONet2D' if otherwise.
    :param label_attr: list-like object of length data.shape[0] or key of '.obs' containing cell labels of interest (clusters, cell types, etc.).
    :param color_attr: list-like object of length data.shape[0], matrix-like object of RGB values, or key of '.obs' containing point-wise color mappings.
    :param trans_attr: list-like object of length data.shape[0] or key of '.obs' of relative numerical values for computing point transparency. \
        Smaller values are more transparent.
    :param trans_fac: Transparency modifier (default:1.5)
    :param trans_th:Minimum transparency Z-score under which points are masked (default:-0.5).
    :param point_size:Size of points in plotly figure (default:3).
    :param stroke_size: Size of points outline (stroke) in plotly figure (default:0.3).
    :param stroke_contrast_fac: Factor by which to lighten (if < 1) darken (if > 1) point outline for contrast (default:1.2).
    :param palette: color palette for labeled data. One of the following \
        list-like object of color values to assign to each plotly trace alphabetically by label. \
        dict of color values with keys corresponding to members of 'plot_labels'.
    :param show_legend: Show legend for labeled data. Ignored if 'label_attr=None'.
    :param hover_text: list-like object of length data.shape[0] pto use for plotly figure hover text. \
        If defaults to point values given by 'label_attr' or point index if 'label_attr=None'
    :param plot_3d: Visualize plot in 3D using 'scatter3D' (default:'FALSE'). \
        If data is AnnData and 'coordinate_attr=None', 'coordinate_attr' defaults to 'ACTIONet3D'. \
        If data is matrix-like, it must have at least 3 columns.
    :param point_order: Numeric list=like object specifying order in which to plot individual points (default:None). \
        If None, points are plotted in random order.
    :param coordinate_attr: If 'data' is AnnData, key of '.obsm' pertaining to plot coordinates.
    :param color_key:If data is AnnData, key of '.obsm' containing point-wise RGB color mappings (default:'denovo_color'). \
    Used only if no other color mapping parameters are given.
    ...

    :return plotly figure 
    """

    if plot_3d:
        coor_dims = 3
        if isinstance(data, AnnData) and coordinate_attr is None:
            coordinate_attr = "ACTIONet3D"
    else:
        coor_dims = 2
        if isinstance(data, AnnData) and coordinate_attr is None:
            coordinate_attr = "ACTIONet2D"

    plot_coors = pu.get_plot_coors(
        data=data,
        coordinate_attr=coordinate_attr,
        scale_coors=True,
        coor_dims=coor_dims
    )
    plot_labels = pu.get_plot_labels(label_attr=label_attr, data=data)

    if plot_labels is None:
        plot_labels = pd.Series("NA", index=range(plot_coors.shape[0]), name="labels", dtype=str)
    else:
        plot_labels = pd.Series(plot_labels, name="labels", dtype=str)

    plot_data = pd.concat([plot_coors, plot_labels], axis=1)
    plot_data["idx"] = range(plot_data.shape[0])

    if label_attr is None or any(elem is not None for elem in [color_attr, trans_attr]):
        if hover_text is None:
            plot_data["text"] = plot_data["idx"]
        else:
            plot_data["text"] = pd.Series(hover_text, dtype=str)

        plot_data["fill"] = pu.get_plot_colors(
            color_attr=color_attr,
            plot_labels=plot_labels,
            data=data,
            color_key=color_key,
            palette=palette,
            return_dict=False
        )
        plot_data["color"] = [lighten_color(c, stroke_size) for c in plot_data["fill"]]

        plot_data["trans"] = pu.get_plot_transparency(
            trans_attr=trans_attr,
            adata=data,
            trans_fac=trans_fac,
            trans_th=trans_th,
            scale=True
        )

        plot_data["fill"] = append_alpha_to_rgb(plot_data["fill"], plot_data["trans"], unzip_colors=True)
        plot_data["color"] = append_alpha_to_rgb(plot_data["color"], plot_data["trans"], unzip_colors=True)

        if point_order is None:
            plot_data = plot_data.sample(frac=1).reset_index(drop=True)
        else:
            plot_data["pidx"] = point_order
            plot_data = plot_data.sort_values(by="pidx").reset_index(drop=True)

        p = pu.make_plotly_scatter_single_trace(
            x=plot_data["x"],
            y=plot_data["y"],
            z=plot_data["z"],
            label_attr=plot_data["labels"],
            cols_point=plot_data["fill"],
            cols_stroke=plot_data["color"],
            point_size=point_size,
            stroke_size=stroke_size,
            hover_text=plot_data["text"],
            plot_3d=plot_3d
        )

    else:
        if hover_text is None:
            plot_data["text"] = plot_data["labels"]
        else:
            plot_data["text"] = pd.Series(hover_text, dtype=str)

        fill_dict = pu.get_plot_colors(
            color_attr=color_attr,
            plot_labels=plot_labels,
            data=data,
            color_key=color_key,
            palette=palette,
            return_dict=True
        )

        stroke_dict = {k: lighten_color(v, stroke_contrast_fac) for (k, v) in fill_dict.items()}

        if point_order is None:
            plot_data = plot_data.sample(frac=1).reset_index(drop=True)
        else:
            plot_data["pidx"] = point_order
            plot_data = plot_data.sort_values(by="pidx").reset_index(drop=True)

        p = pu.make_plotly_scatter_split_trace(
            x=plot_data["x"],
            y=plot_data["y"],
            z=plot_data["z"],
            label_attr=plot_data["labels"],
            fill_dict=fill_dict,
            stroke_dict=stroke_dict,
            show_legend=show_legend,
            hover_text=plot_data["text"],
            plot_3d=plot_3d
        )

    return p


def plot_ACTIONet_gradient(
        adata: AnnData,
        x: Optional[list] = None,
        coordinate_key: Optional[str] = 'ACTIONet2D',
        transparency_key: Optional[str] = None,
        transparency_z_threshold: Optional[float] = -0.5,
        transparancey_factor: Optional[float] = 3,
        alpha_val: Optional[float] = 0.85,
        node_size: Optional[float] = 1,
        add_text: Optional[bool] = True,
        palette: Optional[list] = "Inferno",
        title: Optional[str] = "",
        nonparametric: Optional[bool] = False,
        output_file: Optional[str] = None
) -> go.Figure:
    """
    Projects a given continuous score on the ACTIONet plot
    Parameters
    ----------
    adata:
        ACTIONet output object
    x:
        score vector
    transparancey_key:
        additional continuous attribute to project onto the transparency of nodes
    transparency_z_threshold:
        controls the effect of transparency mapping
    transparancy_factor:
        controls the effect of transparancy mapping
    node_size:
        Size of nodes in the ACTIONet plot
    palette:
        Color palette (named vector or a name for a given known palette)
    coordinate_key:
       Entry in colMaps(ace) containing the plot coordinates (default:'ACTIONet2D')
    alpha_val:
        alpha_val Between [0, 1]. If it is greater than 0, smoothing of scores would be performed
    output_file:
        filename to save plot (optional)
    """

    G = G.astype(dtype=np.float64)

    np.amin(x)

    if log_scale:
        x = np.log1p(x)

    if alpha_val > 0:
        x = _an.compute_network_diffusion_fast(
            G=G,
            X0=sparse.csc_matrix(x)
        )

# def layout_labels(
#         X: np.ndarray,
#         Y: np.ndarray,
#         fig: go.Figure,
#         labels: list,
#         fontsize: Optional[int] = 18,
#         color: Optional[str] = "#0000FF",
#         background: Optional[str] = "#FFFFFF") -> None:
#     texts = []
#     if isinstance(color, tuple):
#         color = [color] * len(labels)
#     for x, y, label, c in zip(X, Y, labels, color):
#         fig.add_annotation(
#             x=x,
#             y=y,
#             text=str(label),
#             bgcolor=background,
#             font=dict(
#                 family="sans serif",
#                 size=fontsize,
#                 color=color))


# def plot_ACTIONet(
#         adata: AnnData,
#         label_key: Optional[str] = None,
#         coordinate_key: Optional[str] = 'ACTIONet2D',
#         transparency_key: Optional[str] = None,
#         transparency_z_threshold: Optional[float] = -0.5,
#         transparency_factor: Optional[float] = 1.5,
#         border_contrast_factor: Optional[float] = 0.1,
#         node_size: Optional[float] = None,
#         add_text: Optional[bool] = True,
#         palette: Optional[list] = None,
#         output_file: Optional[str] = None):
#     # validate supplied plot parameters
#     validate_plot_params(adata, coordinate_key, label_key, transparency_key)
#     coordinates = ut.scale_matrix(adata.obsm[coordinate_key])
#
#     # get mapping of points to colors
#     if label_key is None:
#         v_col = [(r, g, b) for r, g, b in adata.obsm['denovo_color']]
#     else:
#         labels = adata.obs[label_key]
#         unique_labels = sorted(np.unique(labels))
#         if palette is None:
#             if len(unique_labels) <= len(palette_default):
#                 palette = [hex_to_rgb(color) for color in palette_default]
#             else:
#                 palette = [hex_to_rgb(color) for color in palette_default]
#         elif isinstance(palette, str):
#             # get plotly color palette from string specification
#             assert palette.lower() in px.colors.named_colorscales()
#             if palette in px.colors.qualitative.__dict__.keys():
#                 palette = px.colors.qualitative.__dict__[palette]
#             else:
#                 palette = px.colors.sequential.__dict__[palette]
#         # make sure palette colors are rgb values in the 'rgb(r,g,b)' string format
#         assert len(palette) > 0
#         if type(palette[0]) == str:
#             if palette[0].startswith('#'):
#                 palette = [hex_to_rgb(i) for i in palette]
#         else:
#             palette = [pl.colors.label_rgb(i) for i in palette]
#         label_colors = {label: palette[i % len(palette)] for i, label in enumerate(unique_labels)}
#         v_col = [label_colors[label] for label in labels]
#     v_col_for_border = v_col
#     # get darkened colors for plot marker outline
#     v_col_darkened = [pl.colors.label_rgb(adjust_lightness(color, 1 - border_contrast_factor)) for color in v_col_for_border]
#
#     # calculate transparency
#     if transparency_key is not None:
#         transparency = adata.obs[transparency_key].values
#         z = (transparency - np.mean(transparency)) / np.std(transparency, ddof=1)
#         betas = 1 / (1 + np.exp(- transparency_factor * (z - transparency_z_threshold)))
#         betas[z > transparency_z_threshold] = 1
#         betas **= transparency_factor
#     else:
#         betas = [1] * len(v_col_darkened)
#
#     x = coordinates[:, 0]
#     y = coordinates[:, 1]
#     x_min = np.min(x)
#     x_max = np.max(x)
#     y_min = np.min(y)
#     y_max = np.max(y)
#
#     x_min = x_min - (x_max - x_min) / 20
#     x_max = x_max + (x_max - x_min) / 20
#     y_min = y_min - (y_max - y_min) / 20
#     y_max = y_max + (y_max - y_min) / 20
#
#     # plotly scatterplot
#     if node_size is None:
#         node_size = [1] * coordinates.shape[0]
#
#     fig = make_subplots(rows=1, cols=1)
#     fig.append_trace(go.Scatter(x=coordinates[:, 0],
#                                 y=coordinates[:, 1],
#                                 mode='markers',
#                                 marker=dict(opacity=betas,
#                                             color=v_col,
#                                             line=dict(width=1,
#                                                       color=v_col_darkened))), row=1, col=1)
#     fig.update_layout(yaxis_range=[y_min, y_max])
#     fig.update_layout(xaxis_range=[x_min, x_max])
#     # add cluster labels, if specified
#     if add_text and label_key is not None:
#         labels = adata.obs[label_key]
#         unique_labels = sorted(np.unique(labels))
#
#         colors = []
#         centroids = np.zeros((len(unique_labels), coordinates.shape[1]))
#         for i, label in enumerate(unique_labels):
#             label_coordinates = coordinates[labels == label]
#             centroids[i] = stats.trim_mean(label_coordinates, 0.2, axis=0)
#             cur_text_color = pl.colors.label_rgb(adjust_lightness(palette[i % len(palette)], 0.5))
#             colors.append(cur_text_color)
#             fig.add_annotation(x=centroids[i, 0],
#                                y=centroids[i, 1],
#                                text=str(label),
#                                bgcolor='#FFFFFF',
#                                font=dict(
#                                    family="sans serif",
#                                    size=18,
#                                    color=cur_text_color), row=1, col=1)
#     # save to file if requested by user
#     fig.update_layout(showlegend=False)
#     if not (output_file is None):
#         fig.write_image(output_file)
#     # show the figure
#     # fig.show()
#
#     return fig
#
#
# def plot_ACTIONet_interactive(adata: AnnData,
#                               label_key: Optional[str] = None,
#                               transparency_attribute: Optional[str] = None,
#                               transparancy_z_threshold: Optional[float] = -1,
#                               transparancey_factor=1,
#                               node_size: Optional[float] = 1,
#                               palette: Optional[str] = "CPal20",
#                               enrichment_table: Optional[bool] = False,
#                               top_features: Optional[int] = 7,
#                               blacklist_pattern: Optional[str] = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|MALAT1|B2M|GAPDH",
#                               title: Optional[str] = "ACTIONet",
#                               coordinate_slot: Optional[str] = "ACTIONet2D",
#                               plot_3d: bool = False):
#     """
#     Creates an interactive ACTIONet plot with plotly
#     Parameters
#     ---------
#     ace:
#         ACTIONet output object
#     labels:
#         Annotation of interest (clusters, celltypes, etc.) to be projected on the ACTIONet plot
#     transparency_attr:
#         Additional continuous attribute to project onto the transparency of nodes
#     trans_z_threshold, trans_fact:
#         Control the effect of transparency mapping
#     node_size:
#         Size of nodes in the ACTIONet plot
#     palette:
#         Color palette (named vector or a name for a given known palette)
#     enrichment_table:
#         To project the top-ranked features interactively.
#     top_features:
#         Number of features to show per cell
#     blacklist_pattern:
#         List of genes to filter-out
#     title:
#         Main title of the plot
#     coordinate_slot:
#         coordinate_slot in anndata object where visualization should be stored
#     plot_3d:
#         Whether to show the plot in 3D
#     Returns
#     -------
#     Visualized ACTIONet
#     """
#     # validate supplied plot parameters
#     validate_plot_params(adata, coordinate_key, label_key, transparency_key)
#
#     # get the number of variables stored in obs matrix
#     nV = adata.obs.columns.shape[0]
#
#     # determine whether a 3D plot should be generated
#     if ((coordinate_slot == "ACTIONet2D") and (plot_3d is True)):
#         coordinate_slot = "ACTIONet3D"
#     coordinates = ut.scale_matrix(adata.obsm[coordinate_key])
