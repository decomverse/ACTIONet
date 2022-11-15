from random import sample
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.io as pio
import scanpy as sc
from anndata import AnnData
from cycler import Cycler
from matplotlib.axes import Axes
from matplotlib.colors import Colormap
from matplotlib.figure import Figure

from ACTIONet.network.diffusion import diffusion
from ACTIONet.plotting.color import append_alpha_to_rgb, lighten_color
from ACTIONet.plotting.palettes import palette_default

from . import utils as pu

pio.orca.config.use_xvfb = True
pio.orca.config.save()


def validate_plot_params(adata, coordinate_key, label_key, transparency_key):
    if coordinate_key not in adata.obsm.keys():
        raise ValueError(f"Did not find adata.obsm['{coordinate_key}']. " "Please run nt.layoutNetwork() first.")
    if label_key is not None and label_key not in adata.obs.columns:
        raise ValueError(f"Did not find adata.obs['{label_key}'].")
    if transparency_key is not None and transparency_key not in adata.obs.columns:
        raise ValueError(f"Did not find adata.obs['{transparency_key}'].")
    if transparency_key is not None and pd.api.types.is_numeric_dtype(adata.obs[transparency_key].dtype) is False:
        raise ValueError(f"transparency_key must refer to a numeric values, which is not the case for['{transparency_key}'].")


def plot_ACTIONet(
    adata: AnnData,
    annotation: Optional[Union[str, list, pd.Series]],
    projection: sc._utils.Literal["2d", "3d"] = "2d",
    palette: Union[str, Sequence[str], Cycler, None] = palette_default,
    add_outline: Optional[bool] = False,
    frameon: Optional[bool] = False,
    size=5,
    legend_fontsize="small",
    legend_loc: str = "on data",
    return_fig: Optional[bool] = False,
    show: Optional[bool] = False,
    **kwargs,
) -> Union[Figure, Axes, None]:

    if projection == "2d":
        adata.obsm["X_actionet2d"] = adata.obsm["ACTIONet2D"]
        basis = "actionet2d"
    else:
        adata.obsm["X_actionet3d"] = adata.obsm["ACTIONet3D"]
        basis = "actionet3d"

    tmp_key = "__annotations__"
    if annotation is not None:
        if len(annotation) == adata.shape[0]:
            adata.obs[tmp_key] = pd.Series(annotation, index=adata.obs.index.values)
        else:
            adata.obs[tmp_key] = pd.Series(
                adata.obs[annotation].astype("str"),
                index=adata.obs.index.values.astype("str"),
            )

    p = sc.pl.embedding(
        adata,
        basis=basis,
        color=tmp_key,
        projection=projection,
        size=size,
        legend_fontsize=legend_fontsize,
        add_outline=add_outline,
        legend_loc=legend_loc,
        frameon=frameon,
        palette=palette,
        return_fig=return_fig,
        show=show,
        **kwargs,
    )

    _ = adata.obsm.pop("X_actionet2d")
    _ = adata.obs.pop(tmp_key)

    return p


def visualize_markers(
    adata: AnnData,
    genes: Union[str, list, pd.Series, None] = None,
    color_map: Union[Colormap, str, None] = "YlOrRd",
    alpha: float = 0.85,
    **kwargs,
) -> Union[Figure, Axes, None]:

    feature_names = pd.Series([x.decode() if isinstance(x, (bytes, bytearray)) else x for x in list(adata.var.index)])
    adata.var.index = feature_names

    X = adata[:, genes].X
    if alpha != 0:
        X_smooth = diffusion(adata, X, return_raw=True, alpha_val=alpha)
    else:
        X_smooth = X

    if isinstance(X_smooth, np.ndarray):
        if isinstance(genes, str):
            p = plot_ACTIONet(adata, X_smooth[:, 0], color_map=color_map, title=genes, **kwargs)
        elif isinstance(genes, list):
            p = [plot_ACTIONet(adata, X_smooth[:, k], color_map=color_map, title=genes[k], **kwargs) for k in range(X_smooth.shape[1])]

    return p


def archetype_footprint(
    adata: AnnData,
    color_map: Union[Colormap, str, None] = "YlOrRd",
    **kwargs,
) -> Union[Figure, Axes, None]:

    X_smooth = adata.obsm["archetype_footprint"]

    p = [
        plot_ACTIONet(
            adata,
            X_smooth[:, k],
            color_map=color_map,
            title="Archetype %d" % (k + 1),
            **kwargs,
        )
        for k in range(X_smooth.shape[1])
    ]

    return p


def plot_ACTIONet_interactive(
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
    show_legend: Optional[bool] = None,
    hover_text: Union[list, pd.Series, np.ndarray] = None,
    plot_3d: Optional[bool] = False,
    point_order: Union[list, pd.Series, np.ndarray] = None,
    coordinate_key: Optional[str] = None,
    color_key: Optional[str] = "denovo_color",
) -> go.Figure:
    """Creates an interactive ACTIONet plot with plotly
    :param data:AnnData object with coordinates in '.obsm[coordinate_key]' or numeric matrix of X-Y(-Z) coordinates. \
        If data is AnnData, 'coordinate_key' defaults to 'ACTIONet3D' if 'plot_3d=True' or 'ACTIONet2D' if otherwise.
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
        If data is AnnData and 'coordinate_key=None', 'coordinate_key' defaults to 'ACTIONet3D'. \
        If data is matrix-like, it must have at least 3 columns.
    :param point_order: Numeric list=like object specifying order in which to plot individual points (default:None). \
        If None, points are plotted in random order.
    :param coordinate_key: If 'data' is AnnData, key of '.obsm' pertaining to plot coordinates.
    :param color_key:If data is AnnData, key of '.obsm' containing point-wise RGB color mappings (default:'denovo_color'). \
    Used only if no other color mapping parameters are given.
    ...

    :return plotly figure
    """

    if plot_3d:
        coor_dims = 3
        if isinstance(data, AnnData) and coordinate_key is None:
            coordinate_key = "ACTIONet3D"
    else:
        coor_dims = 2
        if isinstance(data, AnnData) and coordinate_key is None:
            coordinate_key = "ACTIONet2D"

    plot_coors = pu.get_plot_coors(data=data, coordinate_key=coordinate_key, scale_coors=True, coor_dims=coor_dims)
    plot_labels = pu.get_plot_labels(label_attr=label_attr, data=data)

    if plot_labels is None:
        plot_labels = pd.Series("NA", index=range(plot_coors.shape[0]), name="labels", dtype=str)
    else:
        plot_labels = pd.Series(plot_labels, name="labels", dtype=str)

    plot_data = pd.concat([plot_coors, plot_labels], axis=1)
    plot_data["idx"] = range(plot_data.shape[0])

    if hover_text is not None:
        plot_data["text"] = pd.Series(hover_text, dtype=str)
    else:
        if label_attr is None:
            plot_data["text"] = plot_data["idx"]
        else:
            plot_data["text"] = plot_data["labels"]

    if point_order is None:
        plot_data["pidx"] = sample(range(plot_data.shape[0]), plot_data.shape[0])
    else:
        plot_data["pidx"] = point_order

    if label_attr is None or any(elem is not None for elem in [color_attr, trans_attr]):
        # if hover_text is None:
        #     plot_data["text"] = plot_data["idx"]
        # else:
        #     plot_data["text"] = pd.Series(hover_text, dtype=str)

        plot_data["fill"] = pu.get_plot_colors(
            color_attr=color_attr,
            plot_labels=plot_labels,
            data=data,
            color_key=color_key,
            palette=palette,
            return_dict=False,
        )
        plot_data["color"] = [lighten_color(c, stroke_size) for c in plot_data["fill"]]

        plot_data["trans"] = pu.get_plot_transparency(
            trans_attr=trans_attr,
            adata=data,
            trans_fac=trans_fac,
            trans_th=trans_th,
            scale=True,
        )

        plot_data["fill"] = append_alpha_to_rgb(plot_data["fill"], plot_data["trans"], unzip_colors=True)
        plot_data["color"] = append_alpha_to_rgb(plot_data["color"], plot_data["trans"], unzip_colors=True)

        # if point_order is None:
        #     plot_data = plot_data.sample(frac=1).reset_index(drop=True)
        # else:
        #     plot_data["pidx"] = point_order
        #     plot_data = plot_data.sort_values(by="pidx").reset_index(drop=True)

        plot_data = plot_data.sort_values(by="pidx").reset_index(drop=True)

        if show_legend is None:
            show_legend = False

        p = pu.make_plotly_scatter_single_trace(
            x=plot_data["x"],
            y=plot_data["y"],
            z=plot_data["z"],
            label_attr=plot_data["labels"],
            cols_point=plot_data["fill"],
            cols_stroke=plot_data["color"],
            point_size=point_size,
            stroke_size=stroke_size,
            show_legend=show_legend,
            hover_text=plot_data["text"],
            plot_3d=plot_3d,
        )

    else:
        # if hover_text is None:
        #     plot_data["text"] = plot_data["labels"]
        # else:
        #     plot_data["text"] = pd.Series(hover_text, dtype=str)

        fill_dict = pu.get_plot_colors(
            color_attr=color_attr,
            plot_labels=plot_labels,
            data=data,
            color_key=color_key,
            palette=palette,
            return_dict=True,
        )

        stroke_dict = {k: lighten_color(v, stroke_contrast_fac) for (k, v) in fill_dict.items()}

        # if point_order is None:
        #     plot_data = plot_data.sample(frac=1).reset_index(drop=True)
        # else:
        #     plot_data["pidx"] = point_order
        #     plot_data = plot_data.sort_values(by="pidx").reset_index(drop=True)

        plot_data = plot_data.sort_values(by="pidx").reset_index(drop=True)

        if show_legend is None:
            show_legend = True

        p = pu.make_plotly_scatter_split_trace(
            x=plot_data["x"],
            y=plot_data["y"],
            z=plot_data["z"],
            label_attr=plot_data["labels"],
            fill_dict=fill_dict,
            stroke_dict=stroke_dict,
            show_legend=show_legend,
            hover_text=plot_data["text"],
            plot_3d=plot_3d,
        )

    return p
