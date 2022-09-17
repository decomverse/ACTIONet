import warnings
from typing import Optional, Union

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from anndata import AnnData
from matplotlib.colors import to_rgb
from scipy import sparse

from ACTIONet.plotting.color import rgb_to_hex
from ACTIONet.plotting.palettes import palette_default
from ACTIONet.tools.utils_public import scale_matrix


def _default_colors(n) -> pd.Series:
    plot_colors = pd.Series("#FF6347", index=range(n), dtype=str)
    return plot_colors


def get_plot_coors(
    data: Union[AnnData, np.ndarray, sparse.spmatrix],
    coordinate_key: Optional[str],
    scale_coors: Optional[bool] = True,
    coor_dims: Optional[int] = 2,
) -> pd.DataFrame:
    if isinstance(data, AnnData):
        coors = data.obsm[coordinate_key]
        if coors.shape[1] < coor_dims:
            err = "data in 'coordinate_key' has < {dims} dimensions".format(dims=coor_dims)
            raise Exception(err)
    else:
        if not isinstance(data, (np.ndarray, sparse.spmatrix)):
            raise ValueError("data must be AnnData, numpy.ndarray, or sparse.spmatrix")
        coors = data

    coors = np.asarray(coors, dtype=np.float64)

    if scale_coors:
        coors = scale_matrix(coors)

    coors = pd.DataFrame(coors)
    if coors.shape[1] < 3:
        coors.insert(2, "z", pd.NA)

    # coors.columns = ['x', 'y', "z"][0:coor_dims]
    coors.columns = ["x", "y", "z"]
    coors.reset_index(drop=True, inplace=True)
    return coors


def get_plot_labels(
    label_attr: Union[str, list, pd.Series, np.ndarray, None],
    data: Union[AnnData, np.ndarray, sparse.spmatrix, None] = None,
) -> Union[pd.Series, None]:
    if label_attr is None:
        return None

    if isinstance(data, AnnData) and isinstance(label_attr, str):
        plot_labels = data.obs[label_attr]
    elif isinstance(label_attr, list) or isinstance(label_attr, np.ndarray) or isinstance(label_attr, pd.Series):
        plot_labels = pd.Series(label_attr, dtype=str, name="labels")

    plot_labels = plot_labels.fillna("NA")
    plot_labels.reset_index(drop=True, inplace=True)

    return plot_labels


def get_plot_colors_from_plot_labels(plot_labels: Optional[Union[list, pd.Series]], n_dim: int, palette: Union[str, list, pd.Series, dict], return_dict: bool):
    plot_labels = pd.Series(plot_labels, dtype=str)
    plot_labels.fillna("NA")

    label_names = sorted(plot_labels.unique())
    num_unique = len(label_names)

    if num_unique == 1:
        plot_colors = _default_colors(n_dim)

    elif isinstance(palette, dict):
        if not all(k in label_names for k in palette.keys()):
            missing_keys = list(set(palette.keys()) - set(label_names))
            err = "keys for {keys} missing from palette".format(keys=missing_keys)
            raise Exception(err)

        if return_dict:
            plot_colors = palette
        else:
            plot_colors = pd.Series(pd.NA, index=range(n_dim), dtype=str)
            for it in palette.items():
                plot_colors.loc[label_names == it[0]] = it[1]

    else:
        if len(palette) < num_unique:
            warnings.warn("Not enough colors in 'palette'. Using default palette.")
            plot_palette = palette_default[0:num_unique]
        else:
            plot_palette = list(palette)[0:num_unique]

        palette_dict = dict(zip(label_names, plot_palette))

        if return_dict:
            plot_colors = palette_dict
        else:
            plot_colors = pd.Series(pd.NA, index=range(n_dim), dtype=str)
            for it in palette_dict.items():
                plot_colors.loc[plot_labels == it[0]] = it[1]
    return plot_colors


def get_plot_colors_from_color_attr(color_attr: Union[str, list, pd.Series, pd.DataFrame, np.ndarray], data: Optional[Union[AnnData, np.ndarray, sparse.spmatrix]], no_data: bool, n_dim: Optional[int]):
    if isinstance(color_attr, (np.ndarray, pd.DataFrame)):
        if color_attr.shape[1] >= 3:
            plot_colors = [rgb_to_hex(color_attr[i, :]) for i in range(color_attr.shape[0])]
            plot_colors = pd.Series(plot_colors, dtype=str)
        else:
            raise Exception("invalid color_attr")

    elif isinstance(color_attr, (list, pd.Series)):
        if no_data:
            raise Exception("'data' must not be None if 'color_attr' is str")
        elif len(color_attr) != n_dim:
            raise Exception("length of 'color_attr' must match 'data'")
        else:
            plot_colors = pd.Series(color_attr, dtype=str)

    elif isinstance(color_attr, str):
        if no_data:
            raise Exception("'data' must not be None if 'color_attr' is str")
        elif isinstance(data, AnnData) and isinstance(color_attr, str):
            plot_colors = data.obs[color_attr]
    else:
        raise Exception("invalid color_attr")
    return plot_colors


def get_plot_colors(
    color_attr: Optional[Union[str, list, pd.Series, pd.DataFrame, np.ndarray]],
    plot_labels: Optional[Union[list, pd.Series]],
    data: Optional[Union[AnnData, np.ndarray, sparse.spmatrix]],
    color_key: Optional[str] = "denovo_color",
    palette: Union[str, list, pd.Series, dict] = palette_default,
    return_dict: Optional[bool] = False,
) -> pd.Series:
    if data is not None:
        n_dim = data.shape[0]
        no_data = False
    else:
        n_dim = None
        no_data = True

    if color_attr is not None:
        return get_plot_colors_from_color_attr(color_attr=color_attr, data=data, no_data=no_data, n_dim=n_dim)
    elif plot_labels is not None:
        return get_plot_colors_from_plot_labels(plot_labels=plot_labels, n_dim=n_dim, palette=palette, return_dict=bool(return_dict))
    else:
        if no_data:
            raise Exception("'data' required for given parameters")
        if isinstance(data, AnnData):
            if color_key is None or color_key not in data.obsm_keys():
                plot_colors = _default_colors(n_dim)
            else:
                plot_colors = pd.Series([tuple(r) for r in data.obsm[color_key]])

    if not return_dict:
        plot_colors = [to_rgb(f) for f in plot_colors]
    return plot_colors


def get_plot_transparency(
    trans_attr: Union[str, list, pd.Series, np.ndarray, None] = None,
    adata: Union[AnnData, None] = None,
    trans_fac: Optional[float] = 1.5,
    trans_th: Optional[float] = -0.5,
    scale: Optional[bool] = True,
) -> Union[pd.Series, int]:
    if trans_attr is None:
        return 1

    if isinstance(trans_attr, str) and not isinstance(adata, AnnData):
        raise Exception("'adata' must be AnnData if 'trans_attr' is str")

    if isinstance(adata, AnnData) and isinstance(trans_attr, str):
        alpha_fac = adata.obs[trans_attr]
        alpha_fac = pd.Series(alpha_fac, dtype=float)

    if scale:
        z = scale_matrix(alpha_fac)
    else:
        z = alpha_fac

    if isinstance(trans_fac, float):
        alpha_val = 1 / (1 + np.exp(-trans_fac * (z - trans_th)))
        alpha_val = np.where(z > trans_th, 1, alpha_val)
        alpha_val = np.power(alpha_val, trans_fac)

    return alpha_val


def make_plotly_scatter_single_trace(
    x: Union[list, np.ndarray, pd.Series],
    y: Union[list, np.ndarray, pd.Series],
    z: Union[list, np.ndarray, pd.Series, None] = None,
    label_attr: Union[list, np.ndarray, pd.Series, None] = None,
    cols_point: Union[str, list, np.ndarray, pd.Series, None] = None,
    cols_stroke: Union[str, list, np.ndarray, pd.Series, None] = None,
    point_size: Optional[float] = 3,
    stroke_size: Optional[float] = 0.3,
    show_legend: Optional[bool] = False,
    hover_text: Union[list, pd.Series, np.ndarray] = None,
    plot_3d: Optional[bool] = False,
) -> go.Figure:
    if hover_text is None:
        if label_attr is None:
            hover_text = list(range(len(x)))
        else:
            hover_text = label_attr

    axis_params = dict(showgrid=False, zeroline=False, visible=False)

    if plot_3d:
        p = go.Figure(
            data=go.Scatter3d(
                x=x,
                y=y,
                z=z,
                marker=dict(
                    color=cols_point,
                    size=point_size,
                    line=dict(color=cols_stroke, width=stroke_size),
                ),
                text=hover_text,
                hoverinfo="text",
                mode="markers",
            ),
            layout=dict(
                scene=dict(xaxis=axis_params, yaxis=axis_params, zaxis=axis_params),
                showlegend=show_legend,
            ),
        )
    else:
        p = go.Figure(
            data=go.Scattergl(
                x=x,
                y=y,
                marker=dict(
                    color=cols_point,
                    size=point_size,
                    line=dict(color=cols_stroke, width=stroke_size),
                ),
                text=hover_text,
                hoverinfo="text",
                mode="markers",
            ),
            layout=dict(
                xaxis=axis_params,
                yaxis=axis_params,
                showlegend=show_legend,
                paper_bgcolor="white",
                plot_bgcolor="white",
            ),
        )

    return p


def make_plotly_scatter_split_trace(
    x: Union[list, np.ndarray, pd.Series],
    y: Union[list, np.ndarray, pd.Series],
    z: Union[list, np.ndarray, pd.Series, None] = None,
    label_attr: Union[list, np.ndarray, pd.Series, None] = None,
    fill_dict: Union[str, list, np.ndarray, pd.Series, None] = None,
    stroke_dict: Union[str, list, np.ndarray, pd.Series, None] = None,
    point_size: Optional[float] = 3,
    stroke_size: Optional[float] = 0.3,
    show_legend: Optional[bool] = True,
    hover_text: Union[list, pd.Series, np.ndarray] = None,
    plot_3d: Optional[bool] = False,
) -> go.Figure:
    plot_data = pd.DataFrame(
        dict(
            x=x,
            y=y,
            z=z,
            labels=label_attr,
        )
    )

    if hover_text is None:
        plot_data["text"] = plot_data["labels"]
    else:
        plot_data["text"] = pd.Series(hover_text, dtype=str)

    trace_names = sorted(plot_data["labels"].unique())

    axis_params = dict(showgrid=False, zeroline=False, visible=False)

    if plot_3d and (isinstance(fill_dict, np.ndarray) or isinstance(fill_dict, list) or isinstance(fill_dict, pd.Series)) and (isinstance(stroke_dict, np.ndarray) or isinstance(stroke_dict, list) or isinstance(stroke_dict, pd.Series)):
        p = go.Figure()

        for n in trace_names:
            sub_data = plot_data[plot_data["labels"] == n]
            p.add_trace(
                go.Scatter3d(
                    x=sub_data["x"],
                    y=sub_data["y"],
                    z=sub_data["z"],
                    marker=dict(
                        color=fill_dict[n],
                        size=point_size,
                        line=dict(color=stroke_dict[n], width=stroke_size),
                    ),
                    text=sub_data["text"],
                    hoverinfo="text",
                    mode="markers",
                    name=n,
                )
            )

        p.update_layout(
            scene=dict(xaxis=axis_params, yaxis=axis_params, zaxis=axis_params),
            showlegend=show_legend,
        )

    else:
        p = go.Figure()

        if (isinstance(fill_dict, np.ndarray) or isinstance(fill_dict, list) or isinstance(fill_dict, pd.Series)) and (isinstance(stroke_dict, np.ndarray) or isinstance(stroke_dict, list) or isinstance(stroke_dict, pd.Series)):
            for n in trace_names:
                sub_data = plot_data[plot_data["labels"] == n]
                p.add_trace(
                    go.Scattergl(
                        x=sub_data["x"],
                        y=sub_data["y"],
                        marker=dict(
                            color=fill_dict[n],
                            size=point_size,
                            line=dict(color=stroke_dict[n], width=stroke_size),
                        ),
                        text=sub_data["text"],
                        hoverinfo="text",
                        mode="markers",
                        name=n,
                    )
                )

        p.update_layout(
            xaxis=axis_params,
            yaxis=axis_params,
            showlegend=show_legend,
            paper_bgcolor="white",
            plot_bgcolor="white",
        )

    return p
