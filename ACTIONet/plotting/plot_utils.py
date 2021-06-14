from typing import Optional, Union
import warnings
import numpy as np
from scipy import sparse
import pandas as pd
from anndata import AnnData
from .color import adjust_lightness, hex_to_rgb, rgb_to_hex
from .palettes import palette_default
from .. import misc_utils as ut
from matplotlib.colors import to_rgb, to_rgba
import plotly.graph_objects as go


def _default_colors(n) -> pd.Series:
    plot_colors = pd.Series("#FF6347", index=range(n), dtype=str)
    return plot_colors


def get_plot_coors(
        data: Union[AnnData, np.ndarray, sparse.spmatrix],
        coordinate_attr: Optional[str],
        scale_coors: Optional[bool] = True,
        coor_dims: Optional[int] = 2
) -> pd.DataFrame:
    if isinstance(data, AnnData):
        coors = data.obsm[coordinate_attr]
        if coors.shape[1] < coor_dims:
            err = "data in 'coordinate_attr' has < {dims} dimensions".format(keys=coor_dims)
            raise Exception(err)
    else:
        if not isinstance(data, (np.ndarray, sparse.spmatrix)):
            raise ValueError(
                "data must be AnnData, numpy.ndarray, or sparse.spmatrix"
            )
        coors = data

    coors = np.asarray(coors, dtype=np.float64)

    if scale_coors:
        coors = ut.scale_matrix(coors)

    coors = pd.DataFrame(coors)
    if coors.shape[1] < 3:
        coors.insert(2, "z", pd.NA)

    # coors.columns = ['x', 'y', "z"][0:coor_dims]
    coors.columns = ['x', 'y', "z"]
    coors.reset_index(drop=True, inplace=True)
    return coors


def get_plot_labels(
        label_attr: Union[str, list, pd.Series, np.ndarray, None],
        data: Union[AnnData, np.ndarray, sparse.spmatrix, None] = None,
) -> pd.Series:
    if label_attr is None:
        return None

    if isinstance(data, AnnData):
        plot_labels = ut.get_attr_or_split_idx(data, attr=label_attr, return_vec=True)
    else:
        plot_labels = label_attr

    plot_labels = pd.Series(plot_labels, dtype=str, name="labels")
    plot_labels = plot_labels.fillna("NA")
    plot_labels.reset_index(drop=True, inplace=True)
    return plot_labels


def get_plot_colors(
        color_attr: Union[str, list, pd.Series, pd.DataFrame, np.ndarray, None],
        plot_labels: Union[list, pd.Series, None],
        data: Union[AnnData, np.ndarray, sparse.spmatrix, None] = None,
        color_key: Optional[str] = "denovo_color",
        palette: Union[str, list, pd.Series, dict] = palette_default,
        return_dict: Optional[bool] = False
) -> pd.Series:
    if data is not None:
        n_dim = data.shape[0]
        no_data = False
    else:
        n_dim = None
        no_data = True

    if color_attr is not None:
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
            else:
                plot_colors = ut.get_attr_or_split_idx(data, attr=color_attr, return_vec=True)

        else:
            raise Exception("invalid color_attr")

    elif plot_labels is not None:

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
                plot_palette = palette[0:num_unique]

            palette_dict = dict(zip(label_names, plot_palette))

            if return_dict:
                plot_colors = palette_dict
            else:
                plot_colors = pd.Series(pd.NA, index=range(n_dim), dtype=str)
                for it in palette_dict.items():
                    plot_colors.loc[plot_labels == it[0]] = it[1]

    else:
        if no_data:
            raise Exception("'data' required for given parameters")

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
        scale: Optional[bool] = True
) -> pd.Series:
    if trans_attr is None:
        return 1

    data_is_AnnData = isinstance(adata, AnnData)
    if data_is_AnnData:
        n_dim = adata.shape[0]

    if isinstance(trans_attr, str):
        if not data_is_AnnData:
            raise Exception("'adata' must be AnnData if 'trans_attr' is str")
        alpha_fac = ut.get_attr_or_split_idx(n_dim, attr=trans_attr, return_vec=True)
    else:
        alpha_fac = pd.Series(trans_attr, dtype=float)

    if scale:
        z = ut.scale_matrix(alpha_fac)
    else:
        z = alpha_fac

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
        hover_text: Union[list, pd.Series, np.ndarray] = None,
        plot_3d: Optional[bool] = False,
) -> go.Figure:

    if hover_text is None:
        if label_attr is None:
            hover_text = list(range(len(x)))
        else:
            hover_text = label_attr

    axis_params = dict(
        showgrid=False,
        zeroline=False,
        visible=False
    )

    if plot_3d:
        p = go.Figure(
            data=go.Scatter3d(
                x=x,
                y=y,
                z=z,
                marker=dict(
                    color=cols_point,
                    size=point_size,
                    line=dict(
                        color=cols_stroke,
                        width=stroke_size
                    )
                ),
                text=hover_text,
                hoverinfo="text",
                mode='markers',
            ),
            layout=dict(
                scene=dict(
                    xaxis=axis_params,
                    yaxis=axis_params,
                    zaxis=axis_params
                ),
                showlegend=False
            )
        )
    else:
        p = go.Figure(
            data=go.Scattergl(
                x=x,
                y=y,
                marker=dict(
                    color=cols_point,
                    size=point_size,
                    line=dict(
                        color=cols_stroke,
                        width=stroke_size
                    )
                ),
                text=hover_text,
                hoverinfo="text",
                mode='markers',
            ),
            layout=dict(
                xaxis=axis_params,
                yaxis=axis_params,
                showlegend=False,
                paper_bgcolor='white',
                plot_bgcolor='white'
            )
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
        plot_3d: Optional[bool] = False
) -> go.Figure:

    plot_data = pd.DataFrame(dict(
        x=x,
        y=y,
        z=z,
        labels=label_attr,
    ))

    if hover_text is None:
        plot_data["text"] = plot_data["labels"]
    else:
        plot_data["text"] = pd.Series(hover_text, dtype=str)

    trace_names = sorted(plot_data["labels"].unique())

    axis_params = dict(
        showgrid=False,
        zeroline=False,
        visible=False
    )

    if plot_3d:
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
                        line=dict(
                            color=stroke_dict[n],
                            width=stroke_size
                        )
                    ),
                    text=sub_data["text"],
                    hoverinfo="text",
                    mode="markers",
                    name=n
                )
            )

        p.update_layout(
            scene=dict(
                xaxis=axis_params,
                yaxis=axis_params,
                zaxis=axis_params
            ),
            showlegend=show_legend
        )

    else:
        p = go.Figure()

        for n in trace_names:
            sub_data = plot_data[plot_data["labels"] == n]
            p.add_trace(
                go.Scattergl(
                    x=sub_data["x"],
                    y=sub_data["y"],
                    marker=dict(
                        color=fill_dict[n],
                        size=point_size,
                        line=dict(
                            color=stroke_dict[n],
                            width=stroke_size
                        )
                    ),
                    text=sub_data["text"],
                    hoverinfo="text",
                    mode="markers",
                    name=n
                )
            )

        p.update_layout(
            xaxis=axis_params,
            yaxis=axis_params,
            showlegend=show_legend,
            paper_bgcolor='white',
            plot_bgcolor='white'
        )

    return p
