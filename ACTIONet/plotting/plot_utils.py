from typing import Optional, Union
import numpy as np
from scipy import sparse
import pandas as pd
from anndata import AnnData
from ._color import adjust_lightness, hex_to_rgb, rgb_to_hex
from ._palettes import palette_default
from .. import _misc_utils as ut

def get_plot_coors(
    X: Union[AnnData, np.ndarray, sparse.spmatrix],
    coordinate_attr: Optional[str],
    scale_coors: Optional[bool] = True,
    coor_dims: Optional[int] = 2
) -> pd.DataFrame:

    if isinstance(X, AnnData):
        coors = X.obsm[coordinate_attr]
    else:
        if not isinstance(X, (np.ndarray, sparse.spmatrix)):
            raise ValueError(
                "X must be AnnData or array"
            )

    coors = np.asarray(coors, dtype=np.float64)

    if scale_coors:
        coors = ut.scale_matrix(coors)

    coors = pd.DataFrame(coors)
    coors.columns = ['x', 'y', "z"][0:coor_dims]

    return coors

def get_plot_labels(
    label_attr: Union[str, list, pd.Series, np.ndarray, None],
    X: Union[AnnData, np.ndarray, sparse.spmatrix, None] = None,
) -> pd.Series:

    if label_attr is None:
        return None

    if isinstance(X, AnnData):
        plot_labels = ut._get_attr_or_split_idx(X, attr=label_attr, return_vec=True)
    else:
        plot_labels = label_attr

    plot_labels = pd.Series(plot_labels, dtype=str, name="labels")
    plot_labels = plot_labels.fillna("NA")

    return plot_labels

def get_plot_colors(
    color_attr: Union[str, list, pd.Series, pd.DataFrame, np.ndarray, None],
    plot_labels: Union[list, pd.Series, np.ndarray, None],
    X: Union[AnnData, np.ndarray, sparse.spmatrix, None] = None,
    color_key: Optional[str] = "denovo_color",
    palette: Union[str, list, pd.Series, np.ndarray] = palette_default
) -> pd.Series:

    n_dim = X.shape[0]

    if color_attr is not None:
        if isinstance(color_attr, (np.ndarray, pd.DataFrame)):

            if color_attr.shape[1] >= 3:
                plot_colors = [rgb_to_hex(color_attr[i, :]) for i in range(color_attr.shape[0])]
                plot_colors = pd.Series(plot_colors, dtype=str)
            else:
                something error

        elif isinstance(color_attr, (list, pd.Series)):
            if len(color_attr) == n_dim:
                # Test this shit
                plot_colors = pd.Series(color_attr, dtype=str)
            else:
                some error

        elif isinstance(color_attr, str):
            if len(color_attr) == n_dim:
                # Test this shit
                plot_colors = ut._get_attr_or_split_idx(X, attr=color_attr, return_vec=True)
            else:
                some error

        else:
            # err = sprintf("Invalid 'color_attr'.\n")
            # stop(err)
            something error

    elif plot_labels is not None:

        if len(color_attr) != n_dim:
            some error

        plot_colors = pd.Series(color_attr, dtype=str)
        plot_colors.fillna("NA")
        # label_names = pd.unique(plot_labels).sort

        # plot_labels = as.character(plot_labels)
        # plot_labels[ is.na(plot_labels)] = "NA"
        label_names = sort(unique(plot_labels))
        num_unique = length(label_names)

        if (num_unique == 1) {
            plot_colors =.default_colors(n_dim)
        } else {

        if (length(palette) == 1) {
            plot_palette = ggpubr::
                get_palette(palette, num_unique)
        } else if (length(palette) < num_unique) {
            plot_palette = CPal_default[1:num_unique]
            msg = sprintf("Not enough colors in 'palette'. Using default palette.\n")
            message(msg)
        } else {

        if (! is.null(names(palette))){
            if (all(label_names % in % names(palette))){
                plot_palette = palette[label_names]
            } else {
                plot_palette = palette[1:num_unique]
            }
        } else {
            plot_palette = palette[1:num_unique]
        }

        }

        names(plot_palette) = label_names
        plot_colors = plot_palette[match(plot_labels, names(plot_palette))]

        }

    else:

    return plot_colors
