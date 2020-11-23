from typing import Optional

import matplotlib as mpl
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
from anndata import AnnData
from scipy import stats

from ._color import adjust_lightness, hex_to_rgba
from .palettes import palette_20, palette_88
from ..tools import scale_matrix

def layout_labels(
    X: np.ndarray,
    Y: np.ndarray,
    labels: list,
    ax: mpl.axes.Axes,
    fontsize: Optional[int] = 1,
    color: Optional[tuple] = (0., 0., 0., 1.),
    background: Optional[tuple] = (1., 1., 1., 1.),
) -> None:
    texts = []
    if isinstance(color, tuple):
        color = [color] * len(labels)
    for x, y, label, c in zip(X, Y, labels, color):
        text = ax.text(x, y, label, color=c)
        text.set_path_effects([
            path_effects.Stroke(linewidth=5, foreground=background),
            path_effects.Normal()
        ])
        texts.append(text)
    adjust_text(texts)

def plot_ACTIONet(
    adata: AnnData,
    label_key: Optional[str] = None,
    coordinate_key: Optional[str] = 'X_ACTIONet2D',
    transparency_key: Optional[str] = None,
    transparency_z_threshold: Optional[float] = -0.5,
    transparency_factor: Optional[float] = 1.5,
    border_contrast_factor: Optional[float] = 0.1,
    node_size: Optional[float] = 0.1,
    add_text: Optional[bool] = True,
    palette: Optional[list] = None,
    ax: Optional[mpl.axes.Axes] = None,
) -> Optional[mpl.axes.Axes]:
    if coordinate_key not in adata.obsm.keys():
        raise ValueError(
            f'Did not find adata.obsm[\'{coordinate_key}\']. '
            'Please run nt.layout_ACTIONet() first.'
        )
    if label_key is not None and label_key not in adata.obs.columns:
        raise ValueError(f'Did not find adata.obs[\'{label_key}\'].')
    if transparency_key is not None and transparency_key not in adata.obs.columns:
        raise ValueError(f'Did not find adata.obs[\'{transparency_key}\'].')

    coordinates = scale_matrix(adata.obsm[coordinate_key])

    if label_key is None:
        v_col = [(r, g, b, 1.0) for r, g, b in adata.obsm['X_denovo_color']]
    else:
        labels = adata.obs[label_key]
        unique_labels = sorted(np.unique(labels))
        if palette is None:
            if len(unique_labels) <= len(palette_20):
                palette = [hex_to_rgba(color) for color in palette_20]
            else:
                palette = [hex_to_rgba(color) for color in palette_88]
        elif isinstance(palette, str):
            cmap = plt.get_cmap(palette)
            palette = [cmap(i) for i in range(cmap.N)]

        label_colors = {label: palette[i % len(palette)] for i, label in enumerate(unique_labels)}
        v_col = [label_colors[label] for label in labels]

    # Transparency
    v_col_darkened = [adjust_lightness(color, 1 - border_contrast_factor) for color in v_col]
    if transparency_key is not None:
        transparency = adata.obs[transparency_key].values
        z = (transparency - np.mean(transparency)) / np.std(transparency, ddof=1)
        betas = 1 / (1 + np.exp(- transparency_factor * (z - transparency_z_threshold)))
        betas[z > transparency_z_threshold] = 1
        betas **= transparency_factor
        v_col_border = [(r, g, b, beta) for (r, g, b, a), beta in zip(v_col_darkened, betas)]
    else:
        v_col_border = v_col_darkened

    x = coordinates[:,0]
    y = coordinates[:,1]
    x_min = np.min(x)
    x_max = np.max(x)
    y_min = np.min(y)
    y_max = np.max(y)

    x_min = x_min - (x_max - x_min) / 20
    x_max = x_max + (x_max - x_min) / 20
    y_min = y_min - (y_max - y_min) / 20
    y_max = y_max + (y_max - y_min) / 20

    permutation = np.random.permutation(adata.shape[0])
    _ax = ax if ax is not None else plt.subplots()[1]
    _ax.scatter(
        x[permutation],
        y[permutation],
        c=np.array(v_col)[permutation],
        s=node_size,
        edgecolors=np.array(v_col_border)[permutation],
    )

    if add_text and label_key is not None:
        labels = adata.obs[label_key]
        unique_labels = sorted(np.unique(labels))

        colors = []
        centroids = np.zeros((len(unique_labels), coordinates.shape[1]))
        for i, label in enumerate(unique_labels):
            label_coordinates = coordinates[labels == label]
            centroids[i] = stats.trim_mean(label_coordinates, 0.2, axis=0)

            colors.append(adjust_lightness(palette[i % len(palette)], 0.5))
        layout_labels(
            centroids[:,0],
            centroids[:,1],
            unique_labels,
            ax=_ax,
            color=colors,
            background=hex_to_rgba('#eeeeee')
        )

    _ax.set_xlim(x_min, x_max)
    _ax.set_ylim(y_min, y_max)
    _ax.set_axis_off()
    plt.show()
