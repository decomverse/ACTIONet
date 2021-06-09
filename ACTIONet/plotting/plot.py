from typing import Optional
import plotly as pl
import plotly.io as pio
import plotly.graph_objs as go
from plotly.graph_objs.scatter import Line
from plotly.subplots import make_subplots
import plotly.express as px
import numpy as np
from adjustText import adjust_text
from anndata import AnnData
from scipy import stats
import pandas as pd
import seaborn as sns
from .color import adjust_lightness, hex_to_rgb, rgb_to_hex
from .palettes import palette_default
from .. import misc_utils as ut

pio.orca.config.use_xvfb = True
pio.orca.config.save()


def validate_plot_params(adata,coordinate_key,label_key,transparency_key):
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


def layout_labels(
    X: np.ndarray,
    Y: np.ndarray,
    fig:go.Figure,
    labels: list,
    fontsize: Optional[int] = 18,
    color: Optional[str] = "#0000FF",
    background: Optional[str] = "#FFFFFF") -> None:
    texts = []
    if isinstance(color, tuple):
        color = [color] * len(labels)
    for x, y, label, c in zip(X, Y, labels, color):
        fig.add_annotation(
            x=x,
            y=y,
            text=str(label),
            bgcolor=background,
            font=dict(
                family="sans serif",
                size=fontsize,
                color=color))


def plot_ACTIONet(
        adata: AnnData,
        label_key: Optional[str] = None,
        coordinate_key: Optional[str] = 'ACTIONet2D',
        transparency_key: Optional[str] = None,
        transparency_z_threshold: Optional[float] = -0.5,
        transparency_factor: Optional[float] = 1.5,
        border_contrast_factor: Optional[float] = 0.1,
        node_size: Optional[float] = None,
        add_text: Optional[bool] = True,
        palette: Optional[list] = None,
        output_file: Optional[str]=None):
    #validate supplied plot parameters
    validate_plot_params(adata,coordinate_key,label_key,transparency_key)
    coordinates = ut.scale_matrix(adata.obsm[coordinate_key])

    #get mapping of points to colors
    if label_key is None:
        v_col = [(r, g, b) for r, g, b in adata.obsm['denovo_color']]
    else:
        labels = adata.obs[label_key]
        unique_labels = sorted(np.unique(labels))
        if palette is None:
            if len(unique_labels) <= len(palette_default):
                palette = [hex_to_rgb(color) for color in palette_default]
            else:
                palette = [hex_to_rgb(color) for color in palette_default]
        elif isinstance(palette, str):
            #get plotly color palette from string specification
            assert palette.lower() in  px.colors.named_colorscales()
            if palette in px.colors.qualitative.__dict__.keys():
                palette=px.colors.qualitative.__dict__[palette]
            else:
                palette=px.colors.sequential.__dict__[palette]
        #make sure palette colors are rgb values in the 'rgb(r,g,b)' string format
        assert len(palette)>0
        if type(palette[0])==str:
            if palette[0].startswith('#'):
                palette=[hex_to_rgb(i) for i in palette]
        else:
            palette=[pl.colors.label_rgb(i) for i in palette]
        label_colors = {label: palette[i % len(palette)] for i, label in enumerate(unique_labels)}
        v_col = [label_colors[label] for label in labels]
    v_col_for_border=v_col
    # get darkened colors for plot marker outline
    v_col_darkened = [pl.colors.label_rgb(adjust_lightness(color, 1 - border_contrast_factor)) for color in v_col_for_border]

    #calculate transparency
    if transparency_key is not None:
        transparency = adata.obs[transparency_key].values
        z = (transparency - np.mean(transparency)) / np.std(transparency, ddof=1)
        betas = 1 / (1 + np.exp(- transparency_factor * (z - transparency_z_threshold)))
        betas[z > transparency_z_threshold] = 1
        betas **= transparency_factor
    else:
        betas=[1]*len(v_col_darkened)


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

    #plotly scatterplot
    if node_size is None:
        node_size=[1]*coordinates.shape[0]


    fig = make_subplots(rows=1, cols=1)
    fig.append_trace(go.Scatter(x=coordinates[:,0],
                                y=coordinates[:,1],
                                mode='markers',
                                marker=dict(opacity=betas,
                                            color=v_col,
                                            line=dict(width=1,
                                                      color=v_col_darkened))),row=1,col=1)
    fig.update_layout(yaxis_range=[y_min,y_max])
    fig.update_layout(xaxis_range=[x_min,x_max])
    #add cluster labels, if specified
    if add_text and label_key is not None:
        labels = adata.obs[label_key]
        unique_labels = sorted(np.unique(labels))

        colors = []
        centroids = np.zeros((len(unique_labels), coordinates.shape[1]))
        for i, label in enumerate(unique_labels):
            label_coordinates = coordinates[labels == label]
            centroids[i] = stats.trim_mean(label_coordinates, 0.2, axis=0)
            cur_text_color=pl.colors.label_rgb(adjust_lightness(palette[i % len(palette)], 0.5))
            colors.append(cur_text_color)
            fig.add_annotation(x=centroids[i,0],
                               y=centroids[i,1],
                               text=str(label),
                               bgcolor='#FFFFFF',
                               font=dict(
                                   family="sans serif",
                                   size=18,
                                   color=cur_text_color),row=1,col=1)
    #save to file if requested by user
    fig.update_layout(showlegend=False)
    if not(output_file is None):
        fig.write_image(output_file)
    #show the figure
    # fig.show()

    return fig


def plot_ACTIONet_interactive(adata: AnnData,
                              label_key: Optional[str]=None,
                              transparency_attribute:Optional[str]=None,
                              transparancy_z_threshold:Optional[float]=-1,
                              transparancey_factor=1,
                              node_size:Optional[float]=1,
                              palette:Optional[str]="CPal20",
                              enrichment_table:Optional[bool]=False,
                              top_features:Optional[int]=7,
                              blacklist_pattern:Optional[str]= "\\.|^RPL|^RPS|^MRP|^MT-|^MT|MALAT1|B2M|GAPDH",
                              title:Optional[str]="ACTIONet",
                              coordinate_slot:Optional[str]="ACTIONet2D",
                              threeD:bool=False):
    """
    Creates an interactive ACTIONet plot with plotly
    Parameters
    ---------
    ace:
        ACTIONet output object
    labels:
        Annotation of interest (clusters, celltypes, etc.) to be projected on the ACTIONet plot
    transparency_attr:
        Additional continuous attribute to project onto the transparency of nodes
    trans_z_threshold, trans_fact:
        Control the effect of transparency mapping
    node_size:
        Size of nodes in the ACTIONet plot
    palette:
        Color palette (named vector or a name for a given known palette)
    enrichment_table:
        To project the top-ranked features interactively.
    top_features:
        Number of features to show per cell
    blacklist_pattern:
        List of genes to filter-out
    title:
        Main title of the plot
    coordinate_slot:
        coordinate_slot in anndata object where visualization should be stored
    threeD:
        Whether to show the plot in 3D
    Returns
    -------
    Visualized ACTIONet
    """
    #validate supplied plot parameters
    validate_plot_params(adata,coordinate_key,label_key,transparency_key)

    #get the number of variables stored in obs matrix
    nV=adata.obs.columns.shape[0]

    #determine whether a 3D plot should be generated
    if((coordinate_slot=="ACTIONet2D") and (threeD is True)):
        coordinate_slot="ACTIONet3D"
    coordinates = ut.scale_matrix(adata.obsm[coordinate_key])


def plot_ACTIONet_gradient(adata: AnnData,
                           x: Optional[list] = None,
                           coordinate_key: Optional[str] = 'ACTIONet2D',
                           transparency_key: Optional[str] = None,
                           transparency_z_threshold: Optional[float] = -0.5,
                           transparancey_factor: Optional[float]=3,
                           alpha_val: Optional[float] = 0.85,
                           node_size: Optional[float] = 1,
                           add_text: Optional[bool] = True,
                           palette: Optional[list] = "Inferno",
                           title:Optional[str]="",
                           nonparametric:Optional[bool]=False,
                           output_file: Optional[str]=None):

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
    coordinates = ut.scale_matrix(adata.obsm[coordinate_key])


# def plot_ACTIONet_test(
#     data: Union[AnnData, pd.DataFrame, np.ndarray],
#     label_attr: Union[str, list, pd.Series, None] = None,
#     color_attr: Union[list, pd.Series, dict, None],
#
#
#     coordinate_attr: Optional[str] = "ACTIONet2D"
#     color_key: Optional[str] = "denovo_color",
#     palette: Union[str, list, pd.Series, np.ndarray] = palette_default
# ) -> :
#
#     something
#
#     return something else
