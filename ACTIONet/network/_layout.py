from typing import Optional

import numpy as np
from anndata import AnnData

import _ACTIONet as _an
from ..tools import scale_matrix

def layout_network(
    adata: AnnData,
    scale: Optional[bool] = True,
    compactness_level: Optional[int] = 50,
    n_epochs: Optional[int] = 500,
    n_threads: Optional[int] = 8,
    copy: Optional[bool] = False
) -> Optional[AnnData]:
    """\
    Network layout

    Embedded the graph into 2D/3D space

    Parameters
    ----------
    adata:
        AnnData object storing the ACTIONet results
    compactness_level:
        Between 0-100
    n_epochs:
        Number of SGD epochs
    n_threads:
        Number of parallel threads
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
        adata : anndata.AnnData
        if `copy=True` returns None or else adds fields to `adata`:

        `.obsm['X_ACTIONred']`
        `.obsm['X_ACTIONet_2D']`
        `.obsm['X_ACTIONet_3D']`
        `.obsm['X_denovo_color']`

        `.uns['obsm_annot']['X_ACTIONred']`
        `.uns['obsm_annot']['X_ACTIONet_2D']`
        `.uns['obsm_annot']['X_ACTIONet_3D']`
        `.uns['obsm_annot']['X_denovo_color']`
    """
    if 'ACTIONet' not in adata.obsp.keys():
        raise ValueError(
            'Did not find adata.obsp[\'ACTIONet\']. '
            'Please run nt.build_network() first.'
        )

    adata = adata.copy() if copy else adata
    G = adata.obsp['ACTIONet']
    if scale:
        initial_coordinates = scale_matrix(adata.obsm['ACTION'])
    else:
        initial_coordinates = adata.obsm['ACTION']

    layout = _an.layout_ACTIONet(G, initial_coordinates.T, compactness_level, n_epochs, n_threads)

    adata.obsm['X_ACTIONred'] = initial_coordinates[:,:3]
    adata.obsm['X_ACTIONet2D'] = layout['coordinates']
    adata.obsm['X_ACTIONet3D'] = layout['coordinates_3D']
    adata.obsm['X_denovo_color'] = layout['colors']

    adata.uns.setdefault('obsm_annot', {}).update({
        'X_ACTIONred': {'type': np.array([b'embedding'], dtype=object)},
        'X_ACTIONet2D': {'type': np.array([b'embedding'], dtype=object)},
        'X_ACTIONet3D': {'type': np.array([b'embedding'], dtype=object)},
        'X_denovo_color': {'type': np.array([b'embedding'], dtype=object)},
    })
    return ACE if copy else None
