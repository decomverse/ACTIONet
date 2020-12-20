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
    copy: Optional[bool] = False,
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

        `.obsm['X_ACTIONet_2D']`
        `.obsm['X_ACTIONet_3D']`
        `.uns['ACTIONet']['colors']`
    """
    if "ACTIONet" not in adata.obsp.keys():
        raise ValueError(
            "Did not find adata.obsp['ACTIONet']. "
            "Please run nt.build_network() first."
        )

    adata = adata.copy() if copy else adata
    G = adata.obsp["ACTIONet"]
    if scale:
        S_r = scale_matrix(adata.obsm["ACTION_S_r"]).T
    else:
        S_r = adata.obsm["ACTION_S_r"].T

    layout = _an.layout_ACTIONet(G, S_r, compactness_level, n_epochs, n_threads)

    adata.obsm["X_ACTIONet_2D"] = layout["coordinates"]
    adata.obsm["X_ACTIONet_3D"] = layout["coordinates_3D"]
    adata.uns["ACTIONet"].update({"colors": layout["colors"]})

    return ACE if copy else None
