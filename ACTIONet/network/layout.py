from typing import Optional, Union

import numpy as np
from scipy import sparse
from anndata import AnnData

import _ACTIONet as _an
from .. import misc_utils as ut


def layout_network(
    adata: Optional[AnnData] = None,
    G: Union[np.ndarray, sparse.spmatrix] = None,
    initial_coordinates: Union[np.ndarray, sparse.spmatrix] = None,
    reduction_key: Optional[str] = "ACTION",
    net_key: Optional[str] = "ACTIONet",
    compactness_level: Optional[int] = 50,
    n_epochs: Optional[int] = 1000,
    layout_algorithm: Optional[str] = "tumap",
    thread_no: Optional[int] = 8,
    seed: Optional[int] = 0,
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False,
) -> Union[AnnData, dict, None]:

    """Network layout, Embeds the graph into 2D/3D space
    :param adata: AnnData object possibly containing 'reduction_key' in '.obsm' and 'net_key' in '.obsp'.
    :param G:Adjacency matrix to use for constructing layout. Required if 'adata=None'.
    :param initial_coordinates: Reduced representation matrix to use for constructing layout. Required if 'adata=None'.
    :param reduction_key: Key of 'adata.obms' containing reduced matrix to use for 'initial_coordinates' in 'layoutNetwork()' (default="ACTION").Ignored if 'adata=None'.
    :param net_key: Key of 'adata.obmp' containing adjacency matrix to use for 'G' in 'layoutNetwork()' (default="ACTIONet").  Ignored if 'adata=None'.
    :param compactness_level: Between 0-100. Ignored if 'layout_algorithm="tumap"'.
    :param layout_algorithm: Algorithm to use for constructing layout: \
    `tumap` (default) \
    t-distributed UMAP \
    `umap` Modified UMAP
    :param: n_epochs: Number of SGD epochs.
    :param thread_no: Number of threads. Defaults to number of threads available - 2.
    :param seed: Random seed
    :param copy: If 'adata' is given, return a copy instead of writing to `adata`
    :param return_raw: If 'adata' is given, return dict of raw 'layoutNetwork()' output instead of storing to 'adata'.
    ...
    :return adata: anndata.AnnData \
    if 'adata' given and `copy=True` returns None or else adds fields to `adata`: \
    `.obsp[net_key_out]`
    :return layout : dict \
    If 'adata=None' or 'return_raw=True', returns dict with 2D/3D coordinates and color values.
    """

    alg_name = layout_algorithm.upper()
    if alg_name not in ["UMAP", "TUMAP"]:
        raise ValueError("'layout_algorithm' must be 'tumap' or 'umap'.")

    if adata is not None:
        if isinstance(adata, AnnData):
            adata = adata.copy() if copy else adata
            initial_coordinates = (
                initial_coordinates
                if initial_coordinates is not None
                else adata.obsm[reduction_key]
            )
            G = G if G is not None else adata.obsp[net_key]
        else:
            raise ValueError("'adata' is not an AnnData object.")
    else:
        if G is None or initial_coordinates is None:
            raise ValueError(
                "'G' and 'initial_coordinates' cannot be NoneType if 'adata=None'."
            )
        if not isinstance(G, (np.ndarray, sparse.spmatrix)):
            raise ValueError("'G' must be numpy.ndarray or sparse.spmatrix.")
        if not isinstance(initial_coordinates, (np.ndarray, sparse.spmatrix)):
            raise ValueError(
                "'initial_coordinates' must be numpy.ndarray or sparse.spmatrix."
            )

    G = G.astype(dtype=np.float64)
    initial_coordinates = ut.scale_matrix(initial_coordinates).T.astype(
        dtype=np.float64
    )

    layout = _an.layoutNetwork(
        G=G,
        initial_position=initial_coordinates,
        algorithm=alg_name,
        compactness_level=compactness_level,
        n_epochs=n_epochs,
        thread_no=thread_no,
        seed=seed,
    )

    if return_raw or adata is None:
        return layout
    else:
        adata.obsm["ACTIONred"] = initial_coordinates[0:3, :].T
        adata.obsm["ACTIONet2D"] = layout["coordinates"]
        adata.obsm["ACTIONet3D"] = layout["coordinates_3D"]
        adata.obsm["denovo_color"] = layout["colors"]
        adata.uns.setdefault("obsm_annot", {}).update(
            {
                "ACTIONred": {"type": np.array([b"embedding"], dtype=object)},
                "ACTIONet2D": {"type": np.array([b"embedding"], dtype=object)},
                "ACTIONet3D": {"type": np.array([b"embedding"], dtype=object)},
                "denovo_color": {"type": np.array([b"embedding"], dtype=object)},
            }
        )
        return adata if copy else None
