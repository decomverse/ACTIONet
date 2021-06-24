from typing import Optional, Union

import numpy as np
from scipy import sparse
from anndata import AnnData

import _ACTIONet as _an
from .. import misc_utils as ut


def layout_network(
        adata: Optional[AnnData] = None,
        G: Union[np.ndarray, sparse.spmatrix] = None,
        S_r: Union[np.ndarray, sparse.spmatrix] = None,
        reduction_name: Optional[str] = "ACTION",
        net_name: Optional[str] = "ACTIONet",
        compactness_level: Optional[int] = 50,
        n_epochs: Optional[int] = 1000,
        layout_alg: Optional[int] = 0,
        thread_no: Optional[int] = 8,
        seed: Optional[int] = 0,
        copy: Optional[bool] = False,
        return_raw: Optional[bool] = False
) -> Union[AnnData, dict, None]:
    """Network layout, Embeds the graph into 2D/3D space
    :param adata: AnnData object possibly containing 'reduction_name' in '.obsm' and 'net_name' in '.obsp'.
    :param G:Adjacency matrix to use for constructing layout. Required if 'adata=None'.
    :param S_r: Reduced representation matrix to use for constructing layout. Required if 'adata=None'.
    :param reduction_name:Key of 'adata.obms' containing reduced matrix to use for 'S_r' in 'layout_ACTIONet()' (default="ACTION").Ignored if 'adata=None'.
    :param net_name: Key of 'adata.obmp' containing adjacency matrix to use for 'G' in 'layout_ACTIONet()' (default="ACTIONet").  Ignored if 'adata=None'.
    :param compactness_level: Between 0-100. Ignored if 'layout_alg=0'.
    :param layout_alg:Algorithm to use for constructing layout: \
    `0` (the default) \
    t-distributed UMAP \
    `1` Modified UMAP
    :param: n_epochs: Number of SGD epochs.
    :param thread_no: Number of threads. Defaults to number of threads available - 2.
    :param seed: Random seed
    :param copy: If 'adata' is given, return a copy instead of writing to `adata`
    :param return_raw: If 'adata' is given, return dict of raw 'layout_ACTIONet()' output instead of storing to 'adata'.
    ...
    :return adata: anndata.AnnData \
    if 'adata' given and `copy=True` returns None or else adds fields to `adata`: \
    `.obsp[net_name_out]`
    :return layout : dict \
    If 'adata=None' or 'return_raw=True', returns dict with 2D/3D coordinates and color values.
    """

    if adata is not None:
        if isinstance(adata, AnnData):
            adata = adata.copy() if copy else adata
            S_r = S_r if S_r is not None else adata.obsm[reduction_name]
            G = G if G is not None else adata.obsp[net_name]
        else:
            raise ValueError("'adata' is not an AnnData object.")
    else:
        if G is None or S_r is None:
            raise ValueError("'G' and 'S_r' cannot be NoneType if 'adata=None'.")
        if not isinstance(G, (np.ndarray, sparse.spmatrix)):
            raise ValueError("'G' must be numpy.ndarray or sparse.spmatrix.")
        if not isinstance(S_r, (np.ndarray, sparse.spmatrix)):
            raise ValueError("'S_r' must be numpy.ndarray or sparse.spmatrix.")

    G = G.astype(dtype=np.float64)
    S_r = ut.scale_matrix(S_r).T.astype(dtype=np.float64)
    layout = _an.layout_ACTIONet(G, S_r, compactness_level, n_epochs, layout_alg, thread_no, seed)

    if return_raw or adata is None:
        return layout
    else:
        adata.obsm["ACTIONred"] = S_r[0:3, :].T
        adata.obsm["ACTIONet2D"] = layout["coordinates"]
        adata.obsm["ACTIONet3D"] = layout["coordinates_3D"]
        adata.obsm["denovo_color"] = layout["colors"]
        adata.uns.setdefault("obsm_annot", {}).update(
            {
                "ACTIONred": {"type": np.array([b'embedding'], dtype=object)},
                "ACTIONet2D": {"type": np.array([b'embedding'], dtype=object)},
                "ACTIONet3D": {"type": np.array([b'embedding'], dtype=object)},
                "denovo_color": {"type": np.array([b'embedding'], dtype=object)},
            })
        return adata if copy else None
