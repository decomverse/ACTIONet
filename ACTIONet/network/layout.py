from typing import Optional, Union

import numpy as np
from anndata import AnnData
from scipy import sparse

import _ACTIONet as _an

from ..tools import utils_public as ut


def layout(
    data: Union[AnnData, np.ndarray, sparse.spmatrix],
    algorithm: Optional[str] = "umap",
    initial_coordinates: Union[np.ndarray, sparse.spmatrix] = None,
    presmooth_network: Optional[bool] = False,
    sim2dist: Optional[int] = 2,
    min_dist: Optional[float] = 1.0,
    spread: Optional[float] = 1.0,
    learning_rate: Optional[float] = 1.0,
    n_epochs: Optional[int] = 250,
    seed: Optional[int] = 0,
    reduction_key: Optional[str] = "ACTION",
    net_key: Optional[str] = "ACTIONet",
    thread_no: Optional[int] = 0,
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False,
) -> Union[AnnData, sparse.spmatrix, None]:
    """Computes knn/k*nn graphs from input data

    Parameters
    ----------
    data : Union[AnnData, np.ndarray, sparse.spmatrix]
        Adjacency matrix of the input graph or AnnData object containing the network.
    algorithm : str, optional
        one of k*nn, by default "k*nn"
    initial_coordinates: Union[np.ndarray, sparse.spmatrix]
        Coordinates to initialize the layout.
    compactness_level: int, optional
        Compactness of the layout (if "umap" is used as the algorithm), by default 50.
    n_epochs: int, optional
        Number of SGD epochs, by default 1000
    seed: int, optional
        Random seed, by default 0.
    reduction_key: str, optional
        if data is an `AnnData` object, in which `obsm` slot to find initial_coordiates, by default "ACTION".
    net_key: str, optional
        if data is an `AnnData` object, in which `obsp` slot to find iput graph, by default "ACTIONet".
    thread_no : Optional[int], optional
        Number of threads, by default 0
    copy : Optional[bool], optional
        If 'data' is AnnData, return a copy instead of writing to `data`, by default False
    return_raw : Optional[bool], optional
        If `return_raw=True` and `data` is AnnData, return sparse adjacency matrix directly, by default False

    Returns
    -------
    Union[AnnData, sparse.spmatrix, None]
        adata: anndata.AnnData if 'adata' given and `copy=True` returns None or else adds fields to `adata`
        G :scipy.sparse.spmatrix. Sparse adjacency matrix encoding ACTIONet if 'return_raw=True'
    """

    alg_name = str(algorithm).lower()
    if alg_name not in ["umap", "tumap"]:
        raise ValueError("'layout_algorithm' must be 'tumap' or 'umap'.")

    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        initial_coordinates = initial_coordinates if initial_coordinates is not None else adata.obsm[reduction_key]
        if net_key in adata.obsp.keys():
            G = adata.obsp[net_key]
        else:
            raise Exception("missing %s in adata.obsp of AnnData" % net_key)
    else:
        G = data

    if G is None or initial_coordinates is None:
        raise ValueError("'G' and 'initial_coordinates' cannot be NoneType.")
    if not isinstance(G, (np.ndarray, sparse.spmatrix)):
        raise ValueError("'G' must be numpy.ndarray or sparse.spmatrix.")
    if not isinstance(initial_coordinates, (np.ndarray, sparse.spmatrix)):
        raise ValueError("'initial_coordinates' must be numpy.ndarray or sparse.spmatrix.")

    G = G.astype(dtype=np.float64)
    initial_coordinates = ut.scale_matrix(initial_coordinates).astype(dtype=np.float64)

    layout = _an.layoutNetwork(
        G=G,
        initial_position=initial_coordinates,
        method=algorithm,
        presmooth_network=presmooth_network,
        sim2dist=sim2dist,
        spread=spread,
        min_dist=min_dist,
        learning_rate=learning_rate,
        n_epochs=n_epochs,
        thread_no=thread_no,
        seed=seed,
    )

    if return_raw or not isinstance(adata, AnnData):
        return layout
    else:
        adata.obsm["ACTIONred"] = initial_coordinates[:, 0:3]
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
