from typing import Optional, Union
from anndata import AnnData
import numpy as np
from scipy import sparse

from .. import misc_utils as ut
import _ACTIONet as _an


def build(
    data: Union[AnnData, np.ndarray, sparse.spmatrix],
    algorithm: Optional[str] = "k*nn",
    distance_metric: Optional[str] = "jsd",
    density: Optional[float] = 1.0,
    k: Optional[int] = 30,
    mutual_edges_only: Optional[bool] = True,
    thread_no: Optional[int] = 0,
    data_key: Optional[str] = "H_stacked",
    net_key: Optional[str] = "ACTIONet",
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False,
) -> Union[AnnData, sparse.spmatrix, None]:
    """Computes knn/k*nn graphs from input data

    Parameters
    ----------
    data : Union[AnnData, np.ndarray, sparse.spmatrix]
        `n_obs` × `n_arch` Matrix or AnnData object containing output of the 'prune_archetypes()'.
    algorithm : Optional[str], optional
        one of k*nn, by default "k*nn"
    distance_metric : Optional[str], optional
        one of jsd, ip, l2, by default "jsd"
    density : Optional[float], optional
        Controls the overall density of constructed network. Larger values results in more retained edges, by default 1.0
    k : Optional[int], optional
        Number of nearest neighbors if knn algorithm is used (ignored, otherwise), by default 30
    mutual_edges_only : Optional[bool], optional
        symmetrization strategy. Whether to use edges that are mutually nearest neighbors or not, by default True
    thread_no : Optional[int], optional
        Number of threads, by default 0
    data_key : Optional[str], optional
        If input data is an AnnData object, it instructs which `obsm` slot contains the data, by default "H_stacked"
    net_key : Optional[str], optional
        If input data is an AnnData object, it instructs which `obsp` slot contains the network, by default "ACTIONet"
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

    data_is_AnnData = isinstance(data, AnnData)

    if data_is_AnnData:
        adata = data.copy() if copy else data
        if data_key in adata.obsm.keys():
            H = adata.obsm[data_key]
        else:
            raise Exception("missing %s in adata.obsm of AnnData" % data_key)
    else:
        H = data

    H = H.T.astype(dtype=np.float64)
    if sparse.issparse(H):
        H = H.toarray()

    G = _an.buildNetwork(
        H=H,
        algorithm=algorithm,
        distance_metric=distance_metric,
        density=density,
        k=k,
        thread_no=thread_no,
        mutual_edges_only=mutual_edges_only,
    )

    if return_raw or not data_is_AnnData:
        return G
    elif data_is_AnnData:
        adata.obsp[net_key] = G
        return adata if copy else None


def layout(
    data: Union[AnnData, np.ndarray, sparse.spmatrix] = None,
    algorithm: Optional[str] = "tumap",
    initial_coordinates: Union[np.ndarray, sparse.spmatrix] = None,
    compactness_level: Optional[int] = 50,
    n_epochs: Optional[int] = 1000,
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

    alg_name = algorithm.upper()
    if alg_name not in ["UMAP", "TUMAP"]:
        raise ValueError("'layout_algorithm' must be 'tumap' or 'umap'.")

    if data is not None:
        if isinstance(data, AnnData):
            adata = data.copy() if copy else data
            initial_coordinates = (
                initial_coordinates
                if initial_coordinates is not None
                else adata.obsm[reduction_key]
            )
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
