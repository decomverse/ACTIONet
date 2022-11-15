from typing import Optional, Union

import numpy as np
from anndata import AnnData
from scipy import sparse

import _ACTIONet as _an


def build(
    data: Union[AnnData, np.ndarray, sparse.spmatrix],
    algorithm: Optional[str] = "k*nn",
    distance_metric: Optional[str] = "jsd",
    density: Optional[float] = 1.0,
    mutual_edges_only: Optional[bool] = True,
    k: Optional[int] = 10,
    H_key: Optional[str] = "H_stacked",
    net_key_out: Optional[str] = "ACTIONet",
    thread_no: Optional[int] = 0,
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False,
) -> Union[AnnData, sparse.spmatrix, None]:
    """Computes knn/k*nn graphs from input data

    Parameters
    ----------
    data : Union[AnnData, np.ndarray, sparse.spmatrix]
        `n_obs` × `n_arch` Matrix or AnnData object containing output of the 'prune_archetypes()'.
    algorithm : Optional[str], optional
        Algorithm to use for network construction. Options are k-nearest neighbors ('knn') and k*-nearest neighbors ('k*nn'), by default "k*nn"
    distance_metric : Optional[str], optional
        one of jsd, ip, l2, by default "jsd"
    density : Optional[float], optional
        Controls the overall density of constructed network. Larger values results in more retained edges, by default 1.0
    mutual_edges_only : Optional[bool], optional
        symmetrization strategy. Whether to use edges that are mutually nearest neighbors or not, by default True
    k : Optional[int], optional
        Number of nearest neighbors if knn algorithm is used (ignored, otherwise), by default 30
    H_key : Optional[str], optional
        If input data is an AnnData object, it instructs which `obsm` slot contains the data, by default "H_stacked"
    net_key_out : Optional[str], optional
        If input data is an AnnData object, it instructs which `obsp` slot contains the network, by default "ACTIONet"
    thread_no : Optional[int], optional
        Number of threads, by default 0
    copy : Optional[bool], optional
        If 'data' is AnnData, return a copy instead of writing to `data`, by default False
    return_raw : Optional[bool], optional
        If `return_raw=True` and `data` is AnnData, return sparse adjacency matrix directly, by default False

    Returns
    -------
    Union[AnnData, sparse.spmatrix, None]
        adata: anndata.AnnData with ACTIONet graph in adata.obsp[net_key_out] if 'adata' given and `copy=True`.`
        G:scipy.sparse.spmatrix. Sparse ACTIONet graph if 'return_raw=True' or 'data' is not 'anndata.AnnData'.
        None: Output in data.obsp[net_key_out] if 'data' is 'anndata.AnnData' and `copy=False`
    """

    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        if H_key in adata.obsm.keys():
            H = adata.obsm[H_key]
        else:
            raise Exception(f"'{H_key}' not in adata.obsm.")
    else:
        adata = None
        H = data

    H = H.T.astype(dtype=np.float64)
    if sparse.issparse(H):
        H = H.toarray()

    G = _an.buildNetwork(
        H=H,
        algorithm=algorithm,
        distance_metric=distance_metric,
        density=density,
        thread_no=thread_no,
        mutual_edges_only=mutual_edges_only,
        k=k,
    )
    G = sparse.spmatrix.tocsc(G)

    if return_raw or not isinstance(adata, AnnData):
        return G
    elif isinstance(adata, AnnData):
        adata.obsp[net_key_out] = G
        return adata if copy else None
    else:
        raise Exception("invalid state encountered")
