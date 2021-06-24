from typing import Optional, Union
from anndata import AnnData
import numpy as np
from scipy import sparse

import _ACTIONet as _an


def build_network(
        data: Union[AnnData, np.ndarray, sparse.spmatrix],
        net_name_out: Optional[str] = "ACTIONet",
        density: Optional[float] = 1.0,
        thread_no: Optional[int] = 0,
        mutual_edges_only: Optional[bool] = True,
        distance_metric: Optional[str] = "jsd",
        nn_approach: Optional[str] = "k*nn",
        k: Optional[int] = None,
        M: Optional[float] = 16,
        ef_construction: Optional[float] = 200,
        ef: Optional[float] = 50,
        copy: Optional[bool] = False,
        return_raw: Optional[bool] = False,
        
        
) -> Union[AnnData, sparse.spmatrix, None]:
    """

    Build ACTIONet

    Computes and returns the ACTIONet graph
    Parameters
    ----------
    data
       `n_obs` × `n_arch` Matrix or AnnData object containing output of the 'prune_archetypes()'.
    net_name_out
        If 'data' is AnnData, store output matrix G in '.obsp' with key 'net_name_out' (default="ACTIONet")
    density
        Controls the overall density of constructed network.
        Larger values results in more retained edges.
    thread_no
        Number of parallel threads used for identifying nearest-neighbors.
        Defaults to available threads on the machine.
    mutual_edges_only
        Whether to return only edges that there is a bi-directional/mutual relationship.
    distance_metric 
        one of jsd, ip, l2 (default=jsd) 
    nn_approach
        one of k*nn, knn (default=k*nn) 
    copy
        If 'data' is AnnData, return a copy instead of writing to `data`.
    return_raw
         If 'return_raw=True' and 'data' is AnnData, return sparse adjacency matrix directly.

    Returns
    -------
    adata : anndata.AnnData
        if 'adata' given and `copy=True` returns None or else adds fields to `adata`:

        `.obsp[net_name_out]`

    G : scipy.sparse.spmatrix
        Sparse adjacency matrix encoding ACTIONet if 'return_raw=True'

    """

    data_is_AnnData = isinstance(data, AnnData)

    if data_is_AnnData:
        adata = data.copy() if copy else data
        if "H_stacked" in adata.obsm.keys():
            H_stacked = adata.obsm["H_stacked"]
        else:
            raise Exception("missing H_stacked key in adata.obsm of AnnData")
    else:
        H_stacked = data.X

    H_stacked = H_stacked.T.astype(dtype=np.float64)
    if sparse.issparse(H_stacked):
        H_stacked = H_stacked.toarray()

    if distance_metric=="jsd":
        H_stacked=H_stacked/np.sum(H_stacked,axis=0)

    G = _an.build_ACTIONet(H_stacked=H_stacked,
                           density=density,
                           thread_no=thread_no,
                           mutual_edges_only=mutual_edges_only,
                           distance_metric=distance_metric,
                           nn_approach=nn_approach,
                           k=k,
                           M=M,
                           ef_construction=ef_construction,
                           ef=ef)

    if return_raw or not data_is_AnnData:
        return G
    elif data_is_AnnData:
        adata.obsp[net_name_out] = G
        return adata if copy else None

