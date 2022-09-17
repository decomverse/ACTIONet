from typing import Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

import _ACTIONet as _an


def cluster(
    data: Union[AnnData, np.ndarray, sparse.spmatrix],
    algorithm: str = "leiden",
    resolution: Optional[float] = 1.0,
    initial_clusters: Optional[Union[np.ndarray, list, pd.Series]] = None,
    initial_clusters_key: Optional[str] = "assigned_archetype",
    final_clusters_key: Optional[str] = "Leiden",
    seed: Optional[int] = 0,
    net_key: Optional[str] = "ACTIONet",
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False,
) -> Union[AnnData, np.ndarray, None]:
    """Computes node centrality scores

    Compute node centralities using different measures

    Parameters
    ----------
    data : Union[AnnData, np.ndarray, sparse.spmatrix]
        Adjacency matrix of the input graph or AnnData object containing the network.
    algorithm: str
        centrality algorithm. Can be "Leiden", "fix", default is "Leiden"
    resolution: float
        Resolution of the Leiden clustering. Larger values results in more clusters.
    initial_clusters:
        Initial clusters.
    initial_clusters_key:
        Key of 'adata.obs' containing the initial clusters (default="assigned_archetype").
        Ignored if data is not an AnnData object.
    net_key:
        Key of 'adata.obsp' containing adjacency matrix (default="ACTIONet").
        Ignored if data is not an AnnData object.
    copy
        If 'adata' is given, return a copy instead of writing to `adata`
    return_raw
        If 'adata' is given, return array of raw node centrality scores instead of storing to 'adata'.

    Returns
    -------
    adata : anndata.AnnData
        if 'adata' given and `copy=True` returns None or else adds fields to `adata`:

        `.obs['node_centrality']`

    node_centrality : np.ndarray
        If 'adata=None' or 'return_raw=True', returns array of node centrality scores for each observation.
    """
    alg_name = algorithm.lower()
    if alg_name not in [
        "leiden",
        "fix",
    ]:
        raise ValueError("'layout_algorithm' must be 'leiden or 'fix'")

    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        initial_clusters = initial_clusters if initial_clusters is not None else adata.obs[initial_clusters_key]
        if net_key in adata.obsp.keys():
            G = adata.obsp[net_key]
        else:
            raise Exception("missing %s in adata.obsp of AnnData" % net_key)
    else:
        G = data

    if G is None:
        raise ValueError("'G' cannot be NoneType.")
    if not isinstance(G, (np.ndarray, sparse.spmatrix)):
        raise ValueError("'G' must be numpy.ndarray or sparse.spmatrix.")

    if not sparse.issparse(G):
        G = sparse.csc_matrix(G)

    if algorithm == "fix":
        final_clusters = initial_clusters
    elif algorithm == "leiden":
        if initial_clusters is None:
            initial_clusters = np.arange(G.shape[0])
        else:
            if isinstance(initial_clusters, pd.Series):
                initial_clusters = np.array(initial_clusters.tolist())
            else:
                initial_clusters = np.array(initial_clusters)

        if 0 <= G.min():
            final_clusters = _an.unsigned_cluster(
                A=G,
                resolution_parameter=resolution,
                initial_clusters=initial_clusters,
                seed=seed,
            )
        else:
            final_clusters = _an.signed_cluster(
                A=G,
                resolution_parameter=resolution,
                initial_clusters=initial_clusters,
                seed=seed,
            )
    else:
        raise ValueError("Clustering algorithm %s not found" % algorithm)

    if return_raw or not isinstance(adata, AnnData):
        return final_clusters
    else:
        adata.obs[final_clusters_key] = final_clusters
        return adata if copy else None
