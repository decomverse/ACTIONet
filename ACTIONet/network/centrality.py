from typing import Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

import _ACTIONet as _an
from ACTIONet.network.diffusion import diffusion


def centrality(
    data: Union[AnnData, np.ndarray, sparse.spmatrix],
    algorithm: Optional[str] = "coreness",
    labels: Union[str, np.ndarray, list, pd.Series] = "assigned_archetype",
    net_key: Optional[str] = "ACTIONet",
    centrality_key: Optional[str] = "node_centrality",
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False,
) -> Union[AnnData, np.ndarray, None]:
    """Computes node centrality scores

    Compute node centralities using different measures

    Parameters
    ----------
    data: Union[AnnData, np.ndarray, sparse.spmatrix]
        Adjacency matrix of the input graph or AnnData object containing the network.
    algorithm: str
        centrality algorithm. Can be "coreness", "localized_coreness", "pagerank", "localized_pagerank", default is "localized_coreness"
        Required if 'adata=None'.
    labels:
        Key of 'adata.obs' containing list-like object of sample labels/scores (default="assigned_archetype").
        Ignored if data is not an AnnData object.
    net_key:
        Key of 'adata.obsp' containing adjacency matrix (default="ACTIONet").
        Ignored if data is not an AnnData object.
    centrality_key:
        Key of 'adata.obsm' to store centrality scores. (default="node_centrality")
        Ignored if `adata=None`
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
    alg_name = str(algorithm).lower()
    if alg_name not in [
        "coreness",
        "pagerank",
        "localized_coreness",
        "localized_pagerank",
    ]:
        raise ValueError("'algorithm' must be 'coreness', 'pagerank', 'localized_coreness', or 'localized_pagerank'.")

    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        labels = adata.obs[labels] if type(labels) == str else labels
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

    if labels is not None:
        if isinstance(labels, pd.Series):
            labels = np.array(labels.tolist())
        elif sparse.issparse(labels):
            labels = np.array(labels)

    if algorithm == "coreness":
        node_centrality = _an.compute_core_number(G)
    elif algorithm == "localized_coreness":
        if labels is None:
            raise ValueError("'labels' cannot be None for localized_coreness")

        node_centrality = _an.compute_archetype_core_centrality(G, labels)
    elif algorithm == "pagerank":
        u = np.ones(G.shape[0])
        scores = u / u.shape[0]
        node_centrality = diffusion(G, algorithm="pagerank", scores=scores)
    elif algorithm == "localized_pagerank":
        if labels is None:
            raise ValueError("'labels' cannot be None for localized_coreness")

        u = labels
        scores = u / u.shape[0]
        node_centrality = diffusion(G, algorithm="pagerank", scores=scores)

    node_centrality = np.array(node_centrality, dtype=np.float64)

    if return_raw or not isinstance(adata, AnnData):
        return node_centrality
    else:
        adata.obs[centrality_key] = node_centrality
        return adata if copy else None
