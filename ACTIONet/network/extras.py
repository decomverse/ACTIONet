from typing import Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

import _ACTIONet as _an


def diffusion(
    data: Union[AnnData, np.ndarray, sparse.spmatrix],
    scores: Optional[Union[np.ndarray, sparse.spmatrix]] = None,
    algorithm: Optional[str] = "pagerank_chebyshev",
    alpha_val: Optional[float] = 0.85,
    max_it: Optional[int] = 5,
    threshold: Optional[float] = 1e-8,
    thread_no: Optional[int] = 0,
    net_key: Optional[str] = "ACTIONet",
    scores_key: Optional[str] = "H_stacked",
    smoothed_scores_key: Optional[str] = "archetype_footprint",
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False,
) -> Union[AnnData, np.ndarray, None]:
    """Computes smoothed scores using network diffusion

    Parameters
    ----------
    data : Union[AnnData, np.ndarray, sparse.spmatrix]
        Adjacency matrix of the input graph or AnnData object containing the network.
    scores : Union[np.ndarray, sparse.spmatrix], optional
        Input scores, by default None
    algorithm : Optional[str], optional
        Diffusion algorithm to use. Can be "pagerank", "pagerank_sym", by default "pagerank"
    alpha_val : Optional[float], optional
        Diffusion parameter. Larger values results in more long-range diffusion, by default 0.85
    max_it : Optional[int], optional
        [description], by default 5
    threshold : Optional[float], optional
        [description], by default 1e-8
    thread_no : Optional[int], optional
        Number of threads to use, by default 0
    net_key : Optional[str], optional
        Key of 'adata.obsp' containing adjacency matrix to use. (default="ACTIONet")
        Ignored if 'adata=None'.
    scores_key : Optional[str], optional
        Key of 'adata.obsm' containing scores. (default="H_stacked")
        Ignored if `adata=None`
    smoothed_scores_key : Optional[str], optional
        Key of 'adata.obsm' to store smoothed scores. (default="archetype_footprint")
        Ignored if `adata=None`
    copy : Optional[bool], optional
        If 'adata' is given, return a copy instead of writing to `adata` (default=False)
    return_raw : Optional[bool], optional
        If 'adata' is given, return array of raw node centrality scores instead of storing to 'adata' (default=False)

    Returns
    -------
    adata : anndata.AnnData
        if 'adata' given and `copy=True` returns None or else adds fields to `adata`:

        `.obsm["archetype_footprint"]`

    smoothed_scores : np.ndarray
        If `adata=None` or `return_raw=True`, returns array of archetype footprint.
    """
    alg_name = algorithm.lower()
    if alg_name not in ["pagerank", "pagerank_sym"]:
        raise ValueError("'layout_algorithm' must be 'pagerank' or 'pagerank_sym'.")

    data_is_AnnData = isinstance(data, AnnData)
    if data_is_AnnData:
        adata = data.copy() if copy else data
        scores = scores if scores is not None else adata.obsm[scores_key]
        if net_key in adata.obsp.keys():
            G = adata.obsp[net_key]
        else:
            raise Exception("missing %s in adata.obsp of AnnData" % net_key)
    else:
        G = data

    if G is None or scores is None:
        raise ValueError("'G' and 'scores' cannot be NoneType.")
    if not isinstance(G, (np.ndarray, sparse.spmatrix)):
        raise ValueError("'G' must be numpy.ndarray or sparse.spmatrix.")
    if not isinstance(scores, (np.ndarray, sparse.spmatrix)):
        raise ValueError("'scores' must be numpy.ndarray or sparse.spmatrix.")

    G = G.astype(dtype=np.float64)
    scores = scores.astype(dtype=np.float64)

    if not sparse.issparse(G):
        G = sparse.csc_matrix(G)

    if sparse.issparse(scores):
        scores = scores.toarray()

    if algorithm == "pagerank":
        smoothed_scores = _an.compute_network_diffusion_Chebyshev(
            G=G,
            X0=scores,
            thread_no=thread_no,
            alpha=alpha_val,
            max_it=max_it,
            res_threshold=threshold,
            norm_type=0,
        )
    elif algorithm == "pagerank_sym":
        smoothed_scores = _an.compute_network_diffusion_Chebyshev(
            G=G,
            X0=scores,
            thread_no=thread_no,
            alpha=alpha_val,
            max_it=max_it,
            res_threshold=threshold,
            norm_type=2,
        )

    smoothed_scores = np.array(smoothed_scores, dtype=np.float64)

    if return_raw or not data_is_AnnData:
        return smoothed_scores
    else:
        adata.obsm[smoothed_scores_key] = smoothed_scores
        return adata if copy else None


def centrality(
    data: Union[AnnData, np.ndarray, sparse.spmatrix],
    annotations: Optional[Union[np.ndarray, list, pd.Series]] = None,
    algorithm: str = "personalized_coreness",
    annotations_key: Optional[str] = "assigned_archetype",
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
    annotations:
        list-like object containing sample annotations/scores of each observation (for localized measures).
    algorithm: str
        centrality algorithm. Can be "coreness", "personalized_coreness", "pagerank", "personalized_pagerank", default is "personalized_coreness"
        Required if 'adata=None'.
    annotations_key:
        Key of 'adata.obs' containing list-like object of sample annotations/scores (default="assigned_archetype").
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
        "pagerank",
        "personalized_pagerank",
        "coreness",
        "personalized_coreness",
    ]:
        raise ValueError(
            "'layout_algorithm' must be 'pagerank or 'coreness', 'personalized_pagerank', and 'personalized_coreness'"
        )

    data_is_AnnData = isinstance(data, AnnData)
    if data_is_AnnData:
        adata = data.copy() if copy else data
        annotations = (
            annotations if annotations is not None else adata.obsm[annotations_key]
        )
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

    if not annotations is None:
        if isinstance(annotations, pd.Series):
            annotations = np.array(annotations.tolist())
        else:
            annotations = annotations.toarray()

    if algorithm == "coreness":
        node_centrality = _an.compute_core_number(G)
    elif algorithm == "personalized_coreness":
        if annotations is None:
            raise ValueError("'annotations' cannot be None for personalized_coreness")

        node_centrality = _an.compute_archetype_core_centrality(G, annotations)
    elif algorithm == "pagerank":
        u = np.ones(G.shape[0])
        scores = u / u.shape[0]
        node_centrality = diffusion(G, algorithm="pagerank", scores=scores)
    elif algorithm == "personalized_pagerank":
        if annotations is None:
            raise ValueError("'annotations' cannot be None for personalized_coreness")

        u = annotations
        scores = u / u.shape[0]
        node_centrality = diffusion(G, algorithm="pagerank", scores=scores)

    node_centrality = np.array(node_centrality, dtype=np.float64)

    if return_raw or not data_is_AnnData:
        return node_centrality
    else:
        adata.obs["node_centrality"] = node_centrality
        return adata if copy else None
