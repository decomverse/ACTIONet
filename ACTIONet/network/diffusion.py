from typing import Optional, Union

import numpy as np
from anndata import AnnData
from scipy import sparse

import _ACTIONet as _an


def diffusion(
    data: Union[AnnData, np.ndarray, sparse.spmatrix],
    scores: Optional[Union[np.ndarray, sparse.spmatrix]] = None,
    algorithm: Optional[str] = "pagerank",
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
        If 'adata' is given, return array of raw diffusion scores instead of storing to 'adata' (default=False)

    Returns
    -------
    adata : anndata.AnnData
        if 'adata' given and `copy=True` returns None or else adds fields to `adata`:

        `.obsm["archetype_footprint"]`

    smoothed_scores : np.ndarray
        If `adata=None` or `return_raw=True`, returns array of archetype footprint.
    """
    alg_name = str(algorithm).lower()
    if alg_name not in ["pagerank", "pagerank_sym"]:
        raise ValueError("'algorithm' must be 'pagerank' or 'pagerank_sym'.")

    if isinstance(data, AnnData):
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

    if isinstance(scores, sparse.spmatrix):
        scores = scores.toarray()

    if algorithm == "pagerank":
        smoothed_scores = _an.compute_network_diffusion_approx(
            G=G,
            X0=scores,
            thread_no=thread_no,
            alpha=alpha_val,
            max_it=max_it,
            res_threshold=threshold,
            norm_type=0,
        )
    elif algorithm == "pagerank_sym":
        smoothed_scores = _an.compute_network_diffusion_approx(
            G=G,
            X0=scores,
            thread_no=thread_no,
            alpha=alpha_val,
            max_it=max_it,
            res_threshold=threshold,
            norm_type=2,
        )

    smoothed_scores = np.array(smoothed_scores, dtype=np.float64)

    if return_raw or not isinstance(adata, AnnData):
        return smoothed_scores
    else:
        adata.obsm[smoothed_scores_key] = smoothed_scores
        return adata if copy else None
