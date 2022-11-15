from typing import Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

import _ACTIONet as _an


def propagate(
    data: Union[AnnData, np.ndarray, sparse.spmatrix],
    labels: Optional[Union[np.ndarray, list, pd.Series]] = None,
    algorithm: str = "lpa",
    lambda_param: float = 1.0,
    iters: int = 3,
    sig_threshold: float = 3.0,
    fixed_samples: Optional[Union[np.ndarray, list, pd.Series]] = [],
    labels_key: Optional[str] = "assigned_archetype",
    net_key: Optional[str] = "ACTIONet",
    updated_labels_key: Optional[str] = None,
    thread_no: int = 0,
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False,
) -> Union[AnnData, np.ndarray, None]:
    """Computes node centrality scores

    Compute node centralities using different measures

    Parameters
    ----------
    data : Union[AnnData, np.ndarray, sparse.spmatrix]
        Adjacency matrix of the input graph or AnnData object containing the network.
    labels:
        list-like object containing sample labels of each observation (for localized measures).
    algorithm: str
        centrality algorithm. Can be "lpa", default is "lpa"
        Required if 'adata=None'.
    lambda_param: int
        LPA lambda parameter, default = 0
    iters: int
        LPA number of iterations
    sig_threshold: float
        LPA significance threshold (z-score) for flipping labels, default = 3
    fixed_samples: Optional[Union[np.ndarray, list, pd.Series]]
        Which samples to keep untouched
    labels_key:
        Key of 'adata.obs' containing list-like object of sample labels (default="assigned_archetype").
        Ignored if data is not an AnnData object.
    net_key:
        Key of 'adata.obsp' containing adjacency matrix (default="ACTIONet").
        Ignored if data is not an AnnData object.
    updated_labels_key:
        key in `adata.obs` to store updated labels.
    copy
        If 'adata' is given, return a copy instead of writing to `adata`
    return_raw
        If 'adata' is given, return array of raw node centrality scores instead of storing to 'adata'.

    Returns
    -------
    adata : anndata.AnnData
        if 'adata' given and `copy=True` returns None or else adds fields to `adata`:

        `.obs[updated_labels_key]`

    node_centrality : np.ndarray
        If 'adata=None' or 'return_raw=True', returns array of node centrality scores for each observation.
    """
    alg_name = algorithm.lower()
    if alg_name not in [
        "lpa",
    ]:
        raise ValueError("'algorithm' must be 'lpa'")

    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
        labels = labels if labels is not None else adata.obs[labels_key]
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
    else:
        G = G.tocsc()

    if (labels is None) and (labels_key is not None):
        labels = adata.obs[labels_key]

    if labels is not None:
        labels_int, uniques = pd.factorize(labels, sort=True)
    else:
        raise ValueError("labels and labels_key cannot both be None")

    if algorithm == "lpa":
        if labels is None:
            raise ValueError("'annotations' cannot be None for personalized_coreness")

        updated_labels_int = _an.run_LPA(
            G=G,
            labels=labels_int,
            lambda_param=lambda_param,
            iters=iters,
            sig_threshold=sig_threshold,
            fixed_labels=fixed_samples,
            thread_no=thread_no,
        )
        updated_labels_int = updated_labels_int.astype(int)
        updated_labels = uniques.values[updated_labels_int]
        updated_labels[updated_labels_int == -1] = None

        updated_labels = pd.Series(updated_labels)

    if return_raw or not isinstance(adata, AnnData):
        return updated_labels
    else:
        adata.obs[updated_labels_key] = updated_labels
        return adata if copy else None
