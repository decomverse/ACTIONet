import random
from typing import Optional, Union

import numpy as np
import pandas as pd
import scipy
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
        raise ValueError("'algorithm' must be 'pagerank' or 'pagerank_sym'.")

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
        label_attr: Union[str, np.ndarray, list, pd.Series] = "assigned_archetype",
        algorithm: Optional[str] = "coreness",
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
    label_attr:
        list-like object or key in .obs containing labels/scores of each observation (for localized measures).
    algorithm: str
        centrality algorithm. Can be "coreness", "localized_coreness", "pagerank", "localized_pagerank", default is "localized_coreness"
        Required if 'adata=None'.
    label_attr:
        Key of 'adata.obs' containing list-like object of sample label_attr/scores (default="assigned_archetype").
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
    alg_name = algorithm.lower()
    if alg_name not in ["coreness", "pagerank", "localized_coreness", "localized_pagerank"]:
        raise ValueError("'algorithm' must be 'coreness', 'pagerank', 'localized_coreness', or 'localized_pagerank'.")

    data_is_AnnData = isinstance(data, AnnData)
    if data_is_AnnData:
        adata = data.copy() if copy else data
        label_attr = (
            label_attr if label_attr is not None else adata.obs[label_attr]
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

    if not label_attr is None:
        if isinstance(label_attr, pd.Series):
            label_attr = np.array(label_attr.tolist())
        elif sparse.issparse(label_attr):
            label_attr = label_attr.toarray()

    if algorithm == "coreness":
        node_centrality = _an.compute_core_number(G)
    elif algorithm == "localized_coreness":
        if label_attr is None:
            raise ValueError("'label_attr' cannot be None for localized_coreness")

        node_centrality = _an.compute_archetype_core_centrality(G, label_attr)
    elif algorithm == "pagerank":
        u = np.ones(G.shape[0])
        scores = u / u.shape[0]
        node_centrality = diffusion(G, algorithm="pagerank", scores=scores)
    elif algorithm == "localized_pagerank":
        if label_attr is None:
            raise ValueError("'label_attr' cannot be None for localized_coreness")

        u = label_attr
        scores = u / u.shape[0]
        node_centrality = diffusion(G, algorithm="pagerank", scores=scores)

    node_centrality = np.array(node_centrality, dtype=np.float64)

    if return_raw or not data_is_AnnData:
        return node_centrality
    else:
        adata.obs[centrality_key] = node_centrality
        return adata if copy else None


def autocorrelation(
        data: Union[AnnData, np.ndarray, sparse.spmatrix],
        scores: Optional[Union[np.ndarray, sparse.spmatrix, list, pd.Series]] = None,
        algorithm: Optional[str] = "Moran",
        normalization_method: Optional[int] = 1,
        perm_no: Optional[int] = 0,
        thread_no: Optional[int] = 0,
        net_key: Optional[str] = "ACTIONet",
        scores_key: Optional[str] = None,
        ) -> np.ndarray:
    """Computes spatial autocorrelation of scores

    Parameters
    ----------
    data : Union[AnnData, np.ndarray, sparse.spmatrix]
        Adjacency matrix of the input graph or AnnData object containing the network.
    scores : Union[np.ndarray, sparse.spmatrix, list, pd.Series], optional
        Input scores, by default None
    algorithm : Optional[str], optional
        Diffusion algorithm to use. Can be "geary", "moran", default is "moran"
    normalization_method: Optional[int]
        How to normalize scores for moran/geary methods.
    perm_no: Optional[int]
        Number of permutations to compute significance of scores., default is 0 (no permutations).    
    thread_no : Optional[int]
        Number of threads to use, by default 0
    net_key : Optional[str], optional
        Key of 'adata.obsp' containing adjacency matrix to use. (default="ACTIONet")
        Ignored if 'adata=None'.
    scores_key : Optional[str], optional
        Key of 'adata.obsm' containing scores. (default="H_stacked")
        Ignored if `adata=None`
    Returns
    -------
    autocorrelations : dictionary
    """
    alg_name = algorithm.lower()
    if alg_name not in ["moran", "geary", "categorical"]:
        raise ValueError("algorithm must be 'moran', 'geary', or 'categorical'.")

    data_is_AnnData = isinstance(data, AnnData)
    if data_is_AnnData:
        adata = data
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

    if sparse.issparse(scores):
        scores = scores.toarray()

    if algorithm == "geary":
        if not sparse.issparse(G):
            auto_out = _an.autocorrelation_Geary(
                    G,
                    scores,
                    normalization_method=normalization_method,
                    perm_no=perm_no,
                    thread_no=thread_no,
                    )
        else:
            auto_out = _an.autocorrelation_Geary_full(
                    G,
                    scores,
                    normalization_method=normalization_method,
                    perm_no=perm_no,
                    thread_no=thread_no,
                    )
    elif algorithm == "moran":
        if not sparse.issparse(G):
            auto_out = _an.autocorrelation_Moran(
                    G,
                    scores,
                    normalization_method=normalization_method,
                    perm_no=perm_no,
                    thread_no=thread_no,
                    )
        else:
            auto_out = _an.autocorrelation_Moran_full(
                    G,
                    scores,
                    normalization_method=normalization_method,
                    perm_no=perm_no,
                    thread_no=thread_no,
                    )
    elif algorithm == "categorical":
        if not scores is None:
            if isinstance(scores, pd.Series):
                scores = np.array(scores.tolist())
            elif sparse.issparse(scores):
                scores = scores.toarray()
        auto_out = assess_categorical_autocorrelation(G, scores, perm_no)

    return auto_out


def compute_phi(A, labels, s0, s1, s2):
    n = len(labels)
    # handle labels passed as list by wrapping them in a pandas series object
    if type(labels) != type(pd.Series):
        labels = pd.Series(labels)
    counts = labels.value_counts()
    categories = counts.index.tolist()
    w = A.data
    pvec = counts / n
    k = len(pvec)
    m1_rawphi = s0 / (n * (n - 1)) * (n**2 * k * (2 - k) - n * sum(1 / pvec))
    Q1 = sum(1 / pvec)
    Q2 = sum(1 / np.power(pvec, 2))
    Q3 = sum(1 / np.power(pvec, 3))
    Q22 = np.sum(np.expand_dims(1 / pvec, axis=1) * np.expand_dims(1 / pvec, axis=0))
    E1 = (np.power(n, 2) * Q22 - n * Q3) / (n * (n - 1))
    E2 = (
            4 * np.power(n, 3) * Q1
            - 4 * np.power(n, 3) * k * Q1
            + np.power(n, 3) * np.power(k, 2) * Q1
            - 2 * (2 * np.power(n, 2) * Q2 - np.power(n, 2) * k * Q2)
            + 2 * n * Q3
            - np.power(n, 2) * Q22
    )
    E2 = E2 / (n * (n - 1) * (n - 2))
    A1 = (
            4 * np.power(n, 4) * np.power(k, 2)
            - 4 * np.power(n, 4) * np.power(k, 3)
            + np.power(n, 4) * np.power(k, 4)
            - (2 * np.power(n, 3) * k * Q1 - np.power(n, 3) * np.power(k, 2) * Q1)
    )
    A2 = (
            4 * np.power(n, 3) * Q1
            - 4 * np.power(n, 3) * k * Q1
            + np.power(n, 3) * np.power(k, 2) * Q1
            - (2 * np.power(n, 2) * Q2 - np.power(n, 2) * k * Q2)
    )
    Apart = A1 - 2 * A2
    B1 = (
            4 * np.power(n, 3) * Q1
            - 4 * np.power(n, 3) * k * Q1
            + np.power(n, 3) * np.power(k, 2) * Q1
            - (2 * np.power(n, 2) * Q2 - np.power(n, 2) * k * Q2)
    )
    B2 = 2 * np.power(n, 2) * Q2 - np.power(n, 2) * k * Q2 - n * Q3
    B3 = np.power(n, 2) * Q22 - n * Q3
    Bpart = B1 - B2 - B3
    C1 = (
            2 * np.power(n, 3) * k * Q1
            - np.power(n, 3) * np.power(k, 2) * Q1
            - np.power(n, 2) * Q22
    )
    C2 = 2 * np.power(n, 2) * Q2 - np.power(n, 2) * k * Q2 - n * Q3
    Cpart = C1 - 2 * C2
    E3 = (Apart - 2 * Bpart - Cpart) / (n * (n - 1) * (n - 2) * (n - 3))
    m2_rawphi = s1 * E1 + (s2 - 2 * s1) * E2 + (np.power(s0, 2) - s2 + s1) * E3
    v_i = labels[A.row]
    v_j = labels[A.col]
    p_i = np.asarray(pvec[v_i])
    p_j = np.asarray(pvec[v_j])
    rawphi = int(
            sum(
                    w
                    * (2 * (v_i.reset_index(drop=True) == v_j.reset_index(drop=True)) - 1)
                    / (p_i * p_j)
                    )
            )
    mean_rawphi = m1_rawphi
    var_rawphi = m2_rawphi - np.power(mean_rawphi, 2)
    phi_z = (rawphi - mean_rawphi) / np.sqrt(var_rawphi)
    phi_logPval = -1 * np.log10(scipy.stats.norm.sf(phi_z))
    dictz = phi_z
    logPval = phi_logPval
    phi = rawphi
    return dictz, logPval, phi


def assess_categorical_autocorrelation(G, labels: list, perm_no: int = 100):
    # if labels is a list of strings, change them to a numerical encoding
    _, labels = np.unique(labels, return_inverse=True)
    labels = labels + 1

    # labels=string_list_to_int_list(labels)
    w = G.data
    s0 = sum(w)
    s1 = sum(4 * w**2) / 2
    s2 = int(sum(np.power(G.sum(axis=1) + G.sum(axis=0).transpose(), 2)))
    G = G.tocoo()
    dictz, logPval, phi = compute_phi(G, labels, s0, s1, s2)
    rand_phis = []
    # set the random seed
    np.random.seed(0)
    random.seed(0)
    for i in range(perm_no):
        rand_labels = random.sample(labels.tolist(), labels.shape[0])
        _, _, rand_phi = compute_phi(G, rand_labels, s0, s1, s2)
        rand_phis.append(rand_phi)
    z = (phi - np.mean(rand_phis)) / np.std(rand_phis)

    out = {"phi": phi, "dictz": dictz, "logPval": logPval, "z": z}

    return out
