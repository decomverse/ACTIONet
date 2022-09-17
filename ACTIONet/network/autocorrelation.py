import random
from typing import Optional, Union

import numpy as np
import pandas as pd
import scipy
from anndata import AnnData
from scipy import sparse

import _ACTIONet as _an


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
    alg_name = str(algorithm).lower()
    if alg_name not in ["moran", "geary", "categorical"]:
        raise ValueError("algorithm must be 'moran', 'geary', or 'categorical'.")

    if isinstance(data, AnnData):
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

    if isinstance(scores, sparse.spmatrix):
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
        if scores is not None:
            if isinstance(scores, pd.Series):
                scores = np.array(scores.tolist())
            elif sparse.issparse(scores):
                scores = scores.toarray()
        auto_out = assess_categorical_autocorrelation(G, scores, int(str(perm_no)))

    return auto_out


def compute_phi(A, labels, s0, s1, s2):
    n = len(labels)
    # handle labels passed as list by wrapping them in a pandas series object
    if type(labels).isinstance(type(pd.Series)) is False:
        labels = pd.Series(labels)
    counts = labels.value_counts()
    w = A.data
    pvec = counts / n
    k = len(pvec)
    m1_rawphi = s0 / (n * (n - 1)) * (n**2 * k * (2 - k) - n * sum(1 / pvec))
    Q1 = sum(1 / pvec)
    Q2 = sum(1 / np.power(pvec, 2))
    Q3 = sum(1 / np.power(pvec, 3))
    Q22 = np.sum(np.expand_dims(1 / pvec, axis=1) * np.expand_dims(1 / pvec, axis=0))
    E1 = (np.power(n, 2) * Q22 - n * Q3) / (n * (n - 1))
    E2 = 4 * np.power(n, 3) * Q1 - 4 * np.power(n, 3) * k * Q1 + np.power(n, 3) * np.power(k, 2) * Q1 - 2 * (2 * np.power(n, 2) * Q2 - np.power(n, 2) * k * Q2) + 2 * n * Q3 - np.power(n, 2) * Q22
    E2 = E2 / (n * (n - 1) * (n - 2))
    A1 = 4 * np.power(n, 4) * np.power(k, 2) - 4 * np.power(n, 4) * np.power(k, 3) + np.power(n, 4) * np.power(k, 4) - (2 * np.power(n, 3) * k * Q1 - np.power(n, 3) * np.power(k, 2) * Q1)
    A2 = 4 * np.power(n, 3) * Q1 - 4 * np.power(n, 3) * k * Q1 + np.power(n, 3) * np.power(k, 2) * Q1 - (2 * np.power(n, 2) * Q2 - np.power(n, 2) * k * Q2)
    Apart = A1 - 2 * A2
    B1 = 4 * np.power(n, 3) * Q1 - 4 * np.power(n, 3) * k * Q1 + np.power(n, 3) * np.power(k, 2) * Q1 - (2 * np.power(n, 2) * Q2 - np.power(n, 2) * k * Q2)
    B2 = 2 * np.power(n, 2) * Q2 - np.power(n, 2) * k * Q2 - n * Q3
    B3 = np.power(n, 2) * Q22 - n * Q3
    Bpart = B1 - B2 - B3
    C1 = 2 * np.power(n, 3) * k * Q1 - np.power(n, 3) * np.power(k, 2) * Q1 - np.power(n, 2) * Q22
    C2 = 2 * np.power(n, 2) * Q2 - np.power(n, 2) * k * Q2 - n * Q3
    Cpart = C1 - 2 * C2
    E3 = (Apart - 2 * Bpart - Cpart) / (n * (n - 1) * (n - 2) * (n - 3))
    m2_rawphi = s1 * E1 + (s2 - 2 * s1) * E2 + (np.power(s0, 2) - s2 + s1) * E3
    v_i = labels[A.row]
    v_j = labels[A.col]
    p_i = np.asarray(pvec[v_i])
    p_j = np.asarray(pvec[v_j])
    rawphi = int(sum(w * (2 * (v_i.reset_index(drop=True) == v_j.reset_index(drop=True)) - 1) / (p_i * p_j)))
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
    labels = [i + 1 for i in labels]

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
        rand_labels = random.sample(labels, len(labels))
        _, _, rand_phi = compute_phi(G, rand_labels, s0, s1, s2)
        rand_phis.append(rand_phi)
    z = (phi - np.mean(rand_phis)) / np.std(rand_phis)

    out = {"phi": phi, "dictz": dictz, "logPval": logPval, "z": z}

    return out
