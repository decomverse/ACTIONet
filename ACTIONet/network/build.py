from typing import Optional, Union

import numpy as np
from anndata import AnnData
from scipy import sparse
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csc_matrix

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

    G = csc_matrix(G)

    if return_raw or not isinstance(adata, AnnData):
        return G
    elif isinstance(adata, AnnData):
        adata.obsp[net_key_out] = G
        return adata if copy else None
    else:
        raise Exception("invalid state encountered")


# Adopted from https://github.com/ldv1/kstar-NN and https://github.com/kfirkfir/k-Star-Nearest-Neighbors
def kStarNN_from_samples(
    X: np.ndarray,
    metric: str = "euclidean",
    L_C: float = 1.0,
    symmetrization: int = 1,
    post_normalize: bool = False,
) -> np.ndarray:
    """
    Computes k*NN graph from input samples.

    Parameters
    ----------
    X : np.ndarray
        Input samples.
    metric : str, optional
        Distance metric to use, by default "euclidean".
    L_C : float, optional
        Scaling factor for the distances, by default 1.0.
    symmetrization : int, optional
        Symmetrization strategy. 1 for OR (either node points to each other), 2 for AND (both nodes point to each other), by default 1.
    post_normalize : bool, optional
        Whether to perform post-normalization, by default False.

    Returns
    -------
    np.ndarray
        Adjacency matrix of the k*NN graph.
    """
    D = squareform(pdist(X, metric=metric))

    adj_mat = kStarNN_from_dists(
        D=D,
        L_C=L_C,
        symmetrization=symmetrization,
        post_normalize=post_normalize,
    )

    return adj_mat


def kStarNN_from_dists(
    D: np.ndarray,
    L_C: float = 1.0,
    symmetrization: int = 1,
    post_normalize: bool = False,
) -> np.ndarray:
    """
    Computes k*NN graph from input distance matrix.

    Parameters
    ----------
    D : np.ndarray
        Input distance matrix.
    L_C : float, optional
        Scaling factor for the distances, by default 1.0.
    symmetrization : int, optional
        Symmetrization strategy. 1 for OR (either node points to each other), 2 for AND (both nodes point to each other), by default 1.
    post_normalize : bool, optional
        Whether to perform post-normalization, by default False.

    Returns
    -------
    np.ndarray
        Adjacency matrix of the k*NN graph.
    """
    num_samples = D.shape[0]

    adj_mat = np.zeros_like(D)
    for i in range(num_samples):
        dists = D[:, i]
        sortIndex = np.argsort(dists)
        beta = np.append(L_C * dists[sortIndex], 10**6)
        lambda_ = beta[0] + 1

        k = 0
        Sum_beta = 0
        Sum_beta_square = 0
        while (lambda_ > beta[k]) and (k < len(beta) - 1):
            k += 1
            Sum_beta += beta[k]
            Sum_beta_square += beta[k] ** 2
            lambda_ = (1 / k) * (Sum_beta + np.sqrt(k + Sum_beta**2 - k * Sum_beta_square))

        w = np.maximum(lambda_ - L_C * dists, 0)
        total_sum = np.sum(w)
        if total_sum > 0:
            adj_mat[:, i] = w / total_sum

    if symmetrization == 1:  # if either node points to each other (OR)
        adj_mat = (adj_mat + adj_mat.T) / 2
    elif symmetrization == 2:  # if both node points to each other (AND)
        adj_mat = np.sqrt(adj_mat * adj_mat.T)  # symmetrize

    if post_normalize:
        rs = 1 / np.sqrt(adj_mat.sum(axis=1)).squeeze()
        cs = 1 / np.sqrt(adj_mat.sum(axis=0)).squeeze()
        adj_mat = rs.T * adj_mat * cs

    return adj_mat
