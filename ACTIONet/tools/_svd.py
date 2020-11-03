from typing import Literal, Optional, Tuple, Union

import numpy as np
from scipy.sparse import issparse, spmatrix

import _ACTIONet as _an

def _irlb(X, dim, n_iters=1000, random_state=0):
    if issparse(X):
        return _an.IRLB_SVD(X, dim, n_iters, random_state)
    return _an.IRLB_SVD_full(X, dim, n_iters, random_state)

def _feng(X, dim, n_iters=5, random_state=0):
    if issparse(X):
        return _an.FengSVD(X, dim, n_iters, random_state)
    return _an.FengSVD_full(X, dim, n_iters, random_state)

def _halko(X, dim, n_iters=5, random_state=0):
    if issparse(X):
        return _an.HalkoSVD(X, dim, n_iters, random_state)
    return _an.HalkoSVD_full(X, dim, n_iters, random_state)

def svd(
    X: Union[np.ndarray, spmatrix],
    dim: int,
    solver: Literal[0, 1, 2] = 0,
    n_iters: Optional[int] = None,
    seed: int = 0
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Randomized singular value decomposition (SVD)

    Parameters
    ----------
    X
        Input matrix (dense or sparse)
    dim
        Target dimensions
    solver
        SVD solver to use:
        `0` (the default)
          randomized SVD used in IRLBA R package
        `1`
          randomized SVD from Halko et al.
        `2`
          randomized SVD from Feng et al.
    n_iters
        Maximum number of iterations. Defaults to 1000 for `irlb` and 5 otherwise.
    seed
        Random seed

    Returns
    -------
    U
        Matrix of left singular vectors
    sigma
        Matrix of singular values
    V
        Matrix of right singular vectors
    """
    if solver == 0:
        result = _irlb(X, dim, n_iters or 1000, seed)
    elif solver == 1:
        result = _feng(X, dim, n_iters or 5, seed)
    elif solver == 2:
        result = _halko(X, dim, n_iters or 5, seed)
    else:
        raise Exception(f'Unknown SVD solver: {solver}')

    return (
        result['U'],
        result['sigma'],
        result['V'],
    )
