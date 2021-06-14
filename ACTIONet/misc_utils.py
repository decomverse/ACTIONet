import string
import random
import numpy as np
import pandas as pd
from scipy.sparse import issparse
from typing import Optional, Union
from anndata import AnnData


def rand_suffix(N):
    str_out = ''.join(random.choices(string.ascii_uppercase + string.ascii_lowercase + string.digits, k=N))
    return str_out


def scale_matrix(
    X: np.ndarray,
    center: Optional[bool] = True,
    scale: Optional[bool] = True,
):
    X = X.astype(dtype=np.float64)
    if center:
        X -= X.mean(axis=0)
    if scale:
        X /= X.std(axis=0, ddof=1)
    return X


def rescale_matrix(
    X: np.ndarray,
    log_scale: Optional[bool] = False,
) -> np.ndarray:

    X = X.astype(dtype=np.float64)
    if issparse(X):
        row_sums = np.array(np.sum(X, axis=1))
        row_sums[row_sums == 0] = 1
        B = np.median(row_sums) * X.multiply(1/row_sums)
    else:
        X = np.array(X)
        row_sums = np.sum(X, axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1
        B = np.median(row_sums) * (X / row_sums)

    if log_scale:
        B = np.log1p(B)
    return B


def double_normalize(
    X: np.ndarray,
    log1p: Optional[bool] = True,
    min_threshold: Optional[int] = 0,
    copy: Optional[bool] = True,
) -> np.ndarray:
    X = X.copy() if copy else X

    X[X < min_threshold] = 0
    if np.max(X) > 100 and log1p:
        X = np.log1p(X)
    np.nan_to_num(X, copy=False, nan=0.0)

    rs = np.sqrt(np.sum(X, axis=1))
    rs[rs == 0] = 1
    D_r = np.diag(np.full(X.shape[0], 1 / rs))

    cs = np.sqrt(np.sum(X, axis=0))
    cs[cs == 0] = 1
    D_c = np.diag(np.full(X.shape[1], 1 / cs))

    X_scaled = D_r @ X @ D_c
    X_scaled = X_scaled / np.max(X_scaled)
    return X_scaled


def get_attr_or_split_idx(
    adata: AnnData,
    attr: Union[str, list],
    groups_use: Union[str, list, None] = None,
    return_vec: Optional[bool] = False,
    d: Optional[int] = 0
) -> [list, dict]:

    if d not in [0, 1]:
        raise ValueError(
            "d must be dim (0 or 1) of adata"
        )

    if isinstance(attr, str) or len(attr) == 1:
        if d == 0:
            split_vec = adata.obs[attr]
        else:
            split_vec = adata.var[attr]

    else:
        if len(attr) != adata.shape[d]:
            raise ValueError(
                "len(attr) does not match shape[{0:d}] of adata".format(d)
            )
        split_vec = attr

    if split_vec is None:
        raise ValueError(
            "Invalid split condition"
        )
    else:
        if isinstance(split_vec, pd.Series):
            split_vec = split_vec.tolist()

    idx = list(range(0, adata.shape[d]))
    idx_out = {}

    if groups_use is not None:
        idx = [i for i, e in enumerate(split_vec) if e in groups_use]
        split_vec = [e for i, e in enumerate(split_vec) if e in groups_use]

        if len(split_vec) == 0:
            raise ValueError(
                "Invalid split condition"
            )

        for g in set(groups_use):
            idx_out[g] = [idx.index(i) for i, e in enumerate(split_vec) if e == g]

    else:
        for g in set(split_vec):
            idx_out[g] = [idx.index(i) for i, e in enumerate(split_vec) if e == g]

    return split_vec if return_vec else idx_out
