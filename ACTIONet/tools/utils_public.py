import random
import string
from typing import Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse


def get_data_or_split(
        adata: AnnData,
        attr: Union[str, list, pd.Series],
        groups_use: Union[str, list, None] = None,
        to_return: Optional[str] = "data",
        d: Optional[int] = 0,
        ) -> Union[list, dict]:
    to_return = to_return.lower()
    if to_return not in ["data", "levels", "split"]:
        raise ValueError("'to_return={type}' must be 'data', 'levels', or 'split'.".format(type=to_return))

    if d not in [0, 1]:
        raise ValueError("d must be dim (0 or 1) of adata")

    if isinstance(attr, str) or len(attr) == 1:
        if d == 0:
            data_vec = adata.obs[attr]
        else:
            data_vec = adata.var[attr]

    else:
        if len(attr) != adata.shape[d]:
            raise ValueError("len(attr) does not match .shape[{dim:d}] of adata".format(dim=d))
        data_vec = attr

    if data_vec is None:
        raise ValueError("Invalid split condition")
    else:
        if isinstance(data_vec, pd.Series):
            data_vec = data_vec.tolist()

    idx = list(range(0, adata.shape[d]))
    idx_dict = {}

    # Ignores 'to_return'. Always returns index dict.
    if groups_use is not None:
        idx = [i for i, e in enumerate(data_vec) if e in groups_use]
        data_vec = [e for i, e in enumerate(data_vec) if e in groups_use]

        if len(data_vec) == 0:
            raise ValueError("Invalid split condition")

        for g in set(groups_use):
            idx_dict[g] = [idx.index(i) for i, e in enumerate(data_vec) if e == g]

    else:
        for g in set(data_vec):
            idx_dict[g] = [idx.index(i) for i, e in enumerate(data_vec) if e == g]

    if to_return == "data":
        return data_vec
    elif to_return == "levels":
        fac_tup = pd.factorize(data_vec, sort=True)
        level_dict = dict({"index": fac_tup[0], "keys": fac_tup[1]})
        return level_dict
    elif to_return == "split":
        return idx_dict


def rand_suffix(N):
    str_out = "".join(
            random.choices(
                    string.ascii_uppercase + string.ascii_lowercase + string.digits, k=N
                    )
            )
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


def normalize_matrix(
        X: Union[np.ndarray, sparse.spmatrix],
        log_transform: Optional[bool] = False,
        scale_factor: Union[str, float, int, None] = None
        ) -> Union[np.ndarray, sparse.spmatrix]:
    X = X.astype(dtype=np.float64)

    if scale_factor != "median" and not isinstance(scale_factor, (int, float)) and scale_factor is not None:
        raise ValueError(f"'scale_factor' must be 'median' or numeric.")

    if sparse.issparse(X):
        row_sums = np.array(np.sum(X, axis=1))
        row_sums[row_sums == 0] = 1
        # B = np.median(row_sums) * X.multiply(1 / row_sums)
        B = X.multiply(1 / row_sums)
    else:
        X = np.array(X)
        row_sums = np.sum(X, axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1
        # B = np.median(row_sums) * (X / row_sums)
        B = X / row_sums

    if scale_factor == "median":
        B = B * np.median(row_sums)
    elif isinstance(scale_factor, (int, float)):
        B = B * scale_factor

    if log_transform:
        B = np.log1p(B)

    if sparse.issparse(B):
        B = sparse.csc_matrix(B)

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
