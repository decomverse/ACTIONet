from typing import Optional

import numpy as np


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


def scale_matrix(
    X: np.ndarray,
    center: Optional[bool] = True,
    scale: Optional[bool] = True,
):
    X = X.copy()
    if center:
        X -= X.mean(axis=0)
    if scale:
        X /= X.std(axis=0, ddof=0)
    return X

def rescale_matrix(
    X: np.ndarray,
    log_scale: Optional[bool] = False,
) -> np.ndarray:
    row_sums = np.sum(X, axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    scaled = np.median(row_sums) * (X / row_sums)
    if log_scale:
        scaled = np.log1p(scaled)
    return scaled
