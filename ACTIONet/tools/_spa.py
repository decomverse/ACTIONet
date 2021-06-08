from typing import Optional, Tuple, Union

import numpy as np
from scipy.sparse import issparse, spmatrix

import _ACTIONet as _an


def spa(
    A: Union[np.ndarray, spmatrix],
    k: int,
) -> Tuple[np.ndarray, np.ndarray]:

    """
    Successive Projection Algorithm (SPA).

    Runs SPA algorithm to solve separable NMF problem.

    Parameters
    ----------
    A
        Matrix matrix
    k
        Number of columns to select

    Returns
    -------
    selected_columns
        Index of k selected columns
    norms
        Residual norm of selected columns
    """

    result = _an.run_SPA_rows_sparse(A, k) if issparse(A) else _an.run_SPA(A, k)

    return (
        result["selected_columns"],
        result["norms"],
    )
