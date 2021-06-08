import numpy as np

import _ACTIONet as _an


def run_simplex_regression(
    A: np.ndarray, B: np.ndarray, computeXtX: bool = False
) -> np.ndarray:

    """
    Simplex-Constrained Regression (AX-B).

    Solves the linear regression problem, subject to coefficients being positive and sum to one.

    Parameters
    ----------
    A:
        Matrix of independent variables (design matrix)
    B:
        Matrix of dependent variables (response variable)

    computeXtX:
        Parameter to simplex regression

    Returns
    -------
    X:
        Coefficient matrix
    """

    return _an.run_simplex_regresion(A, B, computeXtX)
