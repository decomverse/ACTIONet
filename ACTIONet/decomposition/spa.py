"""SPA decomposition for dense matrices.
"""

from typing import Tuple

import numpy as np
from sklearn._config import config_context
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted

import _ACTIONet as _an


def runSPA(
    A: np.ndarray,
    k: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Successive Projection Algorithm (SPA).
    Runs SPA algorithm to solve separable NMF problem.

    Parameters
    ----------
    A : ndarray of shape (n_samples, n_features)
        Input matrix

    k : int
        Number of columns to select

    Returns
    -------
    selected_columns: ndarray of size n_components
        Selected columns from matrix X

    norms: ndarray of size n_components
        Residual norm of the matrix (loss function)
    """
    result = _an.run_SPA(A, k)

    return (
        result["selected_columns"],
        result["norms"],
    )


class SPA(TransformerMixin, BaseEstimator):
    """Dimensionality reduction using Successive Projection Algortihm (convex NMF).
    This transformer performs SPA on input data.

    Find two non-negative matrices (W, H) whose product approximates the non-
    negative matrix X. This factorization can be used for example for
    dimensionality reduction, source separation or topic extraction.

    The objective function is:
       .. math::
            0.5 * ||X - WH||_F^2
    Where:
        :math: H = X[k, ],
        :math: 0 <= W_ij,
        :math: sum_{i} W_ij = 1,

    Parameters
    ----------
    n_components : int, default=None
        Number of components.

    Attributes
    ----------
    components_ : ndarray of shape (n_components, n_features)
        Factorization matrix, sometimes called 'dictionary'.

    n_components_ : int
        The number of components. It is same as the `n_components` parameter
        if it was given. Otherwise, it will be same as the number of
        features.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    selected_samples : ndarray of size n_components
        Selected samples in the optimal solution (0-based)

    residual_norms : ndarray of size n_components
        Residual norm of consecutive basis after adding each selected columns

    See Also
    --------
    ArchetypalAnalysis : Archetypal Analysis that identifies archetypal decomposition of a matrix.
    ACTION : Archetypal analysis initialized with convex NMF.
    ACTIONMR : Multi-resolution extension of the ACTION decomposition.

    References
    ----------
    Fast and Robust Recursive Algorithms for Separable Nonnegative Matrix Factorization
    Nicolas Gillis and Stephen A. Vavasis, 2013 https://arxiv.org/abs/1208.1237

    Examples
    --------
    from ACTIONet.decomposition import SPA
    spa = SPA(n_components=5)
    spa.fit(X)
    SPA(n_components=5)
    print(spa.selected_samples)
    print(spa.components_)
    """

    def __init__(self, n_components=None):
        self.n_components = n_components

    def fit(self, X, y=None, **params):
        """Learn a convex NMF model for the data X.

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        y : Ignored
            Not used, present for API consistency by convention.

        **params : kwargs
            Parameters (keyword arguments) and values passed to
            the fit_transform instance.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        self.fit_transform(X, **params)
        return self

    def fit_transform(self, X, y=None):
        """Learn a convex NMF model for the data X.

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        y : Ignored
            Not used, present for API consistency by convention.

        **params : kwargs
            Parameters (keyword arguments) and values passed to
            the fit_transform instance.

        Returns
        -------
        W : ndarray of shape (n_samples, n_components)
            Transformed data.
        """
        X = self._validate_data(X, accept_sparse=False, dtype=[np.float64, np.float32])

        with config_context(assume_finite=True):
            selected_samples, residual_norms = self._fit_transform(X)

        H = X[selected_samples, :]
        W = _an.run_simplex_regression(H.T, X.T).T

        self.selected_samples = selected_samples
        self.residual_norms = residual_norms
        self.n_components_ = H.shape[0]
        self.components_ = H

        return W

    def _fit_transform(self, X, y=None):
        """Learn a convex NMF model for the data X.

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features)
            Input matrix

        y : Ignored
            Not used, present for API consistency by convention.

        **params : kwargs
            Parameters (keyword arguments) and values passed to
            the fit_transform instance.

        Returns
        -------
        selected_samples: ndarray of size n_components
            Selected samples from matrix X

        residual_norms: ndarray of size n_components
            Residual norm of the matrix (loss function)
        """
        out = _an.run_SPA(X.T, self.n_components)  # Run SPA

        selected_samples = out["selected_columns"].astype(int) - 1
        residual_norms = out["norms"]

        return selected_samples, residual_norms

    def transform(self, X):
        """Transform the data X according to the fitted convex NMF model.

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features)
            input matrix, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        Returns
        -------
        H : ndarray of shape (n_samples, n_components)
            Transformed data.
        """
        check_is_fitted(self)

        X = self._validate_data(X, accept_sparse=False, dtype=[np.float64, np.float32], reset=False)

        H = self.components_
        with config_context(assume_finite=True):
            W = _an.run_simplex_regression(H.T, X.T).T

        return W
