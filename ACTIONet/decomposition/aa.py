"""Archetypal Analysis (robust) for dense matrices.
"""

import numpy as np
from typing import Optional, Tuple

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn._config import config_context
from sklearn.utils.validation import check_is_fitted

import _ACTIONet as _an


def runAA(
    A: np.ndarray,
    k: int,
    max_iter: Optional[int] = 50,
    min_delta: Optional[float] = 1e-16
) -> Tuple[np.ndarray, np.ndarray]:
    """Archetypal Analysis (AA)
    Runs archetypal analysis.

    Parameters
    ----------
    A : ndarray of shape (n_samples, n_features)
        Input matrix

    k : int
        Number of columns to select

    max_iter, min_delta: double
        Define stopping conditions of AA

    Returns
    -------
    C
        Convex matrix of archetype coefficients (# observations x # archetypes)
    H
        Convex matrix of observation coefficients (# archetypes x # observations)
    """

    AA_out = _an.run_AA(A, k, max_iter, min_delta)
    return (
        AA_out["C"],
        AA_out["H"]
    )


class ArchetypalAnalysis(TransformerMixin, BaseEstimator):
    """Dimensionality reduction using archetypal analysis.
    This transformer performs archetypal analysis (AA) on input data.

    The objective function is:
       .. math::
            0.5 * ||X - AZ||_F^2
    Where:
        :math: Z = BX,
        :math: 0 <= B_ij, A_ij,
        :math: \sum_{i} A_ij = 1,
        :math: \sum_{i} B_ij = 1,

    Parameters
    ----------
    n_components : int, default=None
        Number of components.
    n_iter : int, default=100
        Number of iterations for AA solver.
    tol : float, default=1e-6
        Tolerance of the stopping condition for AA solver.


    Attributes
    ----------
    components_ : ndarray of shape (n_components, n_features)
        Archetype matrix (Z in AA, or W in NMF terminology)

    n_components_ : int
        The number of components. It is same as the `n_components` parameter
        if it was given.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    reconstruction_err_ : float
        Loss function for the final solution

    coeff_ : ndarray of shape (n_components, n_samples)
        Coefficient matrix B (convex sample weights) in the archetypal decomposition

    See Also
    --------
    SPA : Successive Projection Algorithm to solve convex NMF
    ACTION : ACTION decomposition (archetypal analysis initialized by convex NMF)
    ACTIONMR : Multi-resolution extension of the ACTION decomposition.

    References
    ----------
    Fast and Robust Archetypal Analysis for Representation Learning
    Chen, et al., 2014 https://arxiv.org/abs/1405.6472

    Examples
    --------
    >>> from decomp import ArchetypalAnalysis
    >>> import numpy as np
    >>> X = np.random.normal(0, 1, (100, 30))
    >>> aa = ArchetypalAnalysis(n_components=5)
    >>> aa.fit(X)
    ArchetypalAnalysis(n_components=5, n_iter=100)
    >>> print(aa.reconstruction_err_)
    >>> print(aa.components_)
    """

    def __init__(self, n_components=None, *, n_iter=100, tol=1e-16, thread_no=0):
        self.n_components = n_components
        self.n_iter = n_iter
        self.tol = tol
        self.thread_no = thread_no

    def fit(self, X, y=None, **params):
        """Learn a AA model for the data X.

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

    def fit_transform(self, X, y=None, Z=None):
        """Learn an AA model for the data X and returns the transformed data.
        This is more efficient than calling fit followed by transform.

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        y : Ignored
            Not used, present for API consistency by convention.

        Z : array-like of shape (n_components, n_features)
            It is used as initial guess for the solution of archetypes (Z).
            It uses SPA to initialize Z if it is not provided.

        Returns
        -------
        A : ndarray of shape (n_samples, n_components)
            Transformed data.
        """
        X = self._validate_data(X, accept_sparse=False)

        if Z == None:
            Xt = X.T
            selected_samples = (
                _an.run_SPA(Xt, self.n_components)["selected_columns"].astype(int) - 1
            )
            Z0 = X[selected_samples, :]
        else:
            Z0 = Z

        with config_context(assume_finite=True):
            A, B, Z = self._fit_transform(X, Z=Z0)

        self.coeff_ = B
        self.components_ = Z
        self.n_components_ = Z.shape[0]
        self.n_features_in_ = Z.shape[1]
        self.reconstruction_err_ = np.linalg.norm(X - np.matmul(A, Z), "fro")

        return A

    def _fit_transform(self, X, y=None, Z=None):
        """Learn an AA model for the data X and returns the transformed data.
        This is more efficient than calling fit followed by transform.

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        y : Ignored
            Not used, present for API consistency by convention.

        Z : ndarray of shape (n_components, n_features)
            It is used as initial guess for the solution of archetypes (Z).

        Returns
        -------
        A: ndarray of shape (n_features, n_samples)
            Loading matrix
        B: ndarray of coefficient matrix (n_components, n_samples)
            Archetype Matrix
        Z: ndarray of shape (n_components, n_features)
            Archetype Matrix
        """
        out = _an.run_AA(X.T, Z.T, self.n_iter, self.tol)

        B, Z, A = out["C"].T, out["W"].T, out["H"].T

        return A, B, Z

    def transform(self, X):
        """Transform the data X according to the fitted AA model.

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        Returns
        -------
        A : ndarray of shape (n_samples, n_components)
            Transformed data.
        """

        check_is_fitted(self)
        X = self._validate_data(
            X, accept_sparse=False, dtype=[np.float64, np.float32], reset=False
        )

        Z = self.components_

        with config_context(assume_finite=True):
            A = _an.run_simplex_regression(Z.T, X.T).T

        return A
