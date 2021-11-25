"""SPA decomposition for dense matrices.
"""
import numpy as np
import _ACTIONet as _an

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn._config import config_context


class SPA(TransformerMixin, BaseEstimator):
    """Dimensionality reduction using Successive Projection Algortihm (convex NMF).
    This transformer performs SPA on input data.

    The objective function is:
       .. math::
            0.5 * ||X - X[, k]H||_{loss}^2
    Where:
        :math: W = X[, k],
        :math: 0 <= H_ij,
        :math: \sum_{i} H_ij = 1,

    Parameters
    ----------
    n_components : int, default=None
        Number of components, if n_components is not set all features
        are kept.

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
    selected_samples : int
        Selected samples in the optimal solution (0-based)

    References
    ----------
    A geometric approach to characterize the functional identity of single cells
    Mohamadi, et al., 2018 https://www.nature.com/articles/s41467-018-03933-2

    Examples
    --------
    >>> from _an.decomp import SPA
    >>> spa = SPA(n_components=5)
    >>> spa.fit(X)
    SPA(n_components=5)
    >>> print(spa.selected_samples)
    >>> print(spa.components_)
    """

    def __init__(self, n_components=2):
        self.n_components = n_components

    def fit(self, X, y=None, **params):
        """Learn a ACTION model for the data X.
        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
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

    def fit_transform(self, X, y=None, W=None, H=None):
        """Learn an ACTION model for the data X and returns the transformed data.
        This is more efficient than calling fit followed by transform.
        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples
            and `n_features` is the number of features.
        y : Ignored
            Not used, present for API consistency by convention.
        Returns
        -------
        W : ndarray of shape (n_samples, n_components)
            Transformed data.
        """
        X = self._validate_data(X, accept_sparse=False)

        with config_context(assume_finite=True):
            out = self._fit_transform(X)

        idx = out["selected_columns"].astype(int) - 1
        W = X[:, idx]
        H = _an.run_simplex_regression(W, X)

        self.n_components_ = H.shape[0]
        self.components_ = H
        self.selected_samples = idx

        return W

    def _fit_transform(self, X, y=None, W=None, H=None, update_H=True):
        """Learn a NMF model for the data X and returns the transformed data.
        Parameters
        ----------
        X : {array-like} of shape (n_samples, n_features)
            Data matrix to be decomposed
        y : Ignored
        Returns
        -------
        selected_indices:
            Selected samples
        """
        # Run SPA
        out = _an.run_SPA(X, self.n_components)
        return out
