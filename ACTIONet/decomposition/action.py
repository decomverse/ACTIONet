"""ACTION decomposition for dense matrices.
"""

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn._config import config_context

import _ACTIONet as _an
from ACTIONet.decomposition import SPA


class ArchetypalAnalysis(TransformerMixin, BaseEstimator):
    """Dimensionality reduction using ACTION decomposition.
    This transformer performs ACTION decomposition on input data.

    The objective function is:
       .. math::
            0.5 * ||X - WH||_F^2
    Where:
        :math: W = XC,
        :math: 0 <= C_ij, H_ij,
        :math: \sum_{i} C_ij = 1,
        :math: \sum_{i} H_ij = 1,

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
        Factorization matrix, sometimes called 'dictionary'.

    n_components_ : int
        The number of components. It is same as the `n_components` parameter
        if it was given. Otherwise, it will be same as the number of
        features.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    coeff : ndarray of shape (n_components, n_features)
        Coefficient matrix C in the AA decomposition.

    References
    ----------
    A geometric approach to characterize the functional identity of single cells
    Mohamadi, et al., 2018 https://www.nature.com/articles/s41467-018-03933-2

    Examples
    --------
    >>> from decomp import ACTION
    >>> action = ACTION(n_components=5)
    >>> action.fit(X)
    ACTION(n_components=5, n_iter=100)
    >>> print(action.err_)
    >>> print(action.components_)
    """

    def __init__(self, n_components=None, *, n_iter=100, tol=1e-16, thread_no=0):
        self.n_components = n_components
        self.n_iter = n_iter
        self.tol = tol
        self.thread_no = thread_no

    def fit(self, X, y=None, **params):
        """Learn a ACTION model for the data X.

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

    def fit_transform(self, X, y=None, W=None, H=None):
        """Learn an ACTION model for the data X and returns the transformed data.
        This is more efficient than calling fit followed by transform.

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        y : Ignored
            Not used, present for API consistency by convention.

        W : array-like of shape (n_samples, n_components)
            It is used as initial guess for the solution of W.

        H : Ignored
            In future versions, it can be used to initialize AA iterations.

        Returns
        -------
        W : ndarray of shape (n_samples, n_components)
            Transformed data.
        """
        X = self._validate_data(X, accept_sparse=False)

        if W == None:
            W0 = SPA(n_components=self.n_components).fit_transform(X)
        else:
            W0 = W

        with config_context(assume_finite=True):
            C, H = self._fit_transform(X, W0=W0)

        W = np.dot(X, C)

        self.coeff = C
        self.n_components_ = H.shape[0]
        self.components_ = H

        return W

    def _fit_transform(self, X, y=None, W=None, H=None):
        """Learn an ACTION model for the data X and returns the transformed data.
        This is more efficient than calling fit followed by transform.

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        y : Ignored
            Not used, present for API consistency by convention.

        W : array-like of shape (n_samples, n_components)
            It is used as initial guess for the solution of W.

        H : Ignored
            In future versions, it can be used to initialize AA iterations.

        Returns
        -------
        C: ndarray of shape (n_samples, n_archetypes)
            Coefficient matrix (C) in archetypal analysis
        H: ndarray of shape (n_archetypes, n_features)
            Loading matrix
        """
        out = _an.run_ACTION(
            S_r=X,
            k_min=self.n_components,
            k_max=self.n_components,
            thread_no=self.thread_no,
            max_it=self.iter_no,
            min_delta=self.tol,
        )

        C = out["C"][self.n_components - 1]
        H = out["H"][self.n_components - 1]

        return C, H
