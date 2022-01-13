"""ACTION decomposition for dense matrices.
"""

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn._config import config_context
from sklearn.utils.validation import check_is_fitted

from anndata import ad

from typing import Optional, Union
from anndata import AnnData

import _ACTIONet as _an


def runACTION(
    data: Union[AnnData, np.ndarray],
    reduction_key: Optional[str] = "ACTION",
    k_min: Optional[int] = 2,
    k_max: Optional[int] = 30,
    thread_no: Optional[int] = 0,
    max_it: Optional[int] = 50,
    min_delta: Optional[float] = 1e-300,
) -> dict:
    """\
    Run ACTION decomposition [Mohammadi, 2018]

    Computes reduced ACTION decomposition.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    reduction_key:
        Key of '.obms' to use as input for ACTION (default="ACTION")
        Ignored if 'data' is not an 'AnnData' object.
    k_min, k_max
        Min. and max. # of archetypes to consider.
    thread_no
        Number of threads. Defaults to number of threads available - 2.
    max_it, min_delta:
        Define stopping conditions of inner AA loop

    Returns
    -------
    ACTION_out
        A dictionary with traces of C and H matrices
    """
    data_is_AnnData = isinstance(data, AnnData)

    if data_is_AnnData:
        if reduction_key not in data.obsm.keys():
            raise ValueError("Did not find adata.obsm['" + reduction_key + "'].")
        else:
            X = data.obsm[reduction_key]
    else:
        X = data

    X = X.T.astype(dtype=np.float64)
    ACTION_out = _an.run_ACTION(X, k_min, k_max, thread_no, max_it, min_delta)

    return ACTION_out


class ACTION(TransformerMixin, BaseEstimator):
    """Dimensionality reduction using ACTION decomposition.
    This transformer performs ACTION decomposition on input data.

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
    ArchetypalAnalysis : Archetypal analysis
    ACTIONMR : Multi-resolution extension of the ACTION decomposition.

    References
    ----------
    A geometric approach to characterize the functional identity of single cells
    Mohamadi, et al., 2018 https://www.nature.com/articles/s41467-018-03933-2

    Examples
    --------
    >>> from decomp import ACTION
    >>> import numpy as np
    >>> X = np.random.normal(0, 1, (100, 30))
    >>> action = ACTION(n_components=5)
    >>> action.fit(X)
    ACTION(n_components=5, n_iter=100)
    >>> print(action.reconstruction_err_)
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

    def fit_transform(self, X, y=None):
        """Learn an ACTION decomposition for the data X and returns the transformed data.
        This is more efficient than calling fit followed by transform.

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        y : Ignored
            Not used, present for API consistency by convention.

        Returns
        -------
        A : ndarray of shape (n_samples, n_components)
            Transformed data.
        """
        X = self._validate_data(X, accept_sparse=False)

        with config_context(assume_finite=True):
            A, B, Z = self._fit_transform(X)

        self.coeff_ = B
        self.components_ = Z
        self.n_components_ = Z.shape[0]
        self.n_features_in_ = Z.shape[1]
        self.reconstruction_err_ = np.linalg.norm(X - np.matmul(A, Z), "fro")

        return A

    def _fit_transform(self, X, y=None):
        """Learn an ACTION model for the data X and returns the transformed data.
        This is more efficient than calling fit followed by transform.

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        y : Ignored
            Not used, present for API consistency by convention.

        Returns
        -------
        A: ndarray of shape (n_features, n_samples)
            Loading matrix
        B: ndarray of coefficient matrix (n_components, n_samples)
            Archetype Matrix
        Z: ndarray of shape (n_components, n_features)
            Archetype Matrix
        """

        out = _an.run_ACTION(
            S_r=X.T,
            k_min=self.n_components,
            k_max=self.n_components,
            thread_no=self.thread_no,
            max_it=self.n_iter,
            min_delta=self.tol,
        )

        B = out["C"][self.n_components - 1].T
        A = out["H"][self.n_components - 1].T
        Z = np.matmul(B, X)

        return A, B, Z

    def transform(self, X):
        """Transform the data X according to the fitted ACTION model.

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

