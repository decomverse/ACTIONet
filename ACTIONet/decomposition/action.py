"""ACTION decomposition for dense matrices.
"""

from typing import Optional, Union

import numpy as np
from anndata import AnnData
from scipy.sparse import issparse, spmatrix
from sklearn._config import config_context
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted

import _ACTIONet as _an


def runACTION(data: Union[AnnData, np.ndarray, spmatrix], reduction_key: Optional[str] = "ACTION", depth: Optional[int] = 10, thread_no: Optional[int] = 0, max_it: Optional[int] = 50, min_delta: Optional[float] = 1e-300, return_raw: Optional[bool] = False, copy: Optional[bool] = False) -> Union[AnnData, dict]:
    """Run Archetypal Analysis
    Parameters
    ----------
    data : Union[AnnData, np.ndarray, sparse.spmatrix]
        Matrix or AnnData object of shape `n_obs` × `n_vars`.
    reduction_key : Optional[str]
        Key for reduction stored in obsm/varm (only if return_raw == False), by default "ACTION"
    depth : int
        Number of archetypes, by default 10
    thread_no
        Number of threads. Defaults to number of threads available - 2.
    max_it : Optional[int]
        Maximum number of AA inner iterations, by default 50
    min_delta : Optional[float]
        Convergence parameter for AA, by default 1e-100
    return_raw : Optional[bool]
        Returns raw output of 'run_ACTION()' as dict, by default False
    copy : Optional[bool]
        If an :class:`~anndata.AnnData` is passed, determines whether a copy is returned. Is ignored otherwise, by default False

    Returns
    -------
    Union[AnnData, dict]
        Result of archetypal analysis.
        If return_raw is False or the input data is not an AnnData object, the output is a dictionary:
        "W": Archetypal matrix,
        "H": Loading matrix,
        "C": Coefficient matrix,
        Otherwise, these values are stored in the input AnnData object obsm/varm with `reduction_key` prefix.
    """

    # cast Optional types to the desired type
    reduction_key = str(reduction_key)
    depth = int(str(depth))

    if isinstance(data, AnnData):
        if reduction_key not in data.obsm.keys():
            raise ValueError("Did not find data.obsm['" + reduction_key + "'].")
        else:
            X = data.obsm[reduction_key]
    else:
        X = data

    if issparse(X):
        X = X.toarray()

    X = X.T.astype(dtype=np.float64)

    ACTION_out = _an.run_ACTION(
        S_r=X,
        k_min=depth,
        k_max=depth,
        thread_no=thread_no,
        max_it=max_it,
        min_delta=min_delta,
    )

    ACTION_out = {
        "C": ACTION_out["C"][depth - 1],
        "H": ACTION_out["H"][depth - 1],
    }
    ACTION_out["W"] = np.matmul(X, ACTION_out["C"])

    if return_raw or not isinstance(data, AnnData):
        return ACTION_out
    else:
        data.obsm[reduction_key + "_" + "C"] = ACTION_out["C"]
        data.obsm[reduction_key + "_" + "H"] = ACTION_out["H"].T
        data.varm[reduction_key + "_" + "W"] = ACTION_out["W"]

        return data if copy else None


class ACTION(TransformerMixin, BaseEstimator):
    """Dimensionality reduction using ACTION decomposition.
    This transformer performs ACTION decomposition on input data.

    The objective function is:
       .. math::
            0.5 * ||X - AZ||_F^2
    Where:
        :math: Z = BX,
        :math: 0 <= B_ij, A_ij,
        :math: sum_{i} A_ij = 1,
        :math: sum_{i} B_ij = 1,

    Parameters
    ----------
    n_components : int, default=None
        Number of components.
    max_it : int, default=100
        Number of iterations for AA solver.
    min_delta : float, default=1e-6
        min_deltaerance of the stopping condition for AA solver.


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
    from decomp import ACTION
    import numpy as np
    X = np.random.normal(0, 1, (100, 30))
    action = ACTION(n_components=5)
    action.fit(X)
    ACTION(n_components=5, max_it=100)
    print(action.reconstruction_err_)
    print(action.components_)
    """

    def __init__(self, n_components=None, *, max_it=100, min_delta=1e-16, thread_no=0):
        self.n_components = n_components
        self.max_it = max_it
        self.min_delta = min_delta
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

        out = runACTION(
            data=X,
            reduction_key=None,
            depth=self.n_components,
            thread_no=self.thread_no,
            max_it=self.max_it,
            min_delta=self.min_delta,
            return_raw=True,
            copy=False,
        )

        B, A, Z = out["C"].T, out["H"].T, out["W"].T

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
        X = self._validate_data(X, accept_sparse=False, dtype=[np.float64, np.float32], reset=False)

        Z = self.components_

        with config_context(assume_finite=True):
            A = _an.run_simplex_regression(Z.T, X.T).T

        return A
