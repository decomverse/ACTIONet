"""Multi-resolution ACTION decomposition for dense matrices.
"""

from typing import Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from natsort import natsorted
from scipy import sparse

from sklearn.base import BaseEstimator, TransformerMixin
from sklearn._config import config_context
from sklearn.utils.validation import check_is_fitted

import _ACTIONet as _an


def prune_archetypes(
    C_trace: list,
    H_trace: list,
    adata: Optional[AnnData] = None,
    min_specificity_z_threshold: Optional[float] = -3,
    min_cells: Optional[int] = 2,
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False,
) -> Union[dict, AnnData, None]:
    """\
    Archetype pruning

    Initial pruning of archetypes

    Parameters
    ----------
    C_trace, H_trace:
        Output of pp.ACTION()
    adata
        Optional AnnData object in which to store output of 'prune_archetypes()'
    min_specificity_z_threshold:
        Controls level of prunning for non-specific archetypes
        (larger values remove more archetypes)
    min_cells:
        Minimum number of influential cells for each archetype
        to be considered nontrivial
    copy
        If 'adata' is given, return a copy instead of writing to `adata`
    return_raw
        If 'return_raw=True' and 'data' is AnnData, return raw output of 'prune_archetypes()'.
    Returns
    -------

    pruned : dict
        dict containing 'C_stacked' and 'H_stacked' matrices

    adata : anndata.AnnData
        if 'adata' given and `copy=True` returns None or else adds fields to `adata`:

        `.obsm['C_stacked']`
        `.obsm['H_stacked']`
        `.uns['obsm_annot']['C_stacked']['type']`
        `.uns['obsm_annot']['H_stacked']['type']`
    """

    if adata is not None:
        if isinstance(adata, AnnData):
            adata = adata.copy() if copy else adata
        else:
            raise ValueError("'adata' is not an AnnData object.")

    pruned = _an.prune_archetypes(
        C_trace, H_trace, min_specificity_z_threshold, min_cells
    )

    pruned["C_stacked"] = sparse.csc_matrix(pruned["C_stacked"])
    pruned["H_stacked"] = sparse.csc_matrix(pruned["H_stacked"].T)

    if return_raw or adata is None:
        return pruned
    else:
        adata.obsm["C_stacked"] = pruned["C_stacked"]
        adata.obsm["H_stacked"] = pruned["H_stacked"]
        adata.uns.setdefault("obsm_annot", {}).update(
            {
                "C_stacked": {"type": np.array([b"internal"], dtype=object)},
                "H_stacked": {"type": np.array([b"internal"], dtype=object)},
            }
        )
        return adata if copy else None


def unify_archetypes(
    adata: Optional[AnnData] = None,
    S_r: Union[np.ndarray, sparse.spmatrix, None] = None,
    C_stacked: Union[np.ndarray, sparse.spmatrix, None] = None,
    H_stacked: Union[np.ndarray, sparse.spmatrix, None] = None,
    reduction_key: Optional[str] = "ACTION",
    violation_threshold: Optional[float] = 0,
    thread_no: Optional[int] = 0,
    copy: Optional[bool] = False,
    return_raw: Optional[bool] = False,
) -> AnnData:
    """\
    Archetype unification

    Aggregates redundant archetypes.

    Parameters
    ----------
    adata:
        AnnData object containing 'reduction_key' in '.obsm' and 'net_key' in '.obsp'.
    S_r:
        Reduced representation matrix to use for unification.
        Required if 'adata=None'.
    C_stacked:
        Matrix containing output 'C_stacked' of 'prune_archetypes()' to use for unification.
        Required if 'adata=None', otherwise retrieved from 'adata.obsm["C_stacked"]'
    H_stacked:
        Matrix containing output 'H_stacked' of 'prune_archetypes()' to use for unification.
        Required if 'adata=None', otherwise retrieved from 'adata.obsm["H_stacked"]'
    reduction_key:
        Key of 'adata.obms' containing reduced matrix to use for 'S_r' in 'unify_archetypes()' (default="ACTION").
        Ignored if 'adata=None'.
    violation_threshold:
        Archetype unification resolution parameter (default=0)
    thread_no:
        Number of threads. Defaults to number of threads available - 2.
    copy
        If 'adata' is given, return a copy instead of writing to `adata`
    return_raw
        If 'adata' is given, return dict of raw 'unify_archetypes()' output instead of storing to 'adata'.
    Returns
    -------
    adata : anndata.AnnData
        if `copy=True` returns None or else adds fields to `adata`:

        `.obsm['C_unified']`
        `.obsm['H_unified']`
        `.uns['ACTION']['archetypes']['unified']`

    unified : dict
        If 'adata=None' or 'return_raw=True', a dict containing 'C_unified' and 'H_unified' matrices
    """

    if adata is not None:
        if isinstance(adata, AnnData):
            adata = adata.copy() if copy else adata
            S_r = S_r if S_r is not None else adata.obsm[reduction_key]
            C_stacked = C_stacked if C_stacked is not None else adata.obsm["C_stacked"]
            H_stacked = H_stacked if H_stacked is not None else adata.obsm["H_stacked"]
        else:
            raise ValueError("'adata' is not an AnnData object.")
    else:
        if S_r is None or C_stacked is None or H_stacked is None:
            raise ValueError("'G' and 'S_r' cannot be NoneType if 'adata=None'.")
        if not isinstance(S_r, (np.ndarray, sparse.spmatrix)):
            raise ValueError("'S_r' must be numpy.ndarray or sparse.spmatrix.")
        if not isinstance(C_stacked, (np.ndarray, sparse.spmatrix)):
            raise ValueError("'C_stacked' must be numpy.ndarray or sparse.spmatrix.")
        if not isinstance(H_stacked, (np.ndarray, sparse.spmatrix)):
            raise ValueError("'H_stacked' must be numpy.ndarray or sparse.spmatrix.")

    S_r = S_r.T.astype(dtype=np.float64)

    if sparse.issparse(C_stacked):
        C_stacked = C_stacked.toarray()

    if sparse.issparse(H_stacked):
        H_stacked = H_stacked.toarray()

    C_stacked = C_stacked.astype(dtype=np.float64)
    H_stacked = H_stacked.T.astype(dtype=np.float64)

    unified = _an.unify_archetypes(
        S_r, C_stacked, H_stacked, violation_threshold, thread_no
    )

    unified["C_unified"] = sparse.csc_matrix(unified["C_unified"])
    unified["H_unified"] = sparse.csc_matrix(unified["H_unified"].T)

    if return_raw or adata is None:
        return unified
    else:
        adata.obsm["C_unified"] = unified["C_unified"]
        adata.obsm["H_unified"] = unified["H_unified"]

        adata.uns.setdefault("obsm_annot", {}).update(
            {
                "C_unified": {"type": np.array([b"internal"], dtype=object)},
                "H_unified": {"type": np.array([b"internal"], dtype=object)},
            }
        )

        groups = unified["assigned_archetype"]
        adata.obs["assigned_archetype"] = pd.Categorical(
            values=groups.astype(int),
            categories=natsorted(map(int, np.unique(groups))),
        )

        return adata if copy else None


class ACTIONMR(TransformerMixin, BaseEstimator):
    """Dimensionality reduction using Multi-resolution ACTION decomposition.
    This transformer performs Multi-resolution ACTION decomposition on input data.

    The objective function is:
       .. math::
            0.5 * ||X - AZ||_F^2
    Where:
        :math: Z = BX,
        :math: 0 <= B_ij, A_ij,
        :math: \sum_{i} A_ij = 1,
        :math: \sum_{i} B_ij = 1,

    The problem is solved for multiple number of archetypes and results are aggregated in a multi-resolution decomposition.

    Parameters
    ----------
    depth : int, default=None
        Maximum number of components/archetypes to consider.

    n_iter : int, default=100
        Number of iterations for AA solver.

    tol : float, default=1e-6
        Tolerance of the stopping condition for AA solver.

    min_specificity_z_threshold : float, default=-3
        Z-score threshold for filtering likely "mixed" archetypes

    min_cells_per_arch: int, default=2
        Minimum number of cells that are used to define an archetype

    unification_violation_threshold : float between [0, 1], default=0
        Amount of redundancy that is tolerated in the unified archetypes. Higher value results in larger number of archetypes
        
    thread_no: int, default=0
        Number of threads to use  

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
    ACTION : ACTION decomposition.

    References
    ----------
    A multiresolution framework to characterize single-cell state landscapes
    Mohamadi, et al., 2020 https://www.nature.com/articles/s41467-020-18416-6

    Examples
    --------
    >>> from decomp import ACTIONMR
    >>> import numpy as np
    >>> X = np.random.normal(0, 1, (1000, 100))
    >>> actionmr = ACTIONMR(depth=30)
    >>> actionmr.fit(X)
    ACTIONMR(depth=30)
    >>> print(actionmr.reconstruction_err_)
    >>> print(actionmr.components_)
    """

    def __init__(
        self,
        depth=None,
        *,
        n_iter=100,
        tol=1e-16,
        min_specificity_z_threshold=-3,
        min_cells_per_arch=2,
        unification_violation_threshold=0,
        thread_no=0,
    ):
        self.depth = depth
        self.n_iter = n_iter
        self.tol = tol
        self.min_specificity_z_threshold = min_specificity_z_threshold
        self.min_cells_per_arch = min_cells_per_arch
        self.unification_violation_threshold = unification_violation_threshold
        self.thread_no = thread_no

    def fit(self, X, y=None, **params):
        """Learn a Multi-resolution ACTION model for the data X.

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
        """Learn an Multi-resolution ACTION decomposition for the data X and returns the transformed data.
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
        """Learn an Multi-resolution ACTION model for the data X and returns the transformed data.
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

        ACTION_out = _an.run_ACTION(
            S_r=X.T,
            k_min=2,
            k_max=self.depth,
            thread_no=self.thread_no,
            max_it=self.n_iter,
            min_delta=self.tol,
        )

        pruned = _an.prune_archetypes(
            ACTION_out["C"],
            ACTION_out["H"],
            self.min_specificity_z_threshold,
            self.min_cells_per_arch,
        )

        unified = _an.unify_archetypes(
            X.T,
            pruned["C_stacked"],
            pruned["H_stacked"],
            self.unification_violation_threshold,
            self.thread_no,
        )

        B = unified["C_unified"].toarray().T.astype(dtype=np.float64)
        A = unified["H_unified"].toarray().T.astype(dtype=np.float64)

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


# X = np.random.normal(0, 1, (1000, 100))
# actionmr = ACTIONMR(depth=20)
# actionmr.fit(X)
# print(actionmr.reconstruction_err_)
# print(actionmr.components_)
