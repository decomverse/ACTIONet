"""Multi-resolution ACTION decomposition for dense matrices.
"""

from typing import Dict, Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from natsort import natsorted
from scipy.sparse import csc_matrix, issparse, spmatrix
from sklearn._config import config_context
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_is_fitted

import _ACTIONet as _an
from ACTIONet.decomposition.utils import prune_archetypes, unify_archetypes


def runACTIONMR(
    data: Union[AnnData, np.ndarray, spmatrix],
    k_min: int = 2,
    k_max: int = 30,
    normalization: Optional[int] = 0,
    max_iter: Optional[int] = 50,
    min_delta: Optional[float] = 1e-100,
    specificity_th: Optional[int] = -3,
    min_cells_per_archetype: Optional[int] = 2,
    unification_backbone_density: Optional[float] = 0.5,
    unification_resolution: Optional[float] = 1.0,
    unification_min_cluster_size: Optional[int] = 3,
    thread_no: Optional[int] = 0,
    reduction_key: Optional[str] = "ACTION",
    return_W: Optional[bool] = False,
    return_raw: Optional[bool] = False,
    copy: Optional[bool] = False,
) -> Union[AnnData, dict]:
    """Run Archetypal Analysis

    Parameters
    ----------
    data : Union[AnnData, np.ndarray, sparse.spmatrix]
        Matrix or AnnData object of shape `n_obs` × `n_vars`.
    k_min : int, default=2
        Minimum number of components/archetypes to consider.
    k_max : int, default=None
        Maximum number of components/archetypes to consider.
    normalization : int default pre-normalization for data matrix, default=1
        It can be 0 (no normalization), or 1 or 2.
    max_iter : Optional[int], optional
        Maximum number of AA inner iterations, by default 50
    min_delta : Optional[float], optional
        Convergence parameter for AA, by default 1e-100
    specificity_th : float, default=-3
        Z-score threshold for filtering likely "mixed" archetypes
    min_cells_per_archetype: int, default=2
        Minimum number of cells that are used to define an archetype
    unification_backbone_density : float Density of the backbone graph, default=0.5
        Higher the value, denser the backbone graph.
    unification_resolution : float Leiden resolution to cluster backbone graph, default=1.0
        Larger values result in potentially more unified archetypes.
    unification_min_cluster_size : int Minimum number of archetypes in each archetype cluster, default=3
        Smaller values allow higher sensitivity at the cost of noisy archetypes.
    thread_no
        Number of threads. Defaults to number of threads available - 2.
    reduction_key : Optional[str], optional
        Key for reduction stored in obsm/varm (only if return_raw == False), by default "ACTION"
    return_W : Optional[bool], optional
        Returns 'W_stacked' and 'W_unified' as part of output.
    return_raw : Optional[bool], optional
        Returns raw output as dict, by default False
    copy : Optional[bool], optional
        If an :class:`~anndata.AnnData` is passed, determines whether a copy is returned. Is ignored otherwise, by default False

    Returns
    -------
    Union[AnnData, dict]
        Result of archetypal analysis.
        If return_raw is False or the input data is not an AnnData object, the output is a dictionary:
            "W_unified": Multi-resolution Archetypal matrix,
            "H_unified": Multi-resolution Loading matrix,
            "C_unified": Multi-resolution Coefficient matrix,
            "W_stacked": Multi-level Archetypal matrix,
            "H_stacked": Multi-level Loading matrix,
            "C_stacked": Multi-level Coefficient matrix,
        Otherwise, these values are stored in the input AnnData object obsm/varm with `reduction_key` prefix.
    """

    if isinstance(data, AnnData):
        data = data.copy() if copy else data

        if reduction_key not in data.obsm.keys():
            # raise ValueError("Did not find data.obsm['" + reduction_key + "'].")
            X = data.X
        else:
            X = data.obsm[reduction_key]
    else:
        X = data

    if issparse(X):
        X = X.toarray()

    X = X.T.astype(dtype=np.float64)

    ACTION_out = _an.run_ACTION(
        S_r=X,
        k_min=k_min,
        k_max=k_max,
        thread_no=thread_no,
        max_it=max_iter,
        min_delta=min_delta,
        normalization=normalization,
    )

    pruned = prune_archetypes(
        C_trace=ACTION_out["C"],
        H_trace=ACTION_out["H"],
        adata=None,
        specificity_th=specificity_th,
        min_cells=min_cells_per_archetype,
        copy=False,
        return_raw=True,
    )

    if isinstance(pruned, Dict):
        unified = unify_archetypes(
            adata=None,
            S_r=X,
            C_stacked=pruned["C_stacked"],
            H_stacked=pruned["H_stacked"],
            backbone_density=unification_backbone_density,
            resolution=unification_resolution,
            min_cluster_size=unification_min_cluster_size,
            normalization=normalization,
            reduction_key=None,
            thread_no=thread_no,
            copy=False,
            return_raw=True,
        )

    if isinstance(pruned, Dict):
        ACTIONMR_out = {
            "C_stacked": csc_matrix(pruned["C_stacked"]),
            "H_stacked": csc_matrix(pruned["H_stacked"]),
            "C_unified": csc_matrix(unified["C_unified"]),
            "H_unified": csc_matrix(unified["H_unified"]),
            "assigned_archetype": pd.Categorical(
                values=unified["assigned_archetype"].astype(int),
                categories=natsorted(map(int, np.unique(unified["assigned_archetype"]))),
            ),
        }

    if return_raw or not isinstance(data, AnnData):
        if return_W:
            ACTIONMR_out["W_stacked"] = np.matmul(X, ACTIONMR_out["C_stacked"].toarray())
            ACTIONMR_out["W_unified"] = np.matmul(X, ACTIONMR_out["C_unified"].toarray())
        return ACTIONMR_out
    else:
        data.obsm["C_stacked"] = ACTIONMR_out["C_stacked"]
        data.obsm["H_stacked"] = ACTIONMR_out["H_stacked"]
        data.obsm["C_unified"] = ACTIONMR_out["C_unified"]
        data.obsm["H_unified"] = ACTIONMR_out["H_unified"]
        data.obs["assigned_archetype"] = ACTIONMR_out["assigned_archetype"]
        if return_W:
            data.varm["W_unified"] = np.matmul(X, ACTIONMR_out["C_unified"].toarray())
            data.varm["W_stacked"] = np.matmul(X, ACTIONMR_out["C_stacked"].toarray())

        return data if copy else None


class ACTIONMR(TransformerMixin, BaseEstimator):
    """Dimensionality reduction using Multi-resolution ACTION decomposition.
    This transformer performs Multi-resolution ACTION decomposition on input data.

    The objective function is:
       .. math::
            0.5 * ||X - AZ||_F^2
    Where:
        :math: Z = BX,
        :math: 0 <= B_ij, A_ij,
        :math: sum_{i} A_ij = 1,
        :math: sum_{i} B_ij = 1,

    The problem is solved for multiple number of archetypes and results are aggregated in a multi-resolution decomposition.

    Parameters
    ----------
    k_min : int, default=2
        Minimum number of components/archetypes to consider.

    k_max : int, default=None
        Maximum number of components/archetypes to consider.

    max_iter : int, default=100
        Number of iterations for AA solver.

    min_delta : float, default=1e-6
        Tolerance of the stopping condition for AA solver.

    specificity_th : float, default=-3
        Z-score threshold for filtering likely "mixed" archetypes

    min_cells_per_archetype: int, default=2
        Minimum number of cells that are used to define an archetype

    unification_th : float between [0, 1], default=0
        Amount of redundancy that is tolerated in the unified archetypes. Higher value results in larger number of archetypes

    thread_no: int, default=0
        Number of threads to use

    Attributes
    ----------
    components_ : ndarray of shape (n_components, n_features)
        Archetype matrix (Z in AA, or W in NMF terminology)

    n_components_ : int
        The estimated number of multi-resolution components

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    reconstruction_err_ : float
        Loss function for the final solution

    coeff_ : ndarray of shape (n_components_, n_samples)
        Multi-resolution coefficient matrix B (convex sample weights)

    stacked_coeffs: ndarray of shape (n_ml_features, n_samples)
        Concatenation of B matrices aross different "k" values, after prunning noisy archetypes
        n_ml_components is the remaining number of archetypes after prunning.

    stacked_loadings: ndarray of shape (n_ml_features, n_samples)
        Concatenation of A matrices aross different "k" values, after prunning noisy archetypes
        n_ml_components is the remaining number of archetypes after prunning.

    assigned_archetype: array of size n_samples
        Discretized sample to archetype assignments


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
    from decomp import ACTIONMR
    import numpy as np
    X = np.random.normal(0, 1, (1000, 100))
    actionmr = ACTIONMR(k_max=30)
    actionmr.fit(X)
    ACTIONMR(k_max=30)
    print(actionmr.reconstruction_err_)
    print(actionmr.components_)
    """

    def __init__(
        self,
        k_min=2,
        k_max=None,
        *,
        max_iter=100,
        min_delta=1e-16,
        specificity_th=-3,
        min_cells_per_archetype=2,
        unification_th=0,
        thread_no=0,
    ):
        self.k_min = k_min
        self.k_max = k_max
        self.max_iter = max_iter
        self.min_delta = min_delta
        self.specificity_th = specificity_th
        self.min_cells_per_archetype = min_cells_per_archetype
        self.unification_th = unification_th
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
            (
                A,
                B,
                Z,
                stacked_coeffs,
                stacked_loadings,
                assigned_archetype,
            ) = self._fit_transform(X)

        self.coeff_ = B
        self.components_ = Z
        self.n_components_ = Z.shape[0]
        self.n_features_in_ = Z.shape[1]
        self.reconstruction_err_ = np.linalg.norm(X - np.matmul(A, Z), "fro")

        self.stacked_coeffs = stacked_coeffs
        self.stacked_loadings = stacked_loadings
        self.assigned_archetype = assigned_archetype

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
        A: ndarray of shape (n_features, n_mr_components)
            Multi-resolution loading matrix
            n_mr_components is the estimated # of multi-resolution archetypes (based on unification_th)
            Computed from stacked_loadings after the unification process to retain nonredundant archetypes.

        B: ndarray of shape (n_mr_components, n_samples)
            Multi-resolution coefficient matrix.
            n_mr_components is the estimated # of multi-resolution archetypes (based on unification_th)
            Computed from stacked_coeffs after the unification process to retain nonredundant archetypes.

        Z: ndarray of shape (n_mr_features, n_features)
            Multi-resolution archetype matrix

        stacked_coeffs: ndarray of shape (n_ml_features, n_samples)
            Concatenation of B matrices aross different "k" values, after prunning noisy archetypes
            n_ml_components is the remaining number of archetypes after prunning.

        stacked_loadings: ndarray of shape (n_features, n_ml_features)
            Concatenation of A matrices aross different "k" values, after prunning noisy archetypes
            n_ml_components is the remaining number of archetypes after prunning.

        assigned_archetype: array of size n_samples
            Discretized sample to archetype assignments
        """

        actionmr_out = runACTIONMR(
            data=X,
            k_min=self.k_min,
            k_max=self.k_max,
            max_iter=self.max_iter,
            min_delta=self.min_delta,
            specificity_th=self.specificity_th,
            min_cells_per_archetype=self.min_cells_per_archetype,
            unification_th=self.unification_th,
            thread_no=self.thread_no,
            reduction_key=None,
            return_W=True,
            return_raw=True,
            copy=False,
        )

        stacked_coeffs = actionmr_out["C_stacked"].T
        stacked_loadings = actionmr_out["H_stacked"].T

        B = actionmr_out["C_unified"].toarray().T.astype(dtype=np.float64)
        A = actionmr_out["H_unified"].toarray().T.astype(dtype=np.float64)
        assigned_archetype = actionmr_out["assigned_archetype"]

        Z = actionmr_out["W_unified"].T.astype(dtype=np.float64)

        return A, B, Z, stacked_coeffs, stacked_loadings, assigned_archetype

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
