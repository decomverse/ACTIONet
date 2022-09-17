from typing import Optional, Union

import harmonypy as hm
import numpy as np
from anndata import AnnData
from scipy.sparse import issparse, spmatrix
from typing_extensions import Literal

import _ACTIONet as _an


def reduce_adata(
    data: Union[AnnData, np.ndarray, spmatrix],
    dim: Optional[int] = 50,
    max_iter: Optional[int] = 1000,
    layer_key: Optional[str] = None,
    reduction_key: Optional[str] = "ACTION",
    svd_solver: Literal[0, 1, 2] = 0,
    seed: Optional[int] = 0,
    return_raw: Optional[bool] = False,
    copy: Optional[bool] = False,
) -> Union[AnnData, np.ndarray, spmatrix, dict]:
    """Kernel Reduction Method [Mohammadi2020].  Computes SVD-reduced form of the kernel matrix.
    :param data: Matrix or AnnData object of shape `n_obs` × `n_vars`.
    :param dim: Target dimension. Defaults to 50, or 1 - minimum dimension size of selected representation.
    :param max_iter: Maximum number of iterations
    :param reduction_key: Prefix of key to use for reduction output
    :param layer_key: Key of 'layers' to use as matrix for filtering
    :param svd_solver:SVD solver to use \
        `0` (the default) \
          randomized SVD used in IRLBA R package \
        `1` \
          randomized SVD from Halko et al. \
        `2` \
          randomized SVD from Feng et al.
    :param seed: Random seed
    :param return_raw: Returns raw output of 'reduce_kernel()' as dict
    :param copy:If an :class:`~anndata.AnnData` is passed, determines whether a copy \
    is returned. Is ignored otherwise.
    ....
    :return ACTION : :class:`~scipy.sparse.spmatrix`, :class:`~numpy.ndarray` \
        If `data` is array-like and `return_raw=False` was passed, \
        this function only returns `ACTION`
    :return adata: anndata.AnnData otherwise if `copy=True` returns None or else adds fields to `adata`  \
        `.obsm['ACTION']` \
             Scaled right singular vectors (reduced cell representations) \
        `.varm['ACTION_V']` \
             Left singular vectors (signifying gene modules) \
        `.uns['ACTION']['params']` \
        `.uns['ACTION']['sigma']` \
        `.varm['ACTION_A']` \
        `.obsm['ACTION_B']` \
    :return raw_output:If 'return_raw=True' returns dict with S_r, V, sigma, A, and B matrices.
    """

    if isinstance(data, AnnData):
        adata = data.copy() if copy else data
    else:
        adata = AnnData(data)

    # ACTIONet C++ library takes cells as columns
    if layer_key is not None:
        X = adata.layers[layer_key]
    else:
        X = adata.X

    adata.uns["metadata"]["default_reduction"] = reduction_key

    # See ACTIONet.h for definitions
    # irlb  = 0
    # halko = 1
    # feng  = 2
    if svd_solver != 0:
        max_iter = 5

    X = X.T.astype(dtype=np.float64)
    if issparse(X):
        if X.getformat() != "csc":
            X = X.tocsc()

        reduced = _an.reduce_kernel(X, dim, max_iter, seed, svd_solver, False)
    else:
        X = np.array(X)
        reduced = _an.reduce_kernel_full(X, dim, max_iter, seed, svd_solver, False)

    reduced["S_r"] = reduced["S_r"].T

    if return_raw:
        reduced["V"] = reduced["V"].T
        return reduced
    elif isinstance(adata, AnnData):
        adata.obsm[reduction_key] = reduced["S_r"]
        adata.uns[reduction_key] = {}
        adata.uns[reduction_key]["sigma"] = reduced["sigma"]
        adata.obsm[str(reduction_key) + "_B"] = reduced["B"]
        adata.varm[str(reduction_key) + "_V"] = reduced["V"]
        adata.varm[str(reduction_key) + "_A"] = reduced["A"]

        adata.uns.setdefault("obsm_annot", {}).update(
            {
                str(reduction_key): {"type": np.array([b"reduction"], dtype=object)},
                str(reduction_key) + "_B": {"type": np.array([b"internal"], dtype=object)},
            }
        )
        adata.uns.setdefault("varm_annot", {}).update(
            {
                str(reduction_key) + "_A": {"type": np.array([b"internal"], dtype=object)},
                str(reduction_key) + "_V": {"type": np.array([b"internal"], dtype=object)},
            }
        )

        return adata if copy else None
    else:
        return reduced["S_r"]


def reduce_and_batch_correct_adata_Harmony(
    data: AnnData,
    batch_key: str,
    dim: Optional[int] = 50,
    max_iter: Optional[int] = 1000,
    layer_key: Optional[str] = None,
    reduction_key: Optional[str] = "ACTION",
    batch_corrected_reduction_key: Optional[str] = "Harmony",
    svd_solver: Literal[0, 1, 2] = 0,
    seed: Optional[int] = 0,
    copy: Optional[bool] = False,
) -> Union[AnnData, np.ndarray, spmatrix, dict]:
    """Kernel Reduction Method [Mohammadi2020].  Computes SVD-reduced form of the kernel matrix.
    :param data: Matrix or AnnData object of shape `n_obs` × `n_vars`.
    :param dim: Target dimension. Defaults to 50, or 1 - minimum dimension size of selected representation.
    :param max_iter: Maximum number of iterations
    :param reduction_key: Prefix of key to use for reduction output
    :param layer_key: Key of 'layers' to use as matrix for filtering
    :param svd_solver:SVD solver to use \
        `0` (the default) \
          randomized SVD used in IRLBA R package \
        `1` \
          randomized SVD from Halko et al. \
        `2` \
          randomized SVD from Feng et al.
    :param seed: Random seed
    :param return_raw: Returns raw output of 'reduce_kernel()' as dict
    :param copy:If an :class:`~anndata.AnnData` is passed, determines whether a copy \
    is returned. Is ignored otherwise.
    ....
    :return ACTION : :class:`~scipy.sparse.spmatrix`, :class:`~numpy.ndarray` \
        If `data` is array-like and `return_raw=False` was passed, \
        this function only returns `ACTION`
    :return adata: anndata.AnnData otherwise if `copy=True` returns None or else adds fields to `adata`  \
        `.obsm['ACTION']` \
             Scaled right singular vectors (reduced cell representations) \
        `.varm['ACTION_V']` \
             Left singular vectors (signifying gene modules) \
        `.uns['ACTION']['params']` \
        `.uns['ACTION']['sigma']` \
        `.varm['ACTION_A']` \
        `.obsm['ACTION_B']` \
    :return raw_output:If 'return_raw=True' returns dict with S_r, V, sigma, A, and B matrices.
    """
    reduce_adata(
        data=data,
        dim=dim,
        max_iter=max_iter,
        layer_key=layer_key,
        reduction_key=reduction_key,
        svd_solver=svd_solver,
        seed=seed,
        return_raw=False,
        copy=False,
    )

    data_mat = data.obsm[reduction_key]
    meta_data = data.obs
    vars_use = batch_key
    eps = 0.01 / np.sqrt(data.shape[1])
    data_mat[data_mat == 0] = eps

    hm_out = hm.run_harmony(data_mat, meta_data, vars_use)
    Z_corr = hm_out.Z_corr.T

    data.obsm[batch_corrected_reduction_key] = Z_corr

    data.uns["metadata"]["default_reduction"] = batch_corrected_reduction_key

    return data if copy else None
