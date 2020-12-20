from typing import Optional, Union
from typing_extensions import Literal

import numpy as np
from anndata import AnnData
from scipy.sparse import issparse, spmatrix

import _ACTIONet as _an


def reduce_kernel(
    data: Union[AnnData, np.ndarray, spmatrix],
    dim: Optional[int] = 50,
    svd_solver: Literal[0, 1, 2] = 0,
    n_iters: Optional[int] = 5,
    seed: Optional[int] = 0,
    prenormalize: Optional[bool] = False,
    return_info: bool = False,
    use_highly_variable: Optional[bool] = None,
    copy: bool = False,
) -> [AnnData, np.ndarray, spmatrix]:
    """\
    Kernel Reduction Method [Mohammadi2020].

    Computes SVD-reduced form of the kernel matrix.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    dim
        Target dimension. Defaults to 50, or 1 - minimum
        dimension size of selected representation.
    svd_solver
        SVD solver to use:
        `0` (the default)
          randomized SVD used in IRLBA R package
        `1`
          randomized SVD from Halko et al.
        `2`
          randomized SVD from Feng et al.
    n_iters
        Maximum number of iterations
    seed
        Random seed
    prenormalise
        Row-normalise the input matrix before SVD
    return_info
        Only relevant when not passing an :class:`~anndata.AnnData`:
        see “**Returns**”.
    use_highly_variable
        Whether to use highly variable genes only, stored in
        `.var['highly_variable']`.
        By default uses them if they have been determined beforehand.
    dtype
        Numpy data type string to which to convert the result.
    copy
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned. Is ignored otherwise.

    Returns
    -------
    ACTION_S_r : :class:`~scipy.sparse.spmatrix`, :class:`~numpy.ndarray`
        If `data` is array-like and `return_info=False` was passed,
        this function only returns `ACTION_S_r`…
    adata : anndata.AnnData
        …otherwise if `copy=True` returns None or else adds fields to `adata`:

        `.obsm['ACTION_S_r']`
             Scaled right singular vectors (reduced cell representations)
        `.varm['ACTION_V']`
             Left singular vectors (signifying gene modules)
        `.uns['ACTION']['params']`
        `.uns['ACTION']['sigma']`
        `.varm['ACTION_A']`
        `.obsm['ACTION_B']`
    """
    data_is_AnnData = isinstance(data, AnnData)
    if data_is_AnnData:
        adata = data.copy() if copy else data
    else:
        adata = AnnData(data)

    if use_highly_variable is True and "highly_variable" not in adata.var.keys():
        raise ValueError(
            "Did not find adata.var['highly_variable']. "
            "Either your data already only consists of highly-variable genes "
            "or consider running `pp.highly_variable_genes` first."
        )
    if use_highly_variable is None:
        use_highly_variable = True if "highly_variable" in adata.var.keys() else False
    adata_comp = (
        adata[:, adata.var["highly_variable"]] if use_highly_variable else adata
    )

    # ACTIONet C++ library takes cells as columns
    X = adata_comp.X.T

    # See ACTIONet.h for definitions
    # irlb  = 0
    # halko = 1
    # feng  = 2
    if svd_solver == 0:
        n_iters = 100 * n_iters
    reduced = (
        _an.reduce_kernel(X, dim, n_iters, seed, svd_solver, prenormalize)
        if issparse(X)
        else _an.reduce_kernel_full(X, dim, n_iters, seed, svd_solver, prenormalize)
    )

    # Note S_r.T
    S_r, V, sigma, A, B = (
        reduced["S_r"].T,
        reduced["V"],
        reduced["sigma"],
        reduced["A"],
        reduced["B"],
    )

    if data_is_AnnData:
        adata.obsm["ACTION_S_r"] = S_r
        adata.uns["ACTION"] = {}
        adata.uns["ACTION"]["params"] = {"use_highly_variable": use_highly_variable}
        adata.uns["ACTION"]["sigma"] = sigma
        adata.obsm["ACTION_B"] = B

        if use_highly_variable:
            adata.varm["ACTION_V"] = np.zeros(shape=(adata.n_vars, dim))
            adata.varm["ACTION_V"][adata.var["highly_variable"]] = V
            adata.varm["ACTION_A"] = np.zeros(shape=(adata.n_vars, dim))
            adata.varm["ACTION_A"][adata.var["highly_variable"]] = A
        else:
            adata.varm["ACTION_V"] = V
            adata.varm["ACTION_A"] = A

        return adata if copy else None
    else:
        if return_info:
            return (S_r, V.T, sigma, A, B)
        else:
            return S_r
