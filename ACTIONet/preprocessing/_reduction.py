from typing import Literal, Optional, Union

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
    copy: bool = False
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

        `.obsm['ACTION']`
             Scaled right singular vectors (reduced cell representations)
        `.varm['ACTION_V']`
        `.varm['ACTION_A']`
        `.obsm['ACTION_B']`

        `.uns['obsm_annot']['ACTION']`
        `.uns['obsm_annot']['ACTION_B']`
        `.uns['varm_annot']['ACTION_V']`
        `.uns['varm_annot']['ACTION_A']`
    """
    data_is_AnnData = isinstance(data, AnnData)
    if data_is_AnnData:
        adata = data.copy() if copy else data
    else:
        adata = AnnData(data)

    # ACTIONet C++ library takes cells as columns
    X = adata.X.T

    # See ACTIONet.h for definitions
    # irlb  = 0
    # halko = 1
    # feng  = 2
    if svd_solver == 0:
        n_iters = 100 * n_iters
    reduced = _an.reduce_kernel(
        X, dim, n_iters, seed, svd_solver, prenormalize
    ) if issparse(X) else _an.reduce_kernel_full(
        X, dim, n_iters, seed, svd_solver, prenormalize
    )

    # Note S_r.T
    S_r, V, sigma, A, B = (
        reduced['S_r'].T, reduced['V'], reduced['sigma'], reduced['A'], reduced['B']
    )

    if data_is_AnnData:
        adata.obsm['ACTION'] = S_r
        adata.varm['ACTION_V'] = V
        adata.varm['ACTION_A'] = A
        adata.obsm['ACTION_B'] = B

        adata.uns.setdefault('obsm_annot', {}).update({
            'ACTION': {'type': np.array([b'reduction'], dtype=object)},
            'ACTION_B': {'type': np.array([b'internal'], dtype=object)},
        })
        adata.uns.setdefault('varm_annot', {
            'ACTION_A': {'type': np.array([b'internal'], dtype=object)},
            'ACTION_V': {'type': np.array([b'internal'], dtype=object)},
        })

        return adata if copy else None
    else:
        if return_info:
            return (S_r, V, A, B)
        else:
            return S_r
