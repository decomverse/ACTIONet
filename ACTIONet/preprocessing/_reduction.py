from typing import Literal, Optional, Union

import numpy as np
from anndata import AnnData
from scipy import sparse

import _ACTIONet as _an
from ..tools import rescale_matrix

def reduce_adata(
    adata: AnnData,
    reduction_key: Optional[str] = 'ACTION',
    normalize: Optional[bool] = True,
    dim: Optional[int] = 50,
    svd_solver: Literal[0, 1, 2] = 0,
    n_iters: Optional[int] = 5,
    seed: Optional[int] = 0,
    copy: bool = False
) -> Optional[AnnData]:
    """\
    Kernel Reduction Method [Mohammadi2020].

    Computes SVD-reduced form of the kernel matrix.

    Parameters
    ----------
    data
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
    copy
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned. Is ignored otherwise.

    Returns
    -------
    adata : anndata.AnnData
        if `copy=True` returns None or else adds fields to `adata`:

        `.obsm[f'{reduction_key}']`
             Scaled right singular vectors (reduced cell representations)
        `.varm[f'{reduction_key}_V']`
        `.varm[f'{reduction_key}_A']`
        `.obsm[f'{reduction_key}_B']`
        `.uns['metadata'][f'{reduction_key}_sigma']`

        `.uns['obsm_annot'][f'{reduction_key}']`
        `.uns['obsm_annot'][f'{reduction_key}_B']`
        `.uns['varm_annot'][f'{reduction_key}_V']`
        `.uns['varm_annot'][f'{reduction_key}_A']`
    """
    adata = adata.copy() if copy else adata

    if normalize and 'counts' not in adata.layers.keys():
        normalized = rescale_matrix(adata.X.toarray() if sparse.issparse(adata.X) else adata.X, log_scale=True)
        old_X = adata.X
        adata.X = sparse.csc_matrix(normalized) if sparse.issparse(adata.X) else normalized
        adata.layers['counts'] = old_X

    # ACTIONet C++ library takes cells as columns
    X = adata.X.T

    # See ACTIONet.h for definitions
    # irlb  = 0
    # halko = 1
    # feng  = 2
    if svd_solver == 0:
        n_iters = 100 * n_iters
    reduced = _an.reduce_kernel(
        X.tocsc(), dim, n_iters, seed, svd_solver, False
    ) if sparse.issparse(X) else _an.reduce_kernel_full(
        X, dim, n_iters, seed, svd_solver, False
    )

    # Note S_r.T
    S_r, V, sigma, A, B = (
        reduced['S_r'].T, reduced['V'], reduced['sigma'], reduced['A'], reduced['B']
    )

    adata.obsm[reduction_key] = S_r
    adata.varm[f'{reduction_key}_V'] = V
    adata.varm[f'{reduction_key}_A'] = A
    adata.obsm[f'{reduction_key}_B'] = B
    adata.uns.setdefault('metadata', {}).update({
        f'{reduction_key}_sigma': sigma
    })

    adata.uns.setdefault('obsm_annot', {}).update({
        reduction_key: {'type': np.array([b'reduction'], dtype=object)},
        f'{reduction_key}_B': {'type': np.array([b'internal'], dtype=object)},
    })
    adata.uns.setdefault('varm_annot', {
        f'{reduction_key}_A': {'type': np.array([b'internal'], dtype=object)},
        f'{reduction_key}_V': {'type': np.array([b'internal'], dtype=object)},
    })

    return adata if copy else None
