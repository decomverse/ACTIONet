from typing import Optional, Union
from typing_extensions import Literal

import numpy as np
from anndata import AnnData
from scipy.sparse import issparse, spmatrix

import _ACTIONet as _an

def reduce_kernel(
    data: Union[AnnData, np.ndarray, spmatrix],
    dim: Optional[int] = 50,
    max_iter: Optional[int] = 10,
    layer_name: Optional[str] = None,
    reduction_name: Optional[str] = "ACTION",
    svd_solver: Literal[0, 1, 2] = 0,
    seed: Optional[int] = 0,
    return_info: Optional[bool] = False,
    use_highly_variable: Optional[bool] = False,
    copy: Optional[bool] = False
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
    max_iter
        Maximum number of iterations
    reduction_name
        Prefix of key to use for reduction output
    layer_name
        Key of 'layers' to use as matrix for filtering
    svd_solver
        SVD solver to use:
        `0` (the default)
          randomized SVD used in IRLBA R package
        `1`
          randomized SVD from Halko et al.
        `2`
          randomized SVD from Feng et al.
    seed
        Random seed
    return_info
        Only relevant when not passing an :class:`~anndata.AnnData`:
        see “**Returns**”.
    use_highly_variable
        Whether to use highly variable genes only, stored in
        `.var['highly_variable']`.
        By default uses them if they have been determined beforehand.
    copy
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned. Is ignored otherwise.

    Returns
    -------
    ACTION : :class:`~scipy.sparse.spmatrix`, :class:`~numpy.ndarray`
        If `data` is array-like and `return_info=False` was passed,
        this function only returns `ACTION`…
    adata : anndata.AnnData
        …otherwise if `copy=True` returns None or else adds fields to `adata`:

        `.obsm['ACTION']`
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

    if use_highly_variable is True:
        if "highly_variable" in adata.var.keys():
            adata_temp = adata[:, adata.var["highly_variable"]]
        else:
            raise ValueError(
                "Did not find adata.var['highly_variable']. "
                "Either your data already only consists of highly-variable genes "
                "or consider running `pp.highly_variable_genes` first."
            )
    else:
        adata_temp = adata

    # ACTIONet C++ library takes cells as columns
    if layer_name is not None:
        X = adata_temp.layers[layer_name].T
    else:
        X = adata_temp.X.T

    # See ACTIONet.h for definitions
    # irlb  = 0
    # halko = 1
    # feng  = 2
    if svd_solver == 0:
        max_iter = 100 * max_iter

    if issparse(X):
        X = X.tocsc()
        reduced = _an.reduce_kernel(X, dim, max_iter, seed, svd_solver, False)
    else:
        X = np.array(X)
        reduced = _an.reduce_kernel_full(X, dim, max_iter, seed, svd_solver, False)


    # Note S_r.T
    S_r, V, sigma, A, B = (
        reduced["S_r"].T,
        reduced["V"],
        reduced["sigma"],
        reduced["A"],
        reduced["B"],
    )

    if data_is_AnnData:
        adata.obsm[reduction_name] = S_r
        adata.uns[reduction_name] = {}
        adata.uns[reduction_name]["params"] = {"use_highly_variable": use_highly_variable}
        adata.uns[reduction_name]["sigma"] = sigma
        adata.obsm[reduction_name + "_B"] = B

        if use_highly_variable:
            adata.varm[reduction_name + "_V"] = np.zeros(shape=(adata.n_vars, dim))
            adata.varm[reduction_name + "_V"][adata.var["highly_variable"]] = V
            adata.varm[reduction_name + "_A"] = np.zeros(shape=(adata.n_vars, dim))
            adata.varm[reduction_name + "_A"][adata.var["highly_variable"]] = A
        else:
            adata.varm[reduction_name + "_V"] = V
            adata.varm[reduction_name + "_A"] = A

        return adata if copy else None
    else:
        if return_info:
            return (S_r, V.T, sigma, A, B)
        else:
            return S_r
