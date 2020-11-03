from typing import Optional, Tuple, Union

import numpy as np
from anndata import AnnData

import _ACTIONet as _an

def run_ACTION(
    adata: AnnData,
    reduction_key: Optional[str] = 'ACTION',
    k_min: Optional[int] = 2,
    k_max: Optional[int] = 30,
    n_threads: Optional[int] = 0,
    max_it: Optional[int] = 50,
    min_delta: Optional[float] = 1e-300
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Run ACTION decomposition.

    Computes reduced ACTION decomposition.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    reduction_key
        Key in `adata.obsm` that contains reduced cell representations
    k_min
        Min. # of archetypes to consider
    k_max
        Max. # of archetypes to consider
    n_threads
        Number of threads. Defaults to number of threads available.
    max_it, min_delta:
        Define stopping conditions of inner AA loop

    Returns
    -------
    C, H
        A tuple with trace of C and H matrices
    """
    if reduction_key not in adata.obsm.keys():
        raise ValueError(
            f'Did not find adata.obsm[\'{reduction_key}\']. '
            'Please run pp.reduce_kernel() first.'
        )
    X = adata.obsm[reduction_key].T

    result = _an.run_ACTION(X, k_min, k_max, n_threads, max_it, min_delta)

    return (result['C'], result['H'])
