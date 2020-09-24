from typing import Optional, Tuple, Union

import numpy as np
from anndata import AnnData

import _ACTIONet as _an

def ACTION(
    data: Union[AnnData, np.ndarray],
    k_min: Optional[int] = 2,
    k_max: Optional[int] = 30,
    n_threads: Optional[int] = 0,
    max_it: Optional[int] = 50,
    min_delta: Optional[float] = 1e-300
) -> Tuple[np.ndarray, np.ndarray]:
    """\
    Run ACTION decomposition [Mohammadi2018]_.

    Computes reduced ACTION decomposition.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    k_max
        Max. # of archetypes to consider
    n_threads
        Number of threads. Defaults to number of threads available.
    max_it, min_delta:
        Define stopping conditions of inner AA loop

    Returns
    -------
    C, H
        A dictionary with trace of C and H matrices
    """
    data_is_AnnData = isinstance(data, AnnData)
    if data_is_AnnData:
        if 'ACTION_S_r' not in data.obsm.keys():
            raise ValueError(
                'Did not find adata.obsm[\'ACTION_S_r\']. '
                'Please run pp.reduce_kernel() first.'
            )
        X = data.obsm['ACTION_S_r'].T
    else:
        X = data.T

    result = _an.run_ACTION(X, k_min, k_max, n_threads, max_it, min_delta)
    
    return (result['C'], result['H'])