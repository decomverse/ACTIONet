from typing import Optional, Tuple, Union

import numpy as np
from anndata import AnnData

import _ACTIONet as _an


def ACTION(
    data: Union[AnnData, np.ndarray],
    reduction_name: Optional[str] = "ACTION",
    k_min: Optional[int] = 2,
    k_max: Optional[int] = 30,
    thread_no: Optional[int] = 0,
    max_it: Optional[int] = 50,
    min_delta: Optional[float] = 1e-300,
) -> dict:
    """\
    Run ACTION decomposition [Mohammadi, 2018]

    Computes reduced ACTION decomposition.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    reduction_name:
        Key of '.obms' to use as input for ACTION (default="ACTION")
        Ignored if 'data' is not an 'AnnData' object.
    k_min, k_max
        Min. and max. # of archetypes to consider.
    thread_no
        Number of threads. Defaults to number of threads available - 2.
    max_it, min_delta:
        Define stopping conditions of inner AA loop

    Returns
    -------
    ACTION_out
        A dictionary with traces of C and H matrices
    """
    data_is_AnnData = isinstance(data, AnnData)

    if data_is_AnnData:
        if reduction_name not in data.obsm.keys():
            raise ValueError("Did not find adata.obsm['" + reduction_name + "'].")
        else:
            X = data.obsm[reduction_name]
    else:
        X = data

    X = X.T.astype(dtype=np.float64)
    ACTION_out = _an.run_ACTION(X, k_min, k_max, thread_no, max_it, min_delta)

    return ACTION_out
