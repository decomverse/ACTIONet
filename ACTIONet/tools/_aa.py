from typing import Optional, Tuple

import numpy as np

import _ACTIONet as _an


def aa(
    A: np.ndarray,
    k: int,
    max_iter: Optional[int] = 50,
    min_delta: Optional[float] = 1e-16,
) -> Tuple[np.ndarray, np.ndarray]:
    """\
    Archetypal Analysis (AA)

    Runs SPA algorithm to solve separable NMF problem.

    Parameters
    ----------
    A:
        Input matrix
    W0:
        Initial estimate of archetypes with k (# archetype) columns
    max_it, min_delta:
        Define stopping conditions

    Returns
    -------
    C:
        Convex matrix of archetype coefficients (#observations x # archetypes)
    W:
        Matrix of archetypes (# features x # archetypes)
    H:
        Convex matrix of observation coefficients (# archetypes x # observations)    
        
    """
    result = _an.run_AA(A, k)
    return (
        AA_out["C"],
        AA_out["H"],
    )
