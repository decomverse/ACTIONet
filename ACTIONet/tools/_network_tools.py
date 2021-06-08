from typing import Optional, Union, Tuple
from typing_extensions import Literal

import numpy as np
from anndata import AnnData

import _ACTIONet as _an


def run_LPA(
        G: np.ndarray,
        labels: np.ndarray,
        lambda_val: Optional[float] = 1,
        max_iter: Optional[int] = 3,
        sig_threshold: Optional[int] = 3,
        fixed_labels: Optional[Union[np.ndarray, list]] = np.empty(0)
) -> [np.ndarray]:

    result = _an.run_LPA(G, labels, lambda_val, max_iter, sig_threshold, fixed_labels)
    return result
