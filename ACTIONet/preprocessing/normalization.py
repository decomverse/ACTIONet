from typing import Optional

import numpy as np
from scipy import sparse
from anndata import AnnData
from ..tools import misc_utils as ut


def normalize(
    adata: AnnData,
    layer_key: Optional[str] = None,
    layer_key_out: Optional[str] = None,
    log_scale: Optional[bool] = True,
    copy: Optional[bool] = False,
) -> Optional[AnnData]:

    adata = adata.copy() if copy else adata

    if layer_key is not None:
        S = ut.rescale_matrix(adata.layers[layer_key], log_scale=log_scale)
    else:
        S = ut.rescale_matrix(adata.X, log_scale=log_scale)

    S = sparse.csc_matrix(S)

    if layer_key_out is not None:
        adata.layers[layer_key_out] = S
    else:
        adata.X = S

    return adata if copy else None
