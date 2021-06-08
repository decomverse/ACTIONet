from typing import Optional

import numpy as np
from scipy import sparse
from anndata import AnnData
from .. import misc_utils as ut


def normalize_adata(
        adata: AnnData,
        layer_name: Optional[str] = None,
        layer_name_out: Optional[str] = None,
        log_scale: Optional[bool] = True,
        copy: Optional[bool] = False
) -> Optional[AnnData]:

    adata = adata.copy() if copy else adata

    if layer_name is not None:
        S = ut.rescale_matrix(adata.layers[layer_name], log_scale=log_scale)
    else:
        S = ut.rescale_matrix(adata.X, log_scale=log_scale)

    S = sparse.csc_matrix(S)

    if layer_name_out is not None:
        adata.layers[layer_name_out] = S
    else:
        adata.X = S

    return adata if copy else None
