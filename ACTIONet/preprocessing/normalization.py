from typing import Optional, Union

from anndata import AnnData
import numpy as np

from .. import tools as tl


def normalize(
    adata: AnnData,
    layer_key: Optional[str] = None,
    layer_key_out: Optional[str] = None,
    log_transform: Optional[bool] = True,
    scale_factor: Union[str, float, int, None] = 1e3,
    copy: Optional[bool] = False,
) -> Optional[AnnData]:
    adata = adata.copy() if copy else adata

    if (
        "norm_method" in adata.uns["metadata"].keys()
    ):  # Already normalized? leave it alone!
        return adata if copy else None

    adata = adata.copy() if copy else adata

    if layer_key is None and "default_assay" in adata.uns["metadata"].keys():
        layer_key = adata.uns["metadata"]["default_assay"]

    if layer_key is not None:
        if layer_key not in adata.layers.keys():
            raise ValueError("Did not find adata.layers['" + layer_key + "']. ")
        S = adata.layers[layer_key]
    else:
        S = adata.X

    S = tl.normalize_matrix(S, log_transform=log_transform, scale_factor=scale_factor)

    adata.uns["metadata"]["norm_method"] = b"default"

    if layer_key_out is not None:
        adata.uns["metadata"]["default_assay"] = layer_key_out
        adata.layers[layer_key_out] = S
    else:
        adata.uns["metadata"]["default_assay"] = None
        adata.X = S

    return adata if copy else None
