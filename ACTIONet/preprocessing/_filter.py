from typing import Optional

import numpy as np
import scanpy as sc
from anndata import AnnData

def filter_adata(
    adata: AnnData,
    min_cells_per_feature: Optional[int] = None,
    min_features_per_cell: Optional[int] = None,
    min_umis_per_cell: Optional[int] = None,
    max_umis_per_cell: Optional[int] = None,
    copy: Optional[bool] = False,
) -> Optional[AnnData]:
    """Filter AnnData by cells or genes.

    Filtering is done iteratively until convergence.

    Parameters
    ----------
    adata
        Annotated data matrix
    min_cells_per_feature
        Minimum number of cells per feature
    min_features_per_cell
        Minimum number of (non-zero) features per cell
    min_umis_per_cell
        Minimum number of UMIs per cell
    max_umis_per_cell
        Maximum number of UMIs per cell
    copy
        Return a copy instead of writing to `adata`

    Returns
    -------
    adata : anndata.AnnData
        if `copy=True` returns None or else filtered `adata`:
    """
    adata = adata.copy() if copy else adata
    original_dims = adata.shape
    previous_dims = None

    while previous_dims != adata.shape:
        previous_dims = adata.shape

        if min_umis_per_cell is not None or max_umis_per_cell is not None or min_features_per_cell is not None:
            sc.pp.filter_cells(
                adata,
                min_counts=min_umis_per_cell,
                max_counts=max_umis_per_cell,
                min_genes=min_features_per_cell,
            )
        if min_cells_per_feature is not None:
            sc.pp.filter_genes(
                adata,
                min_cells=min_cells_per_feature * original_dims[0] if min_cells_per_feature > 0 and min_cells_per_feature < 1 else min_cells_per_feature
            )

    return adata if copy else None
