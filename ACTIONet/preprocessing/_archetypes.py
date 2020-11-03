from typing import Literal, Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from natsort import natsorted
from scipy import sparse

import _ACTIONet as _an

def prune_archetypes(
    adata: AnnData,
    C_trace: list,
    H_trace: list,
    min_specificity_z_threshold: Optional[float] = -3,
    min_cells: Optional[int] = 2,
    copy: Optional[bool] = False
) -> Optional[AnnData]:
    """
    Archetype pruning

    Initial pruning of archetypes

    Parameters
    ----------
    adata:
        Current AnnData object storing the ACTIONet results
    C_trace, H_trace:
        Outputs of pp.run_ACTION()
    min_specificity_z_threshold:
        Controls level of prunning for non-specific archetypes
        (larger values remove more archetypes)
    min_cells:
        Minimum number of influential cells for each archetype
        to be considered nontrivial
    copy
        Determines whether a copy of `adata` is returned.
    Returns
    -------
    adata : anndata.AnnData
        if `copy=True` returns None or else adds fields to `adata`:

        `.obsm['H_stacked']`
        `.obsm['C_stacked']`

        `.uns['obsm_annot']['H_stacked']`
        `.uns['obsm_annot']['C_stacked']`
    """
    pruned = _an.prune_archetypes(C_trace, H_trace, min_specificity_z_threshold, min_cells)

    adata = adata.copy() if copy else adata
    adata.obsm['C_stacked'] = sparse.csc_matrix(pruned['C_stacked'])
    adata.obsm['H_stacked'] = sparse.csc_matrix(pruned['H_stacked'].T)
    adata.uns.setdefault('obsm_annot', {}).update({
        'H_stacked': {'type': np.array([b'internal'], dtype=object)},
        'C_stacked': {'type': np.array([b'internal'], dtype=object)},
    })

    return adata if copy else None

def unify_archetypes(
    adata: AnnData,
    sensitivity: Optional[float] = 1.0,
    normalization_type: Literal[1, 3] = 1,
    edge_threshold: Optional[float] = 0.5,
    copy: Optional[bool] = False
) -> AnnData:
    """\
    Archetype unification

    Aggregates redundant archetypes.

    Parameters
    ----------
    adata
        Current AnnData object storing the ACTIONet results
    sensitivity
    normalization_type
    edge_threshold
    copy
        Determines whether a copy of `adata` is returned.
    Returns
    -------
        adata : anndata.AnnData
        if `copy=True` returns None or else adds fields to `adata`:

        `.obs['assigned_archetype']`
        `.obsm['H_unified']`
        `.obsm['C_unified']`

        `.uns['obsm_annot']['H_unified']`
        `.uns['obsm_annot']['C_unified']`
    """
    # Check for ACTION_S_r and ACTION_H_stacked
    if 'ACTION' not in adata.obsm.keys():
        raise ValueError(
            'Did not find adata.obsm[\'ACTION\']. '
            'Please run pp.reduce_kernel() first.'
        )
    if 'C_stacked' not in adata.obsm.keys() or 'H_stacked' not in adata.obsm.keys():
        raise ValueError(
            'Did not find adata.obsm[\'C_stacked\'] or adata.obsm[\'H_stacked\']. '
            'Please run pp.prune_archetypes() first.'
        )

    adata = adata.copy() if copy else adata
    S_r = adata.obsm['ACTION'].T
    C = adata.obsm['C_stacked']
    if sparse.issparse(C):
        C = C.toarray()
    H = adata.obsm['H_stacked'].T
    if sparse.issparse(H):
        H = H.toarray()
    unified = _an.unify_archetypes(
        S_r,
        C,
        H,
        sensitivity,
        normalization_type,
        edge_threshold,
    )

    adata.obs['assigned_archetype'] = unified['assigned_archetypes']
    adata.obsm['C_unified'] = sparse.csc_matrix(unified['C_unified'])
    adata.obsm['H_unified'] = sparse.csc_matrix(unified['H_unified'].T)
    adata.uns.setdefault('obsm_annot', {}).update({
        'C_unified': {'type': np.array([b'internal'], dtype=object)},
        'H_unified': {'type': np.array([b'internal'], dtype=object)},
    })

    return adata if copy else None
