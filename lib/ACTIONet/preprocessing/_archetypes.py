from typing import Optional, Union

import numpy as np
import pandas as pd
from anndata import AnnData
from natsort import natsorted
from scipy.sparse import issparse, spmatrix

import _ACTIONet as _an

def prune_archetypes(
    adata: AnnData,
    C_trace: list,
    H_trace: list,
    min_specificity_z_threshold: Optional[float] = -3,
    min_cells_per_archetype: Optional[int] = 2,    
    copy: Optional[bool] = False    
) -> Optional[AnnData]:
    """\
    Archetype pruning

    Initial pruning of archetypes

    Parameters
    ----------
    adata:
        Current AnnData object storing the ACTIONet results    
    C_trace, H_trace:
        Output of pp.ACTION()
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

        `.obsm['ACTION_C_stacked']`
        `.obsm['ACTION_H_stacked']`
        `.uns['ACTION']['pruned']`
    """
    if 'ACTION' not in adata.uns.keys():
        raise ValueError(
            'Did not find adata.uns[\'ACTION\']. '
            'Please run pp.ACTION() first.'
        )

    pruned = _an.prune_archetypes(C_trace, H_trace, min_specificity_z_threshold)
    
    adata = adata.copy() if copy else adata
    adata.obsm['ACTION_C_stacked'] = pruned['C_stacked']
    adata.obsm['ACTION_H_stacked'] = pruned['H_stacked'].T
    adata.uns['ACTION'].setdefault('archetypes', {}).update({
        'pruned': pruned['selected_archs']
    })
    
    return adata if copy else None

def unify_archetypes(
    adata: AnnData,
    min_edge_weight: Optional[float] = 0.5,
    min_coreness: Optional[int] = 0,
    resolution: Optional[float] = 1.0,
    min_repeat: Optional[int] = 2,
    n_threads: Optional[int] = 0,
    alpha: Optional[float] = 0.05,
    beta: Optional[float] = 0.5,
    outlier_threshold: Optional[float] = 0.5,
    min_points: Optional[int] = 5,
    min_cluster_size: Optional[int] = 5,
    cond_threshold: Optional[float] = 10.0,
    normalization_type: Optional[int] = 3,
    preprocess_adj: Optional[bool] = False,
    reduce_G: Optional[bool] = False,
    method_type: Optional[int] = 0,
    sensitivity: Optional[float] = 0.95,
    copy: Optional[bool] = False
) -> AnnData:
    """\
    Archetype unification

    Aggregates redundant archetypes.

    Parameters
    ----------
    adata:
        Current AnnData object storing the ACTIONet results

    copy
        Determines whether a copy of `adata` is returned. 
    Returns
    -------
        adata : anndata.AnnData
        if `copy=True` returns None or else adds fields to `adata`:

        `.obsm['ACTION_C_unified']`
        `.obsm['ACTION_H_unified']`
        `.uns['ACTION']['unified']`
    """
    # Check for ACTION_S_r and ACTION_H_stacked
    if 'ACTION_S_r' not in adata.obsm.keys():
        raise ValueError(
            'Did not find adata.obsm[\'ACTION_S_r\']. '
            'Please run pp.reduce_kernel() first.'
        )
    if 'ACTION_C_stacked' not in adata.obsm.keys() or 'ACTION_H_stacked' not in adata.obsm.keys():
        raise ValueError(
            'Did not find adata.obsm[\'ACTION_C_stacked\'] or adata.obsm[\'ACTION_H_stacked\']. '
            'Please run pp.prune_archetypes() first.'
        )
    
    adata = adata.copy() if copy else adata
    unified = _an.unify_archetypes(
        adata.obsm['ACTION_S_r'].T,
        adata.obsm['ACTION_C_stacked'],
        adata.obsm['ACTION_H_stacked'].T,
        min_edge_weight,
        min_coreness,
        resolution,
        min_repeat,
        n_threads,
        alpha,
        beta,
        outlier_threshold,
        min_points,
        min_cluster_size,
        cond_threshold,
        normalization_type,
        preprocess_adj,
        reduce_G,
        method_type,
        sensitivity,
    )
    
    adata.obsm['ACTION_C_unified'] = unified['C_unified']
    adata.obsm['ACTION_H_unified'] = unified['H_unified'].T
    adata.uns['ACTION'].setdefault('archetypes', {}).update({
        'unified': unified['selected_archetypes']
    })

    groups = unified['assigned_archetypes']
    adata.obs['ACTION'] = pd.Categorical(
        values=groups.astype('U'),
        categories=natsorted(map(str, np.unique(groups))),
    )

    return adata if copy else None
