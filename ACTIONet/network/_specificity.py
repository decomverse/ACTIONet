from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import issparse

import _ACTIONet as _an

def _compute_archetype_specificity(S, H):
    if issparse(S):
        return _an.compute_archetype_feature_specificity(S.tocsc(), H)
    return _an.compute_archetype_feature_specificity_full(S, H)

def _compute_cluster_specificity(S, assignments):
    if issparse(S):
        return _an.compute_cluster_feature_specificity(S.tocsc(), assignments)
    return _an.compute_cluster_feature_specificity_full(S, assignments)

def compute_archetype_feature_specificity(
    adata: AnnData,
    copy: Optional[bool] = False
) -> AnnData:
    """\
    Computes Feature (i.e., gene) specificity of archetypes

    Uses Archetype footprints to estimate markers (soft clustering)
    Parameters
    ----------
    adata
        Current AnnData object storing the ACTIONet results
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
        adata : anndata.AnnData
        if `copy=True` returns None or else adds fields to `adata`:

        `.varm['unified_feature_profile']`
        `.varm['unified_feature_specificity']`

        `.uns['varm_annot']['unified_feature_profile']`
        `.uns['varm_annot']['unified_feature_specificity']`
    """
    if 'archetype_footprint' not in adata.obsm.keys():
        raise ValueError(
            'Did not find adata.obsm[\'archetype_footprint\']. '
            'Please run nt.compute_network_diffusion() first.'
        )

    adata = adata.copy() if copy else adata
    S = adata.X.T
    H = adata.obsm['archetype_footprint'].T

    specificity = _compute_archetype_specificity(S, H)
    adata.varm['unified_feature_profile'] = specificity['archetypes']
    adata.varm['unified_feature_specificity'] = specificity['upper_significance']
    adata.uns.setdefault('varm_annot', {}).update({
        'unified_feature_profile': {'type': np.array([b'internal'], dtype=object)},
        'unified_feature_specificity': {'type': np.array([b'reduction'], dtype=object)},
    })

    return adata if copy else None


def compute_cluster_feature_specificity(
    adata: AnnData,
    cluster_key: Optional[str] = 'leiden',
    copy: Optional[bool] = False
) -> Optional[AnnData]:
    """\
    Computes Feature (i.e., gene) specificity of clusters

    Uses cluster membership vector to estimate markers (disjoint clustering)
    Parameters
    ----------
    adata
        Current AnnData object storing the ACTIONet results
    cluster_key
        Key in `adata.obs` that holds the clustering variable
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
        adata : anndata.AnnData
        if `copy=True` returns None or else adds fields to `adata`:

        `.varm[f'{cluster_key}_profile']`
        `.varm[f'{cluster_key}_feature_specificity']`
    """
    if cluster_key not in adata.obs.keys():
        raise ValueError(f'Did not find adata.obs[\'{cluster_key}\'].')
    adata = adata.copy() if copy else adata

    S = adata.X.T
    assignments = adata.obs[cluster_key]
    if isinstance(assignments, pd.Series):
        assignments = pd.factorize(assignments)[0]

    specificity = _compute_cluster_specificity(S, assignments)

    adata.varm[f'{cluster_key}_profile'] = specificity['archetypes']
    adata.varm[f'{cluster_key}_upper_significance'] = specificity['upper_significance']

    return adata if copy else None
