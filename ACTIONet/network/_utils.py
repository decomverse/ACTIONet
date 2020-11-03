from typing import Optional

import numpy as np
from anndata import AnnData
from scipy import sparse

import _ACTIONet as _an
from ..tools import scale_matrix

def compute_archetype_core_centrality(
    adata: AnnData,
    copy: Optional[bool] = False
) -> AnnData:
    """
    Computes node centrality scores

    Uses graph core-ness to compute node centralities

    Parameters
    ----------
    adata
        AnnData object storing the ACTIONet results
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
        if `copy=True` returns None or else adds fields to `adata`:

        `.obs['node_centrality']`
    """
    if 'ACTIONet' not in adata.obsp.keys():
        raise ValueError(
            'Did not find adata.obsp[\'ACTIONet\']. '
            'Please run nt.build_ACTIONet() first.'
        )
    if 'assigned_archetype' not in adata.obs.keys():
        raise ValueError(
            'Did not find adata.obs[\'assigned_archetype\']. '
            'Please run pp.unify_archetypes() first.'
        )

    adata = adata.copy() if copy else adata
    G = adata.obsp['ACTIONet']
    assignments = adata.obs['assigned_archetype']

    scores = _an.compute_archetype_core_centrality(G, assignments)
    adata.obs['node_centrality'] = scores

    return adata if copy else None

def compute_network_diffusion(
    adata: AnnData,
    alpha: Optional[float] = 0.85,
    n_threads: Optional[int] = 0,
    copy: Optional[bool] = False,
):
    """

    Parameters
    ----------
    adata
        AnnData object storing the ACTIONet results
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
        if `copy=True` returns None or else adds fields to `adata`:

        `.obsm['archetype_footprint']`

        `.uns['obsm_annot']['archetype_footprint']`
    """
    if 'H_unified' not in adata.obsm.keys():
        raise ValueError(
            'Did not find adata.obsm[\'H_unified\']. '
            'Please run pp.unify_archetypes() first.'
        )
    if 'ACTIONet' not in adata.obsp.keys():
        raise ValueError(
            'Did not find adata.obsp[\'ACTIONet\']. '
            'Please run nt.built_network() first.'
        )

    adata = adata.copy() if copy else adata
    H = adata.obsm['H_unified']
    G = adata.obsp['ACTIONet']
    archetype_footprint = _an.compute_network_diffusion(
        G, H, alpha=alpha, thread_no=n_threads
    )
    adata.obsm['archetype_footprint'] = archetype_footprint
    adata.uns.setdefault('obsm_annot', {}).update({
        'archetype_footprint': {'type': np.array([b'generic'], dtype=object)},
    })

    return adata if copy else None

def construct_backbone(
    adata: AnnData,
    scale: Optional[bool] = True,
    network_density: Optional[float] = 1.0,
    mutual_edges_only: Optional[bool] = True,
    layout_compactness: Optional[int] = 50,
    layout_epochs: Optional[int] = 100,
    footprint_alpha: Optional[float] = 0.85,
    n_threads: Optional[int] = 0,
    copy: Optional[bool] = False,
):
    if 'ACTIONet' not in adata.obsp.keys():
        raise ValueError(
            'Did not find adata.obsp[\'ACTIONet\']. '
            'Please run nt.built_network() first.'
        )

    adata = adata.copy() if copy else adata
    if 'archetype_footprint' not in adata.obsm.keys():
        compute_network_diffusion(
            adata,
            alpha=footprint_alpha,
            n_threads=n_threads
        )

    archetype_footprint = adata.obsm['archetype_footprint']
    if scale:
        W = np.exp(scale_matrix(archetype_footprint))
    else:
        W = np.exp(archetype_footprint)

    arch_vis_out = _an.transform_layout(
        sparse.csc_matrix(W),
        coor2D=adata.obsm['X_ACTIONet2D'].T,
        coor3D=adata.obsm['X_ACTIONet3D'].T,
        colRGB=adata.obsm['X_denovo_color'].T,
        n_epochs=layout_epochs,
        compactness_level=layout_compactness,
        thread_no=n_threads,
    )
    arch_G = _an.compute_full_sim(archetype_footprint)
    np.fill_diagonal(arch_G, 0)
    backbone = {
        'G': arch_G,
        'coordinates': arch_vis_out['coordinates'].T,
        'coordinates_3D': arch_vis_out['coordinates_3D'].T,
    }
    adata.uns['backbone'] = backbone

    return adata if copy else None
