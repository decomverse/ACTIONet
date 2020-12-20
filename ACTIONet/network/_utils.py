from typing import Optional

import numpy as np
from anndata import AnnData

import _ACTIONet as _an


def compute_archetype_core_centrality(
    adata: AnnData, key: Optional[str] = "ACTIONet", copy: Optional[bool] = False
) -> AnnData:
    """
    Computes node centrality scores

    Uses graph core-ness to compute node centralities

    Parameters
    ----------
    adata
        AnnData object storing the ACTIONet results
    key
        `adata.obsp` key that stores the ACTIONet connectivities
    copy
        Return a copy instead of writing to adata.

    Returns
    -------
        None, if copy is False, ACE: AnnData, if copy is True.
        "node_centrality" is to ACE.obs.
    """
    if "ACTIONet" not in adata.obsp.keys():
        raise ValueError(
            "Did not find adata.obsp['ACTIONet']. "
            "Please run nt.build_network() first."
        )
    if "ACTION" not in adata.obs.keys():
        raise ValueError(
            "Did not find adata.obs['ACTION']. "
            "Please run pp.unify_archetypes() first."
        )

    adata = adata.copy() if copy else adata
    G = adata.obsp["ACTIONet"]
    assignments = adata.obs["ACTION"]

    scores = _an.compute_archetype_core_centrality(G, assignments)
    adata.obs["ACTIONet_centrality"] = scores

    return adata if copy else None


def compute_network_diffusion(
    adata: AnnData,
    archetypes_key: Optional[str] = "ACTION_H_unified",
    alpha: Optional[float] = 0.85,
    n_threads: Optional[int] = 0,
    copy: Optional[bool] = False,
):
    if archetypes_key not in adata.obsm.keys():
        raise ValueError(f"Did not find adata.obsm['{archetypes_key}'].")
    if "ACTIONet" not in adata.obsp.keys():
        raise ValueError(
            "Did not find adata.obsp['ACTIONet']. "
            "Please run nt.built_network() first."
        )

    adata = adata.copy() if copy else adata
    H = adata.obsm[archetypes_key].T
    G = adata.obsp["ACTIONet"]
    archetype_footprint = _an.compute_network_diffusion(
        G, H, alpha=alpha, thread_no=n_threads
    )
    adata.obsm["ACTION_archetype_footprint"] = archetype_footprint

    return adata if copy else None


def construct_backbone(
    adata: AnnData,
    archetypes_key: Optional[str] = "ACTION_H_unified",
    footprint_key: Optional[str] = "ACTION_archetype_footprint",
    scale: Optional[bool] = True,
    network_density: Optional[float] = 1.0,
    mutual_edges_only: Optional[bool] = True,
    layout_compactness: Optional[int] = 50,
    layout_epochs: Optional[int] = 100,
    footprint_alpha: Optional[float] = 0.85,
    n_threads: Optional[int] = 0,
    copy: Optional[bool] = False,
):
    if "ACTIONet" not in adata.obsp.keys():
        raise ValueError(
            "Did not find adata.obsp['ACTIONet']. "
            "Please run nt.built_network() first."
        )

    adata = adata.copy() if copy else adata
    if footprint_key not in adata.obsm.keys():
        compute_network_diffusion(
            adata,
            archetypes_key=archetypes_key,
            alpha=footprint_alpha,
            n_threads=n_threads,
        )

    archetype_footprint = adata.obsm[footprint_key]
    if scale:
        W = np.exp(
            (archetype_footprint - np.mean(archetype_footprint, axis=0))
            / np.std(archetype_footprint, axis=0, ddof=1)
        )
    else:
        W = np.exp(archetype_footprint)

    arch_vis_out = _an.transform_layout(
        W,
        coor2D=adata.obsm["X_ACTIONet_2D"].T,
        coor3D=adata.obsm["X_ACTIONet_3D"].T,
        colRGB=adata.uns["ACTIONet"]["colors"],
        n_epochs=layout_epochs,
        compactness_level=layout_compactness,
        thread_no=n_threads,
    )
    arch_G = _an.compute_full_sim(archetype_footprint)
    backbone = {
        "G": arch_G,
        "coordinates": arch_vis_out["coordinates"],
        "coordinates_3D": arch_vis_out["coordinates_3D"],
    }
    adata.uns["ACTIONet_backbone"] = backbone

    return adata if copy else None
