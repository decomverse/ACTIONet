from typing import Optional, Tuple, Union

from anndata import AnnData

from . import network as nt
from . import preprocessing as pp
from . import tools as tl

def run_ACTIONet(
    adata: AnnData,
    reduction_key: Optional[str] = 'ACTION',
    k_max: Optional[int] = 30,
    min_cells_per_archetype: Optional[int] = 2,
    min_specificity_z_threshold: Optional[int] = -3,
    network_density: Optional[int] = 1,
    mutual_edges_only: Optional[bool] = True,
    layout_compactness: Optional[int] = 50,
    layout_epochs: Optional[int] = 500,
    layout_in_parallel: Optional[bool] = True,
    unification_threshold: Optional[float] = 0.5,
    unification_magnitude: Optional[int] = 3,
    unification_z_threshold: Optional[float] = -3.0,
    footprint_alpha: Optional[float] = 0.85,
    max_iter_ACTION: Optional[int] = 50,
    n_threads: Optional[int] = 0,
    copy: Optional[bool] = False,
) -> Optional[AnnData]:
    """A wrapper function to call all main functions of the ACTIONet

    Parameters
    ----------
    adata
    reduction_key
    k_max
    min_cells_per_archetype
    min_specificity_z_threshold
    network_density
    mutual_edges_only
    layout_compactness
    layout_epochs
    layout_in_parallel
    unification_sensitivity
    footprint_alpha
    max_iter_ACTION
    n_threads
    copy

    Returns
    -------
    adata : anndata.AnnData
        if `copy=True` returns None or else adds fields to `adata`:
    """
    adata = adata.copy() if copy else adata

    # Run ACTION
    C, H = pp.run_ACTION(
        adata,
        reduction_key=reduction_key,
        k_min=2,
        k_max=k_max,
        n_threads=n_threads,
        max_it=max_iter_ACTION,
        min_delta=1e-300
    )

    # Prune nonspecific and/or unreliable archetypes
    pp.prune_archetypes(
        adata,
        C,
        H,
        min_specificity_z_threshold=min_specificity_z_threshold,
        min_cells=min_cells_per_archetype
    )

    # Build ACTIONet
    nt.build_ACTIONet(
        adata,
        density=network_density,
        n_threads=n_threads,
        mutual_edges_only=mutual_edges_only
    )

    # Layout ACTIONet
    if not layout_in_parallel:
        nt.layout_ACTIONet(
            adata,
            reduction_key=reduction_key,
            compactness_level=layout_compactness,
            n_epochs=layout_epochs,
            n_threads=1,
        )
    else:
        nt.layout_ACTIONet(
            adata,
            reduction_key=reduction_key,
            compactness_level=layout_compactness,
            n_epochs=layout_epochs,
            n_threads=n_threads,
        )

    # Identiy equivalent classes of archetypes and group them together
    pp.unify_archetypes(
        adata,
        z_threshold=unification_z_threshold,
        cor_threshold=unification_threshold,
        magnification=unification_magnitude,
    )

    # Use graph core of global and induced subgraphs to infer centrality/quality of
    # each cell
    nt.compute_archetype_core_centrality(adata)

    # Smooth archetype footprints
    nt.compute_network_diffusion(adata, alpha=footprint_alpha, n_threads=n_threads)

    # Compute gene specificity for each archetype
    nt.compute_archetype_feature_specificity(adata)

    nt.construct_backbone(
        adata,
        network_density=network_density,
        mutual_edges_only=mutual_edges_only,
        layout_compactness=layout_compactness,
        layout_epochs=int(layout_epochs/5),
        n_threads=1,
    )

    return adata if copy else None

def reconstruct_ACTIONet(
    adata: AnnData,
    reduction_key: Optional[str] = 'ACTION',
    network_density: Optional[int] = 1,
    mutual_edges_only: Optional[bool] = True,
    layout_compactness: Optional[int] = 50,
    layout_epochs: Optional[int] = 500,
    layout_in_parallel: Optional[bool] = False,
    n_threads: Optional[int] = 0,
    copy: Optional[bool] = False,
) -> Optional[AnnData]:
    """Reconstructs the ACTIONet graph with the new parameters (uses prior decomposition)

    Parameters
    ----------
    adata
    reduction_key
    network_density
    mutual_edges_only
    layout_compactness
    layout_epochs
    layout_in_parallel
    n_threads
    copy

    Returns
    -------
    adata : anndata.AnnData
        if `copy=True` returns None or else adds fields to `adata`:
    """
    adata = adata.copy() if copy else adata

    # re-build ACTIONet
    nt.build_ACTIONet(
        adata,
        density=network_density,
        n_threads=n_threads,
        mutual_edges_only=mutual_edges_only
    )

    # layout ACTIONet
    if not layout_in_parallel:
        nt.layout_ACTIONet(
            adata,
            reduction_key=reduction_key,
            compactness_level=layout_compactness,
            n_epochs=layout_epochs,
            n_threads=1,
        )
    else:
        nt.layout_ACTIONet(
            adata,
            reduction_key=reduction_key,
            compactness_level=layout_compactness,
            n_epochs=layout_epochs,
            n_threads=n_threads,
        )

    nt.construct_backbone(
        adata,
        network_density=network_density,
        mutual_edges_only=mutual_edges_only,
        layout_compactness=layout_compactness,
        layout_epochs=int(layout_epochs/5),
        n_threads=1,
    )

    return adata if copy else None

def rerun_layout(
    adata: AnnData,
    reduction_key: Optional[str] = 'ACTION',
    network_density: Optional[int] = 1,
    mutual_edges_only: Optional[bool] = True,
    layout_compactness: Optional[int] = 50,
    layout_epochs: Optional[int] = 1000,
    n_threads: Optional[int] = 1,
    copy: Optional[bool] = False,
) -> Optional[AnnData]:
    """Rerun layout on the ACTIONet graph with new parameters

    Parameters
    ----------
    adata
    network_density
    mutual_edges_only
    layout_compactness
    layout_epochs
    n_threads
    copy

    Returns
    -------
    adata : anndata.AnnData
        if `copy=True` returns None or else adds fields to `adata`:
    """
    adata = adata.copy() if copy else adata

    # re-layout ACTIONet
    nt.layout_ACTIONet(
        adata,
        reduction_key=reduction_key,
        compactness_level=layout_compactness,
        n_epochs=layout_epochs,
        n_threads=n_threads,
    )
    nt.construct_backbone(
        adata,
        network_density=network_density,
        mutual_edges_only=mutual_edges_only,
        layout_compactness=layout_compactness,
        layout_epochs=int(layout_epochs/5),
        n_threads=1,
    )

    return adata if copy else None

def rerun_archetype_aggregation(
    adata: AnnData,
    network_density: Optional[int] = 1,
    mutual_edges_only: Optional[bool] = True,
    layout_compactness: Optional[int] = 50,
    layout_epochs: Optional[int] = 500,
    unification_threshold: Optional[float] = 0.5,
    unification_magnitude: Optional[int] = 3,
    unification_z_threshold: Optional[float] = -3.0,
    footprint_alpha: Optional[float] = 0.85,
    n_threads: Optional[int] = 0,
    copy: Optional[bool] = False,
) -> Optional[AnnData]:
    adata = adata.copy() if copy else adata

    pp.unify_archetypes(
        adata,
        z_threshold=unification_z_threshold,
        cor_threshold=unification_threshold,
        magnification=unification_magnitude,
    )
    nt.compute_archetype_core_centrality(adata)
    nt.compute_network_diffusion(adata, alpha=footprint_alpha, n_threads=n_threads)
    nt.compute_archetype_feature_specificity(adata)
    nt.construct_backbone(
        adata,
        network_density=network_density,
        mutual_edges_only=mutual_edges_only,
        layout_compactness=layout_compactness,
        layout_epochs=int(layout_epochs/5),
        n_threads=1,
    )

    return adata if copy else None
